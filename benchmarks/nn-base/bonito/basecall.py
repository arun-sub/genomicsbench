"""
Bonito Basecaller
"""

import os, sys, time
from glob import glob
from warnings import warn
from logging import getLogger
from os.path import realpath, splitext, dirname
from multiprocessing import Process, Queue, Lock, cpu_count
from datetime import timedelta
from collections import OrderedDict
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from ont_fast5_api.fast5_interface import get_fast5_file
from bonito_cuda_runtime import CuModel

import numpy as np
from tqdm import tqdm
import torch
import toml
from scipy.signal import find_peaks
import torch.nn as nn
from torch import sigmoid
from torch.jit import script
from torch.autograd import Function
from torch.nn import ReLU, LeakyReLU
from torch.nn import Module, ModuleList, Sequential, Conv1d, BatchNorm1d, Dropout
from fast_ctc_decode import beam_search, viterbi_search


################################################################################
# model
################################################################################
@script
def swish_jit_fwd(x):
    return x * sigmoid(x)


@script
def swish_jit_bwd(x, grad):
    x_s = sigmoid(x)
    return grad * (x_s * (1 + x * (1 - x_s)))


class SwishAutoFn(Function):

    @staticmethod
    def forward(ctx, x):
        ctx.save_for_backward(x)
        return swish_jit_fwd(x)

    @staticmethod
    def backward(ctx, grad):
        x = ctx.saved_tensors[0]
        return swish_jit_bwd(x, grad)


class Swish(Module):
    """
    Swish Activation function

    https://arxiv.org/abs/1710.05941
    """
    def forward(self, x):
        return SwishAutoFn.apply(x)


activations = {
    "relu": ReLU,
    "swish": Swish,
}


class Model(Module):
    """
    Model template for QuartzNet style architectures

    https://arxiv.org/pdf/1910.10261.pdf
    """
    def __init__(self, config):
        super(Model, self).__init__()
        if 'qscore' not in config:
            self.qbias = 0.0
            self.qscale = 1.0
        else:
            self.qbias = config['qscore']['bias']
            self.qscale = config['qscore']['scale']

        self.config = config
        self.stride = config['block'][0]['stride'][0]
        self.alphabet = config['labels']['labels']
        self.features = config['block'][-1]['filters']
        self.encoder = Encoder(config)
        self.decoder = Decoder(self.features, len(self.alphabet))

    def forward(self, x):
        encoded = self.encoder(x)
        return self.decoder(encoded)

    def decode(self, x, beamsize=5, threshold=1e-3, qscores=False, return_path=False):
        if beamsize == 1 or qscores:
            seq, path  = viterbi_search(x, self.alphabet, qscores, self.qscale, self.qbias)
        else:
            seq, path = beam_search(x, self.alphabet, beamsize, threshold)
        if return_path: return seq, path
        return seq


class Encoder(Module):
    """
    Builds the model encoder
    """
    def __init__(self, config):
        super(Encoder, self).__init__()
        self.config = config

        features = self.config['input']['features']
        activation = activations[self.config['encoder']['activation']]()
        encoder_layers = []

        for layer in self.config['block']:
            encoder_layers.append(
                Block(
                    features, layer['filters'], activation,
                    repeat=layer['repeat'], kernel_size=layer['kernel'],
                    stride=layer['stride'], dilation=layer['dilation'],
                    dropout=layer['dropout'], residual=layer['residual'],
                    separable=layer['separable'],
                )
            )

            features = layer['filters']

        self.encoder = Sequential(*encoder_layers)

    def forward(self, x):
        return self.encoder(x)


class TCSConv1d(Module):
    """
    Time-Channel Separable 1D Convolution
    """
    def __init__(self, in_channels, out_channels, kernel_size, stride=1, padding=0, dilation=1, groups=1, bias=False, separable=False):

        super(TCSConv1d, self).__init__()
        self.separable = separable

        if separable:
            self.depthwise = Conv1d(
                in_channels, in_channels, kernel_size=kernel_size, stride=stride,
                padding=padding, dilation=dilation, bias=bias, groups=in_channels
            )

            self.pointwise = Conv1d(
                in_channels, out_channels, kernel_size=1, stride=stride,
                dilation=dilation, bias=bias, padding=0
            )
        else:
            self.conv = Conv1d(
                in_channels, out_channels, kernel_size=kernel_size,
                stride=stride, padding=padding, dilation=dilation, bias=bias
            )

    def forward(self, x):
        if self.separable:
            x = self.depthwise(x)
            x = self.pointwise(x)
        else:
            x = self.conv(x)
        return x


class Block(Module):
    """
    TCSConv, Batch Normalisation, Activation, Dropout
    """
    def __init__(self, in_channels, out_channels, activation, repeat=5, kernel_size=1, stride=1, dilation=1, dropout=0.0, residual=False, separable=False):

        super(Block, self).__init__()

        self.use_res = residual
        self.conv = ModuleList()

        _in_channels = in_channels
        padding = self.get_padding(kernel_size[0], stride[0], dilation[0])

        # add the first n - 1 convolutions + activation
        for _ in range(repeat - 1):
            self.conv.extend(
                self.get_tcs(
                    _in_channels, out_channels, kernel_size=kernel_size,
                    stride=stride, dilation=dilation,
                    padding=padding, separable=separable
                )
            )

            self.conv.extend(self.get_activation(activation, dropout))
            _in_channels = out_channels

        # add the last conv and batch norm
        self.conv.extend(
            self.get_tcs(
                _in_channels, out_channels,
                kernel_size=kernel_size,
                stride=stride, dilation=dilation,
                padding=padding, separable=separable
            )
        )

        # add the residual connection
        if self.use_res:
            self.residual = Sequential(*self.get_tcs(in_channels, out_channels))

        # add the activation and dropout
        self.activation = Sequential(*self.get_activation(activation, dropout))

    def get_activation(self, activation, dropout):
        return activation, Dropout(p=dropout)

    def get_padding(self, kernel_size, stride, dilation):
        if stride > 1 and dilation > 1:
            raise ValueError("Dilation and stride can not both be greater than 1")
        return (kernel_size // 2) * dilation

    def get_tcs(self, in_channels, out_channels, kernel_size=1, stride=1, dilation=1, padding=0, bias=False, separable=False):
        return [
            TCSConv1d(
                in_channels, out_channels, kernel_size,
                stride=stride, dilation=dilation, padding=padding,
                bias=bias, separable=separable
            ),
            BatchNorm1d(out_channels, eps=1e-3, momentum=0.1)
        ]

    def forward(self, x):
        _x = x
        for layer in self.conv:
            _x = layer(_x)
        if self.use_res:
            _x += self.residual(x)
        return self.activation(_x)


class Decoder(Module):
    """
    Decoder
    """
    def __init__(self, features, classes):
        super(Decoder, self).__init__()
        self.layers = Sequential(Conv1d(features, classes, kernel_size=1, bias=True))

    def forward(self, x):
        x = self.layers(x)
        return nn.functional.log_softmax(x.transpose(1, 2), dim=2)

################################################################################
# util
################################################################################
def load_model(dirname, device, weights=None, half=False, chunksize=0, use_rt=False):
    """
    Load a model from disk
    """
    dirname = os.path.join(os.getcwd(), dirname)

    device = torch.device(device)
    config = os.path.join(dirname, 'config.toml')
    weights = os.path.join(dirname, 'weights_%s.tar' % weights)
    model = Model(toml.load(config))

    state_dict = torch.load(weights, map_location=device)
    new_state_dict = OrderedDict()
    for k, v in state_dict.items():
        name = k.replace('module.', '')
        new_state_dict[name] = v

    model.load_state_dict(new_state_dict)

    if use_rt:
        model = CuModel(model.config, chunksize, new_state_dict)

    if half: model = model.half()
    model.eval()
    model.to(device)
    return model


def half_supported():
    """
    Returns whether FP16 is support on the GPU
    """
    return torch.cuda.get_device_capability()[0] >= 7


def chunk(raw_data, chunksize, overlap):
    """
    Convert a read into overlapping chunks before calling
    """
    if chunksize > 0 and raw_data.shape[0] > chunksize:
        num_chunks = raw_data.shape[0] // (chunksize - overlap) + 1
        tmp = torch.zeros(num_chunks * (chunksize - overlap)).type(raw_data.dtype)
        tmp[:raw_data.shape[0]] = raw_data
        return tmp.unfold(0, chunksize, chunksize - overlap).unsqueeze(1)
    return raw_data.unsqueeze(0).unsqueeze(0)


def stitch(predictions, overlap):
    """
    Stitch predictions together with a given overlap
    """
    if predictions.shape[0] == 1:
        return predictions.squeeze(0)
    stitched = [predictions[0, 0:-overlap]]
    for i in range(1, predictions.shape[0] - 1): stitched.append(predictions[i][overlap:-overlap])
    stitched.append(predictions[-1][overlap:])
    return np.concatenate(stitched)


def mean_qscore_from_qstring(qstring):
    """
    Convert qstring into a mean qscore
    """
    if len(qstring) == 0: return 0.0
    err_probs = [10**((ord(c) - 33) / -10) for c in qstring]
    mean_err = np.mean(err_probs)
    return -10 * np.log10(max(mean_err, 1e-4))


def get_raw_data(filename):
    """
    Get the raw signal and read id from the fast5 files
    """
    with get_fast5_file(filename, 'r') as f5_fh:
        for res in f5_fh.get_reads():
            yield Read(res, filename)


class Read:

    def __init__(self, read, filename):

        self.read_id = read.read_id
        self.run_id = read.get_run_id().decode()
        self.filename = os.path.basename(read.filename)

        read_attrs = read.handle[read.raw_dataset_group_name].attrs
        channel_info = read.handle[read.global_key + 'channel_id'].attrs

        self.offset = int(channel_info['offset'])
        self.sampling_rate = channel_info['sampling_rate']
        self.scaling = channel_info['range'] / channel_info['digitisation']

        self.mux = read_attrs['start_mux']
        self.channel = channel_info['channel_number'].decode()
        self.start = read_attrs['start_time'] / self.sampling_rate
        self.duration = read_attrs['duration'] / self.sampling_rate

        # no trimming
        self.template_start = self.start
        self.template_duration = self.duration

        raw = read.handle[read.raw_dataset_name][:]
        scaled = np.array(self.scaling * (raw + self.offset), dtype=np.float32)
        self.signal = norm_by_noisiest_section(scaled)


def norm_by_noisiest_section(signal, samples=100, threshold=6.0):
    """
    Normalise using the medmad from the longest continuous region where the
    noise is above some threshold relative to the std of the full signal.
    """
    threshold = signal.std() / threshold
    noise = np.ones(signal.shape)

    for idx in np.arange(signal.shape[0] // samples):
        window = slice(idx * samples, (idx + 1) * samples)
        noise[window] = np.where(signal[window].std() > threshold, 1, 0)

    # start and end low for peak finding
    noise[0] = 0; noise[-1] = 0
    peaks, info = find_peaks(noise, width=(None, None))

    if len(peaks):
        widest = np.argmax(info['widths'])
        med, mad = med_mad(signal[info['left_bases'][widest]: info['right_bases'][widest]])
    else:
        med, mad = med_mad(signal)
    return (signal - med) / mad


def med_mad(x, factor=1.4826):
    """
    Calculate signal median and median absolute deviation
    """
    med = np.median(x)
    mad = np.median(np.absolute(x - med)) * factor
    return med, mad




################################################################################
# io
################################################################################
logger = getLogger('bonito')

def summary_file():
    """
    Return the filename to use for the summary tsv.
    """
    if sys.stdout.isatty():
        return 'summary.tsv'
    return '%s_summary.tsv' % splitext(realpath('/dev/fd/1'))[0]


def write_summary_header(fd, sep='\t'):
    """
    Write the summary tsv header.
    """
    fields = [
        'filename',
        'read_id',
        'run_id',
        'channel',
        'mux',
        'start_time',
        'duration',
        'template_start',
        'template_duration',
        'sequence_length_template',
        'mean_qscore_template',
    ]
    fd.write('%s\n' % sep.join(fields))
    fd.flush()


def write_summary_row(fd, read, seqlen, qscore, sep='\t'):
    """
    Write a summary tsv row.
    """
    fields = [str(field) for field in [
        read.filename,
        read.read_id,
        read.run_id,
        read.channel,
        read.mux,
        read.start,
        read.duration,
        read.template_start,
        read.template_duration,
        seqlen,
        qscore,
    ]]

    fd.write('%s\n' % sep.join(fields))
    fd.flush()

def write_fasta(header, sequence, fd=sys.stdout):
    """
    Write a fasta record to a file descriptor.
    """
    fd.write(">%s\n" % header)
    fd.write("%s\n" % sequence)
    fd.flush()


def write_fastq(header, sequence, qstring, fd=sys.stdout):
    """
    Write a fastq record to a file descriptor.
    """
    fd.write("@%s\n" % header)
    fd.write("%s\n" % sequence)
    fd.write("+\n")
    fd.write("%s\n" % qstring)
    fd.flush()


class PreprocessReader(Process):
    """
    Reader Processor that reads and processes fast5 files
    """
    def __init__(self, directory, maxsize=5):
        super().__init__()
        self.directory = directory
        self.queue = Queue(maxsize)

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.stop()

    def run(self):
        for fast5 in tqdm(glob("%s/*fast5" % self.directory), ascii=True, ncols=100, leave=False):
            for read in get_raw_data(fast5):
                self.queue.put(read)
        self.queue.put(None)

    def stop(self):
        self.join()


class DecoderWriterPool:
   """
   Simple pool of decoder writers
   """
   def __init__(self, model, procs=4, **kwargs):
       self.lock = Lock()
       self.queue = Queue()
       self.procs = procs if procs else cpu_count()
       self.decoders = []

       with open(summary_file(), 'w') as summary:
           write_summary_header(summary)

       for _ in range(self.procs):
           decoder = DecoderWriter(model, self.queue, self.lock, **kwargs)
           decoder.start()
           self.decoders.append(decoder)

   def stop(self):
       for decoder in self.decoders: self.queue.put(None)
       for decoder in self.decoders: decoder.join()

   def __enter__(self):
       return self

   def __exit__(self, exc_type, exc_val, exc_tb):
       self.stop()


class DecoderWriter(Process):
    """
    Decoder Process that writes output records to stdout
    """
    def __init__(self, model, queue, lock, fastq=False, beamsize=5):
        super().__init__()
        self.queue = queue
        self.lock = lock
        self.model = model
        self.fastq = fastq
        self.beamsize = beamsize

    def run(self):
        while True:
            job = self.queue.get()
            if job is None: return
            read, predictions = job

            # convert logprobs to probs
            predictions = np.exp(predictions.astype(np.float32))

            sequence, path = self.model.decode(
                predictions, beamsize=self.beamsize, qscores=True, return_path=True
            )
            sequence, qstring = sequence[:len(path)], sequence[len(path):]
            mean_qscore = mean_qscore_from_qstring(qstring)

            if not self.fastq:  # beam search
                qstring = '*'
                sequence, path = self.model.decode(
                    predictions, beamsize=self.beamsize, qscores=False, return_path=True
                )

            if sequence:
                with self.lock, open(summary_file(), 'a') as summary:
                    if self.fastq:
                        write_fastq(read.read_id, sequence, qstring)
                    else:
                        write_fasta(read.read_id, sequence)
                    write_summary_row(summary, read, len(sequence), mean_qscore)
            else:
                logger.warn("> skipping empty sequence %s", read.read_id)

################################################################################
# basecalling 
################################################################################

def main(args):

    sys.stderr.write("> loading model\n")

    model = load_model(
        args.model_directory, args.device, weights=args.weights,
        half=args.half, chunksize=args.chunksize, use_rt=args.cudart,
    )

    samples = 0
    num_reads = 0
    max_read_size = 4e6
    dtype = np.float16 if args.half else np.float32
    reader = PreprocessReader(args.reads_directory)
    writer = DecoderWriterPool(model, beamsize=args.beamsize, fastq=args.fastq)

    t0 = time.perf_counter()
    sys.stderr.write("> calling\n")

    with writer, reader, torch.no_grad():

        while True:

            read = reader.queue.get()
            if read is None:
                break

            if len(read.signal) > max_read_size:
                sys.stderr.write("> skipping long read %s (%s samples)\n" % (read.read_id, len(read.signal)))
                continue

            num_reads += 1
            samples += len(read.signal)

            raw_data = torch.tensor(read.signal.astype(dtype))
            chunks = chunk(raw_data, args.chunksize, args.overlap)

            posteriors_ = model(chunks.to(args.device)).cpu().numpy()
            posteriors = stitch(posteriors_, args.overlap // model.stride // 2)

            writer.queue.put((read, posteriors[:raw_data.shape[0]]))

    duration = time.perf_counter() - t0

    sys.stderr.write("> completed reads: %s\n" % num_reads)
    sys.stderr.write("> duration: %s\n" % timedelta(seconds=np.round(duration)))
    sys.stderr.write("> samples per second %.1E\n" % (samples / duration))
    sys.stderr.write("> done\n")


def argparser():
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        add_help=False
    )
    parser.add_argument("model_directory")
    parser.add_argument("reads_directory")
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--weights", default="0", type=str)
    parser.add_argument("--beamsize", default=5, type=int)
    parser.add_argument("--chunksize", default=0, type=int)
    parser.add_argument("--overlap", default=0, type=int)
    parser.add_argument("--half", action="store_true", default=half_supported())
    parser.add_argument("--fastq", action="store_true", default=False)
    parser.add_argument("--cudart", action="store_true", default=False)
    return parser

if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()
    main(args)
