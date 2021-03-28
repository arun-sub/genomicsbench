import os
import sys
from time import time
import numpy as np
import deepdish as dd
import shared.param as param
from clair.model import Clair
from argparse import ArgumentParser


def prediction(args, m):

    print("Begin predicting...")
    prediction_output = []
    input_mini_match = dd.io.load(args.input_fn)
    output_mini_match = dd.io.load(args.output_fn)
    time_counter = {"Load_mini_batch": [],
                    "Model_prediction": [],
                    "Write_batch_to_output": []}

    begin_time = time()
    for i in range(len(input_mini_match)):
        mini_batch = input_mini_match[i]
        X, _ = mini_batch
        tmp_time = time()
        m.predict(X)
        cost_time = time() - tmp_time
        #print(cost_time)
        time_counter["Model_prediction"].append(round(cost_time, 4))
        prediction_output.append(m.prediction)

    end_time = time() - begin_time

    comp = []
    #for i in range(len(input_mini_match)):
    #    print(prediction_output[i][0], output_mini_match[i][0])
    #    comp.append(np.all(np.round(prediction_output[i][0], 3) == np.round(output_mini_match[i][0], 3)))

    #print(comp)
    #if False not in comp:
    #    print("My_prediction function is correct, which takes %.4f s" % end_time)
    #else:
    #    print("My_prediction function is wrong, which takes %.4f s" % end_time)
    #dd.io.save("time_counter_my_prediction.h5", time_counter)
    print("Time taken: %.4f s" % end_time)

def Run(args):

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"

    if args.threads is None:
        if args.tensor_fn == "PIPE":
            param.NUM_THREADS = 4
    else:
        param.NUM_THREADS = args.threads
        param.NUM_THREADS -= 1
        if param.NUM_THREADS < 1:
            param.NUM_THREADS = 1

    m = Clair()
    m.init()
    m.restore_parameters(os.path.abspath(args.chkpnt_fn))

    prediction(args, m)


def main():
    parser = ArgumentParser(description="Call variants using a trained model and tensors of candididate variants")

    parser.add_argument('--input_fn', type=str, default="prediction_input.h5",
                        help="input file")

    parser.add_argument('--output_fn', type=str, default="prediction_output.h5",
                        help="output file")

    parser.add_argument('--tensor_fn', type=str, default="PIPE",
                        help="Tensor input, use PIPE for standard input")

    parser.add_argument('--chkpnt_fn', type=str, default=None,
                        help="Input a checkpoint for testing")

    parser.add_argument('--call_fn', type=str, default=None,
                        help="Output variant predictions")

    parser.add_argument('--bam_fn', type=str, default="bam.bam",
                        help="BAM file input, default: %(default)s")

    parser.add_argument('--qual', type=int, default=None,
                        help="If set, variant with equal or higher quality will be marked PASS, or LowQual otherwise, optional")

    parser.add_argument('--sampleName', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the VCF file")

    parser.add_argument('--showRef', action='store_true',
                        help="Show reference calls, optional")

    parser.add_argument('--debug', action='store_true',
                        help="Debug mode, optional")

    parser.add_argument('--ref_fn', type=str, default=None,
                        help="Reference fasta file input, optional, print contig tags in the VCF header if set")

    parser.add_argument('--threads', type=int, default=None,
                        help="Number of threads, optional")


    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    Run(args)


if __name__ == "__main__":
    main()
