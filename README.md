<p align="center"><img src="https://github.com/arun-sub/punnet/blob/master/img/GenomicsBenchLogo-Colored.png" width="750"></p>

# About

A benchmark suite covering the major steps in short and long-read genome sequence analysis pipelines such as basecalling, sequence mapping, de-novo assembly, variant calling and polishing.

## Download

* Latest source code

```bash
git clone --recursive https://github.com/arun-sub/genomicsbench.git
```

* Input datasets

```bash
wget https://genomicsbench.eecs.umich.edu/input-datasets.tar.gz
```

## Prequisites

* RHEL/Fedora system prerequisites

```bash
sudo yum -y install $(cat rhel.prerequisites)
```
* Debian system prerequisites

```bash
sudo apt-get install $(cat debian.prerequisites)
```

## Python setup (optional: only needed for GPU benchmarks)

To run Python-based benchmarks nn-base and nn-variant, follow the steps below:

* Download and install miniconda from [this](https://docs.conda.io/projects/continuumio-conda/en/latest/user-guide/install/download.html) link.

* Follow the steps below to set up a conda environment:

```bash
# make sure channels are added in conda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# create conda environment named "genomicsbench"
conda create -n genomicsbench -c bioconda clair python==3.6.8
conda activate genomicsbench
conda install deepdish

pip install --upgrade pip
pip install -r requirements.txt
pypy3 -m ensurepip
pypy3 -m pip install --no-cache-dir intervaltree==3.0.2
```

## Compile

* CPU benchmarks

```bash
make -j<num_threads>
```

Notes: 

- MKLROOT and MKL_IOMPS_DIR variables need to be set in Makefile to run `grm`
- VTUNE_HOME variable needs to be set if you want to run any VTune based analyses

* GPU benchmarks

Set CUDA_LIB=/usr/local/cuda or to the path of the local CUDA installation in Makefile. 
Also ensure environment variables PATH and LD_LIBRARY_PATH include the path to CUDA binaries and libraries.

```bash
make -j<num_threads> gpu
```

## Running

* CPU benchmarks

```bash
cd scripts
chmod +x ./run_cpu.sh
./run_cpu.sh <path to input dataset folder> <input size to run: small | large>
```

* GPU benchmarks

```bash
cd scripts
chmod +x ./run_gpu.sh
./run_gpu.sh <path to input dataset folder> <input size to run: small | large>
```

## Citation

If you use GenomicsBench or find GenomicsBench useful, please cite this work:

> **Arun Subramaniyan, Yufeng Gu, Timothy Dunn, Somnath Paul, Md. Vasimuddin, Sanchit Misra, David Blaauw, Satish Narayanasamy, Reetuparna Das. *GenomicsBench: A Benchmark Suite for Genomics*, In IEEE International Symposium on Performance Analysis of Systems and Software (ISPASS), 2021 (to appear)**

```
@inproceedings{genomicsbench,
    title={GenomicsBench: A Benchmark Suite for Genomics}},
    author={Subramaniyan, Arun and Gu, Yufeng and Dunn, Timothy and Paul, Somnath and Vasimuddin, Md. and Misra, Sanchit and Blaauw, David and Narayanasamy, Satish and Das, Reetuparna},
    booktitle={Proceedings of the IEEE International Symposium on Performance Analysis of Systems and Software (ISPASS)},
    year={2021},
}
```

## Issues and bug reporting

GenomicsBench is under active development and we appreciate any feedback and suggestions from the community. Feel free to raise an issue or submit a pull request on Github. For assistance in using GenomicsBench, please contact: Arun Subramaniyan (arunsub@umich.edu), Yufeng Gu (yufenggu@umich.edu), Timothy Dunn (timdunn@umich.edu)

## Licensing

Each benchmark is individually licensed according to the tool it is extracted from.

## Acknowledgement

This work was supported in part by Precision Health at the University of Michigan, by the Kahn foundation, by the NSF under the CAREER-1652294 award and the Applications Driving Architectures (ADA) Research Center, a JUMP Center co-sponsored by SRC and DARPA.

