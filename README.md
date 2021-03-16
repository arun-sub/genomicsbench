<p align="center"><img src="https://github.com/arun-sub/punnet/blob/master/img/GenomicsBenchLogo-Colored.png" width="750"></p>

# GenomicsBench

A benchmark suite for genomics

## Download

```bash
git clone --recursive https://github.com/arun-sub/genomicsbench.git
```

## Compilation (for CPU benchmarks)

* RHEL/Fedora system prerequisites

```bash
sudo yum -y install $(cat rhel.prerequisites)
```
* Debian system prerequisites

```bash
sudo apt-get install $(cat debian.prerequisites)
```

* Compile

`make -j<num_threads>`

## Python setup (for GPU benchmarks)

To run Python-based benchmarks nn-base and nn-variant, follow the steps below:

* Download and install miniconda from [this](https://docs.conda.io/projects/continuumio-conda/en/latest/user-guide/install/download.html) link.

* Follow the steps below to set up a conda environment:

```bash
# make sure channels are added in conda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# create conda environment named "genomicsbench"
conda create -n genomicsbench -c bioconda clair
conda activate genomicsbench
conda install deepdish

pip install --upgrade pip
pip install -r requirements.txt
pypy3 -m ensurepip
pypy3 -m pip install --no-cache-dir intervaltree==3.0.2
```

## Citation

If you use GenomicsBench or find GenomicsBench useful, please cite this work:

> **Arun Subramaniyan, Yufeng Gu, Timothy Dunn, Somnath Paul, Md. Vasimuddin, Sanchit Misra, David Blaauw, Satish Narayanasamy, Reetuparna Das. *GenomicsBench: A Benchmark Suite for Genomics*, In IEEE International Symposium on Performance Analysis of Systems and Software (ISPASS), 2021 (to appear)**

## Issues and bug reporting

GenomicsBench is under active development and we appreciate any feedback and suggestions from the community. Feel free to raise an issue or submit a pull request on Github. For assistance in using GenomicsBench, contact the authors: Arun Subramaniyan (arunsub@umich.edu), Yufeng Gu (yufenggu@umich.edu), Timothy Dunn (timdunn@umich.edu).

