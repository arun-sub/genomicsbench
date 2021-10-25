<p align="center"><img src="https://github.com/arun-sub/punnet/blob/master/img/GenomicsBenchLogo-Colored.png" width="750"></p>

# About

A benchmark suite covering the major steps in short and long-read genome sequence analysis pipelines such as basecalling, sequence mapping, de-novo assembly, variant calling and polishing.

# Usage

1) Clone this repository and scalar branch:
`git clone -b scalar --recursive https://github.com/arun-sub/genomicsbench.git`

2) Change directory to the specific benchmark:
`cd genomicsbench/scalar-benchmarks/`
`cd bsw # chain, fmi, phmm`

3) Compile the code:
`make`

4) Run using small or large dataset:
`sbatch run_test small # large`

Look at the individual folders in `scalar-benchmarks` for more information.
