# Neural Network-based Variant Calling (NN-VARIANT)

## Based on [Clair3](https://github.com/HKU-BAL/Clair3)
Symphonizing pileup and full-alignment for high-performance long-read variant calling

If you find `NN-VARIANT` useful, please cite:

```
@article{Zheng2021,
  doi = {10.1101/2021.12.29.474431},
  url = {https://doi.org/10.1101/2021.12.29.474431},
  year = {2021},
  month = dec,
  publisher = {Cold Spring Harbor Laboratory},
  author = {Zhenxian Zheng and Shumin Li and Junhao Su and Amy Wing-Sze Leung and Tak-Wah Lam and Ruibang Luo},
  title = {Symphonizing pileup and full-alignment for deep learning-based long-read variant calling}
}
```

## Installation

NN-Variant requires the following Python3 packages (requirements.txt):

| Package           | Version  |
| ----------------- | -------- |
| Python            | 3.6.10   |
| TF                | 2.1.0    |
| pypy              | 3.6      |
| intervaltree      | 3.0.2    |
| mpmath            | 1.2.1    |
| tensorflow-addons | 0.11.2   |
| tables            | 3.6.1    |
| pigz              | 4.4      |
| parallel          | 20191122 |
| zstd              | 1.4.4    |
| samtools          | 1.1.0    |
| whatshap          | 1.0      |
| ensurepip         |          |

## Execution

```
./Clair3/callVar.sh --bam_fn=<file.bam> --ref_fn=<file.fa> --threads=<nthreads> --platform={ont,hifi,ilmn}| \
	          --model_path=\"r941_prom_hac_g360+g422\" \
	          --bed_fn=region.bed \
              --output=. \
              --chunk_size=\$((${contig_len}/\${OMP_NUM_THREADS}+1)) \
              \$PYPY"
```

