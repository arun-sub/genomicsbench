# Wavefront Alignment (WFA)

## 1. INTRODUCTION

### 1.1 What is WFA?

The wavefront alignment (WFA) algorithm is an exact gap-affine algorithm that takes advantage of  
homologous regions between the sequences to accelerate the alignment process. As opposed to 
traditional dynamic programming algorithms that run in quadratic time, the WFA runs in time O(ns),
proportional to the read length n and the alignment score s, using O(s^2) memory. Moreover, the WFA
exhibits simple data dependencies that can be easily vectorized, even by the automatic features of 
modern compilers, for different architectures, without the need to adapt the code.

This library implements the WFA and the WFA-Adapt algorithms for gap-affine penalties. It also 
provides support functions to display and verify the results. Moreover, it implements a benchmarking
tool that serves to evaluate the performance of these two algorithms, together with other 
high-performance alignment methods (checkout branch `benchmark`). The library can be executed   
through the benchmarking tool for evaluation purposes or can be integrated into your code by calling
the WFA functions.

If you are interested in benchmarking WFA with other algorithms implemented or integrated into the
WFA library, checkout branch `benchmark`.

### 1.2 Introduction to benchmarking WFA. Simple tests

The WFA includes the benchmarking tool *align-benchmark* to test performance of. This tool takes as 
input a dataset containing pairs of sequences (i.e., pattern and text) to align. Patterns are 
preceded by the '>' symbol and texts by the '<' symbol. Example:

```
>ATTGGAAAATAGGATTGGGGTTTGTTTATATTTGGGTTGAGGGATGTCCCACCTTCGTCGTCCTTACGTTTCCGGAAGGGAGTGGTTAGCTCGAAGCCCA
<GATTGGAAAATAGGATGGGGTTTGTTTATATTTGGGTTGAGGGATGTCCCACCTTGTCGTCCTTACGTTTCCGGAAGGGAGTGGTTGCTCGAAGCCCA
>CCGTAGAGTTAGACACTCGACCGTGGTGAATCCGCGACCACCGCTTTGACGGGCGCTCTACGGTATCCCGCGATTTGTGTACGTGAAGCAGTGATTAAAC
<CCTAGAGTTAGACACTCGACCGTGGTGAATCCGCGATCTACCGCTTTGACGGGCGCTCTACGGTATCCCGCGATTTGTGTACGTGAAGCGAGTGATTAAAC
[...]
```

You can either generate a custom dataset of your own, or use the *generate-dataset* tool to generate
a random dataset. For example, the following command generates a dataset named 'sample.dataset.seq' 
of 5M pairs of 100 bases with an alignment error of 5% (i.e., 5 mismatches, insertions or deletions 
per alignment).

```
$> ./bin/generate_dataset -n 5000000 -l 100 -e 0.05 -o sample.dataset.seq
```

Once you have the dataset ready, you can run the *align-benchmark* tool to benchmark the performance 
of a specific pairwise alignment method. For example, the WFA algorithm:

```
$> ./bin/align_benchmark -i sample.dataset.seq -a gap-affine-wfa
...processed 10000 reads (benchmark=125804.398 reads/s;alignment=188049.469 reads/s)
...processed 20000 reads (benchmark=117722.406 reads/s;alignment=180925.031 reads/s)
[...]
...processed 5000000 reads (benchmark=113844.039 reads/s;alignment=177325.281 reads/s)
[Benchmark]
=> Total.reads            5000000
=> Time.Benchmark        43.92 s  (    1   call,  43.92  s/call {min43.92s,Max43.92s})
  => Time.Alignment      28.20 s  ( 64.20 %) (    5 Mcalls,   5.64 us/call {min438ns,Max47.05ms})
```

The *align-benchmark* tool will finish and report overall benchmark time (including reading the 
input, setup, checking, etc.) and the time taken by the algorithm (i.e., *Time.Alignment*).

## 2. AUTHORS

  Santiago Marco-Sola \- santiagomsola@gmail.com     

## 3. REPORTING BUGS

Feedback and bug reporting it's highly appreciated. 
Please report any issue or suggestion on github, or by email to the main developer (santiagomsola@gmail.com).

## 4. LICENSE

`WFA` uses the same license as [WFA](https://github.com/smarco/WFA).

## 5. CITATION

**Santiago Marco-Sola, Juan Carlos Moure, Miquel Moreto, Antonio Espinosa**. ["Fast gap-affine pairwise alignment using the wavefront algorithm."](https://doi.org/10.1093/bioinformatics/btaa777) Bioinformatics, 2020.
