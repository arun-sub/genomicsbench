#!/bin/bash

reads=10

mkdir -p ./output
/usr/local/cuda/bin/nvprof \
    --metrics all \
    python3 bonito/basecall.py \
        models/bonito_dna_r941 \
        data/$reads \
        --half \
        --fastq \
        > output/bonito.fastq 

# python3 bonito/basecall.py \
#     models/bonito_dna_r941 \
#     data/$reads \
#     --half \
#     --fastq \
#     --chunksize 3000 \
#     --cudart \
#     > output/bonito.fastq
