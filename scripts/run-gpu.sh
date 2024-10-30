#!/bin/bash

set -e

usage() {
	echo -e "\n Usage $0 <INPUTS_DIR> <INPUT_SIZE>\n\n Example: $0 [../input-datasets] [small | large]\n"
}

if [[ ( $# == "--help" ) || ( $# == "-h" ) ]]
then
	usage
	exit 0
fi

if [[ $# -lt 1 ]]
then
	usage
	exit 1
fi

INPUTS_DIR=$1
INPUTS_SIZE=$2

if [[ ( $INPUTS_SIZE == "small" ) ]]
then

	echo "Running nn-base"
	python ../benchmarks/nn-base/bonito/basecall.py ../benchmarks/nn-base/models/bonito_dna_r941 $INPUTS_DIR/nn-base/small --device cuda:0 --fastq > $INPUTS_DIR/nn-base/small/out-small.fastq

	echo "Running nn-variant"
	../benchmarks/nn-variant/Clair3/callVar.sh --bam_fn="$INPUTS_DIR/nn-variant/HG002_GRCh38_ONT-UL_GIAB_20200122_chr20_0_10000000.phased.bam" --ref_fn="$INPUTS_DIR/nn-variant/hg38_chr20.fa" --threads=1 --platform="ont" --model_path="$INPUTS_DIR/nn-variant/models/r941_prom_hac_g360+g422" --bed_fn="$INPUTS_DIR/nn-variant/small/region.bed" --output="$INPUTS_DIR/nn-variant/output-small"

	echo "Running abea"
	../benchmarks/abea/f5c eventalign -b $INPUTS_DIR/abea/small/1000reads.bam -g $INPUTS_DIR/abea/humangenome.fa -r $INPUTS_DIR/abea/1000reads.fastq -B 3.7M > $INPUTS_DIR/abea/small/events.tsv

else

	echo "Running nn-base"
	python ../benchmarks/nn-base/bonito/basecall.py ../benchmarks/nn-base/models/bonito_dna_r941 $INPUTS_DIR/nn-base/large --device cuda:0 --fastq > $INPUTS_DIR/nn-base/large/out-large.fastq

	echo "Running nn-variant"
	../benchmarks/nn-variant/Clair3/callVar.sh --bam_fn="$INPUTS_DIR/nn-variant/HG002_GRCh38_ONT-UL_GIAB_20200122_chr20_0_10000000.phased.bam" --ref_fn="$INPUTS_DIR/nn-variant/hg38_chr20.fa" --threads=1 --platform="ont" --model_path="$INPUTS_DIR/nn-variant/models/r941_prom_hac_g360+g422" --bed_fn="$INPUTS_DIR/nn-variant/large/region.bed" --output="$INPUTS_DIR/nn-variant/output-large"

	echo "Running abea"
	../benchmarks/abea/f5c eventalign -b $INPUTS_DIR/abea/large/10000reads.bam -g $INPUTS_DIR/abea/humangenome.fa -r $INPUTS_DIR/abea/10000reads.fastq -B 3.7M > $INPUTS_DIR/abea/large/events.tsv

fi
