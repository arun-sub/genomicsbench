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

if [[ ( $INPUTS_SIZE == "large" ) ]]
then
	echo "Running fmi"
	../benchmarks/fmi/fmi $INPUTS_DIR/fmi/broad $INPUTS_DIR/fmi/large/SRR7733443_10m_1.fastq 512 19 1

	echo "Running bsw"
	../benchmarks/bsw/bsw -pairs $INPUTS_DIR/bsw/large/bandedSWA_SRR7733443_1m_input.txt -t 1 -b 512

	echo "Running phmm"
	export LD_LIBRARY_PATH=../tools/GKL/build/native:$LD_LIBRARY_PATH
	../benchmarks/phmm/phmm -f $INPUTS_DIR/phmm/large/large.in -t 1

	echo "Running dbg"
	../benchmarks/dbg/dbg $INPUTS_DIR/dbg/large/ERR194147-mem2-chr22.bam chr22:0-50818468 $INPUTS_DIR/dbg/large/Homo_sapiens_assembly38.fasta 1

	echo "Running chain"
	../benchmarks/chain/chain -i $INPUTS_DIR/chain/large/c_elegans_40x.10k.in -o $INPUTS_DIR/chain/large/c_elegans_40x.10k.out

	echo "Running fast-chain"
	../benchmarks/fast-chain/chain -i $INPUTS_DIR/chain/large/c_elegans_40x.10k.in -o $INPUTS_DIR/fast-chain/large/c_elegans_40x.10k.out

	echo "Running poa"
	../benchmarks/poa/poa -s $INPUTS_DIR/poa/large/input.fasta -t 1

	echo "Running kmer-cnt"
	../benchmarks/kmer-cnt/kmer-cnt --reads $INPUTS_DIR/kmer-cnt/large/Loman_E.coli_MAP006-1_2D_50x.fasta --config ../tools/Flye/flye/config/bin_cfg/asm_raw_reads.cfg --threads 1 --debug

	echo "Running pileup"
	../benchmarks/pileup/pileup $INPUTS_DIR/pileup/large/HG002_prom_R941_guppy360_2_GRCh38_ch20.bam chr20:1-64444167 1 > $INPUTS_DIR/pileup/large/pileup.txt

	echo "Running grm"
	export LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/2021.1.1/lib/intel64:/opt/intel/oneapi/compiler/2021.1.2/linux/compiler/lib/intel64_lin:$LD_LIBRARY_PATH
	../benchmarks/grm/2.0/build_dynamic/plink2 --maf 0.01 --pgen $INPUTS_DIR/grm/large/chr1_phase3.pgen --pvar $INPUTS_DIR/grm/large/chr1_phase3.pvar --psam $INPUTS_DIR/grm/large/phase3_corrected.psam --make-grm-bin --out $INPUTS_DIR/grm/large/grm --threads 1

	echo "Running wfa"
	../benchmarks/wfa/bin/align_benchmark -i $INPUTS_DIR/bsw/large/bandedSWA_SRR7733443_1m_input.txt -o checksum.file -t 1
else

	echo "Running fmi"
	../benchmarks/fmi/fmi $INPUTS_DIR/fmi/broad $INPUTS_DIR/fmi/small/SRR7733443_1m_1.fastq 512 19 1

	echo "Running bsw"
	../benchmarks/bsw/bsw -pairs $INPUTS_DIR/bsw/small/bandedSWA_SRR7733443_100k_input.txt -t 1 -b 512

	echo "Running phmm"
	export LD_LIBRARY_PATH=../tools/GKL/build/native:$LD_LIBRARY_PATH
	../benchmarks/phmm/phmm -f $INPUTS_DIR/phmm/small/5m.in -t 1

	echo "Running dbg"
	../benchmarks/dbg/dbg $INPUTS_DIR/dbg/small/ERR194147-mem2-chr22.bam chr22:16000000-16500000 $INPUTS_DIR/dbg/large/Homo_sapiens_assembly38.fasta 1

	echo "Running chain"
	../benchmarks/chain/chain -i $INPUTS_DIR/chain/small/in-1k.txt -o $INPUTS_DIR/chain/small/out-1k.txt

	echo "Running fast-chain"
	../benchmarks/fast-chain/chain -i $INPUTS_DIR/chain/small/in-1k.txt -o $INPUTS_DIR/fast-chain/small/out-1k.txt

	echo "Running poa"
	../benchmarks/poa/poa -s $INPUTS_DIR/poa/small/input-1000.fasta -t 1

	echo "Running kmer-cnt"
	../benchmarks/kmer-cnt/kmer-cnt --reads $INPUTS_DIR/kmer-cnt/small/Loman_E.coli_MAP006-1_2D_50x_1000.fasta --config ../tools/Flye/flye/config/bin_cfg/asm_raw_reads.cfg --threads 1 --debug

	echo "Running pileup"
	../benchmarks/pileup/pileup $INPUTS_DIR/pileup/small/saureus.bam tig00000061:1-1499707 1 > $INPUTS_DIR/pileup/small/pileup.txt

	echo "Running grm"
	export LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/2021.1.1/lib/intel64:/opt/intel/oneapi/compiler/2021.1.2/linux/compiler/lib/intel64_lin:$LD_LIBRARY_PATH
	../benchmarks/grm/2.0/build_dynamic/plink2 --maf 0.01 --pgen $INPUTS_DIR/grm/small/chr22_phase3.pgen --pvar $INPUTS_DIR/grm/small/chr22_phase3.pvar --psam $INPUTS_DIR/grm/small/phase3_corrected.psam --make-grm-bin --out $INPUTS_DIR/grm/small/grm --threads 1

	echo "Running wfa"
	../benchmarks/wfa/bin/align_benchmark -i $INPUTS_DIR/bsw/small/bandedSWA_SRR7733443_100k_input.txt -o checksum.file -t 1
fi
