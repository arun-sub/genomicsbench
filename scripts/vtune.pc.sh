#!/bin/bash

set -e

usage() {
	echo -e "\n Usage $0 <INPUTS_DIR> <OUTPUTS_DIR>\n"
}


vtune_pc() {
	vtune -collect-with runsa -start-paused -r $1 -knob event-config=INST_RETIRED.ANY,CPU_CLK_UNHALTED.THREAD,MEM_INST_RETIRED.ALL_LOADS:sa=100000,MEM_LOAD_RETIRED.L1_HIT:sa=100000,MEM_LOAD_RETIRED.L2_HIT:sa=100,MEM_LOAD_RETIRED.L3_HIT:sa=100,MEM_LOAD_RETIRED.L1_MISS:sa=100,MEM_LOAD_RETIRED.L2_MISS:sa=100,MEM_LOAD_RETIRED.L3_MISS:sa=100,OFFCORE_RESPONSE:request=DEMAND_DATA_RD:response=L3_MISS.SNOOP_MISS:sa=1000,OFFCORE_RESPONSE:request=DEMAND_DATA_RD:response=L3_MISS.SNOOP_HIT_NO_FWD:sa=1000,OFFCORE_RESPONSE:request=DEMAND_CODE_RD:response=L3_MISS.SNOOP_MISS:sa=1000,OFFCORE_RESPONSE:request=DEMAND_CODE_RD:response=L3_MISS.SNOOP_HIT_NO_FWD:sa=1000,OFFCORE_REQUESTS.L3_MISS_DEMAND_DATA_RD:sa=1000,OFFCORE_REQUESTS.DEMAND_DATA_RD:sa=1000,OFFCORE_REQUESTS:ALL_DATA_RD:sa=1000,OFFCORE_RESPONSE:request=OTHER:response=L3_MISS_LOCAL_DRAM.SNOOP_HIT_NO_FWD:sa=1000,OFFCORE_RESPONSE:request=OTHER:response=L3_MISS_LOCAL_DRAM.SNOOP_MISS:sa=1000,LOAD_HIT_PRE.SW_PF:sa=100,BR_INST_RETIRED.ALL_BRANCHES,BR_MISP_RETIRED.ALL_BRANCHES,DTLB_LOAD_MISSES.MISS_CAUSES_A_WALK,CYCLE_ACTIVITY.STALLS_L3_MISS,CYCLE_ACTIVITY.STALLS_MEM_ANY -- $2
}

if [[ ( $# == "--help" ) || ( $# == "-h" ) ]]
then
	usage
	exit 0
fi

if [[ $# -lt 2 ]]
then
	usage
	exit 1
fi

INPUTS_DIR=$1
OUTPUTS_DIR=$2

echo "Running fmi"
vtune_pc $OUTPUTS_DIR/fmi_pc "../benchmarks/fmi/fmi $INPUTS_DIR/fmi/broad $INPUTS_DIR/fmi/large/SRR7733443_10m_1.fastq 512 19 1"

echo "Running bsw"
vtune_pc $OUTPUTS_DIR/bsw_pc "../benchmarks/bsw/bsw -pairs $INPUTS_DIR/bsw/large/banded_SRR7733443_1m_input.txt -t 1 -b 512"

echo "Running phmm"
export LD_LIBRARY_PATH=../tools/GKL/build/native:$LD_LIBRARY_PATH
vtune_pc $OUTPUTS_DIR/phmm_pc "../benchmarks/phmm/phmm -f $INPUTS_DIR/phmm/large/large.in -t 1"

echo "Running dbg"
vtune_pc $OUTPUTS_DIR/dbg_pc "../benchmarks/dbg/dbg $INPUTS_DIR/dbg/large/ERR194147-mem2-chr22.bam chr22:0-50818468 $INPUTS_DIR/dbg/large/Homo_sapiens_assembly38.fasta 1"

echo "Running chain"
vtune_pc $OUTPUTS_DIR/chain_pc "../benchmarks/chain/chain -i $INPUTS_DIR/chain/large/c_elegans_40x.10k.in -o $INPUTS_DIR/chain/large/c_elegans_40x.10k.out"

echo "Running poa"
vtune_pc $OUTPUTS_DIR/poa_pc "../benchmarks/poa/poa -s $INPUTS_DIR/poa/large/input.fasta -t 1"

echo "Running kmer-cnt"
vtune_pc $OUTPUTS_DIR/kmer-cnt_pc "../benchmarks/kmer-cnt/kmer-cnt --reads $INPUTS_DIR/kmer-cnt/large/Loman_E.coli_MAP006-1_2D_50x.fasta --config ../tools/Flye/flye/config/bin_cfg/asm_raw_reads.cfg --threads 1 --debug"

echo "Running pileup"
vtune_pc $OUTPUTS_DIR/pileup_pc "../benchmarks/pileup/pileup $INPUTS_DIR/pileup/large/HG002_prom_R941_guppy360_2_GRCh38_ch20.bam chr20:1-64444167 1"

echo "Running grm"
export LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/2021.1.1/lib/intel64:/opt/intel/oneapi/compiler/2021.1.2/linux/compiler/lib/intel64_lin:$LD_LIBRARY_PATH
vtune_pc $OUTPUTS_DIR/grm_pc "../benchmarks/grm/2.0/build_dynamic/plink2 --maf 0.01 --pgen $INPUTS_DIR/grm/large/chr1_phase3.pgen --pvar $INPUTS_DIR/grm/large/chr1_phase3.pvar --psam $INPUTS_DIR/grm/large/phase3_corrected.psam --make-grm-bin --out $INPUTS_DIR/grm/large/grm --threads 1"

