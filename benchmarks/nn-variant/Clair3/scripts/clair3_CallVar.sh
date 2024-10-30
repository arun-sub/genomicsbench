#!/bin/bash
SCRIPT_NAME=$(basename "$0")
Usage="Usage: ./${SCRIPT_NAME} --bam_fn=BAM --ref_fn=REF --output=OUTPUT_DIR --threads=THREADS --platform=PLATFORM --model_path=MODEL_PREFIX [--bed_fn=BED] [options]"
# INFO: whole calling workflow of clair3

set -e
ARGS=`getopt -o b:f:t:m:p:o:r::c::s::h::g \
-l bam_fn:,ref_fn:,threads:,model_path:,platform:,output:,\
bed_fn::,vcf_fn::,ctg_name::,sample_name::,help::,qual::,samtools::,python::,pypy::,parallel::,whatshap::,chunk_num::,chunk_size::,var_pct_full::,var_pct_phasing::,\
snp_min_af::,indel_min_af::,ref_pct_full::,pileup_only::,fast_mode::,gvcf::,print_ref_calls::,haploid_precise::,haploid_sensitive::,include_all_ctgs::,\
no_phasing_for_fa::,pileup_model_prefix::,fa_model_prefix::,call_snp_only::,remove_intermediate_dir::,enable_phasing::,enable_long_indel:: -n 'run_clair3.sh' -- "$@"`

if [ $? != 0 ] ; then echo"No input. Terminating...">&2 ; exit 1 ; fi
eval set -- "${ARGS}"

while true; do
   case "$1" in
    -b|--bam_fn ) BAM_FILE_PATH="$2"; shift 2 ;;
    -f|--ref_fn ) REFERENCE_FILE_PATH="$2"; shift 2 ;;
    -t|--threads ) THREADS="$2"; shift 2 ;;
    -m|--model_path ) MODEL_PATH="$2"; shift 2 ;;
    -p|--platform ) PLATFORM="$2"; shift 2 ;;
    -o|--output ) OUTPUT_FOLDER="$2"; shift 2 ;;
    --bed_fn ) BED_FILE_PATH="$2"; shift 2 ;;
    --vcf_fn ) VCF_FILE_PATH="$2"; shift 2 ;;
    --ctg_name ) CONTIGS="$2"; shift 2 ;;
    --sample_name ) SAMPLE="$2"; shift 2 ;;
    --chunk_num ) CHUNK_NUM="$2"; shift 2 ;;
    --chunk_size ) CHUNK_SIZE="$2"; shift 2 ;;
    --qual ) QUAL="$2"; shift 2 ;;
    --samtools ) SAMTOOLS="$2"; shift 2 ;;
    --python ) PYTHON="$2"; shift 2 ;;
    --pypy ) PYPY="$2"; shift 2 ;;
    --parallel ) PARALLEL="$2"; shift 2 ;;
    --whatshap ) WHATSHAP="$2"; shift 2 ;;
    --var_pct_full ) PRO="$2"; shift 2 ;;
    --ref_pct_full ) REF_PRO="$2"; shift 2 ;;
    --var_pct_phasing ) PHASING_PCT="$2"; shift 2 ;;
    --pileup_only ) PILEUP_ONLY="$2"; shift 2 ;;
    --fast_mode ) FAST_MODE="$2"; shift 2 ;;
    --call_snp_only ) SNP_ONLY="$2"; shift 2 ;;
    --print_ref_calls ) SHOW_REF="$2"; shift 2 ;;
    --gvcf ) GVCF="$2"; shift 2 ;;
    --snp_min_af ) SNP_AF="$2"; shift 2 ;;
    --indel_min_af ) INDEL_AF="$2"; shift 2 ;;
    --pileup_model_prefix ) PILEUP_PREFIX="$2"; shift 2 ;;
    --fa_model_prefix ) FA_PREFIX="$2"; shift 2 ;;
    --haploid_precise ) HAP_PRE="$2"; shift 2 ;;
    --haploid_sensitive ) HAP_SEN="$2"; shift 2 ;;
    --include_all_ctgs ) INCLUDE_ALL_CTGS="$2"; shift 2 ;;
    --no_phasing_for_fa ) NO_PHASING="$2"; shift 2 ;;
    --remove_intermediate_dir ) RM_TMP_DIR="$2"; shift 2 ;;
    --enable_phasing ) ENABLE_PHASING="$2"; shift 2 ;;
    --enable_long_indel ) ENABLE_LONG_INDEL="$2"; shift 2 ;;

    -- ) shift; break; ;;
    -h|--help ) print_help_messages; break ;;
    * ) print_help_messages; exit 0 ;;
   esac
done


SHELL_FOLDER=$(cd "$(dirname "$0")";pwd)
CLAIR3="${SHELL_FOLDER}/../clair3.py"

if [ ${BED_FILE_PATH} = "EMPTY" ] ; then BED_FILE_PATH= ; fi
RETRIES=4

PILEUP_CHECKPOINT_PATH="${MODEL_PATH}/${PILEUP_PREFIX}"
echo "PILEUP_CHECKPOINT_PATH ${PILEUP_CHECKPOINT_PATH}"
FULL_ALIGNMENT_CHECKPOINT_PATH="${MODEL_PATH}/${FA_PREFIX}"
LOG_PATH="${OUTPUT_FOLDER}/log"
TMP_FILE_PATH="${OUTPUT_FOLDER}/tmp"
SPLIT_BED_PATH="${TMP_FILE_PATH}/split_beds"
PILEUP_VCF_PATH="${TMP_FILE_PATH}/pileup_output"
GVCF_TMP_PATH="${TMP_FILE_PATH}/gvcf_tmp_output"
PHASE_OUTPUT_PATH="${TMP_FILE_PATH}/phase_output"
FULL_ALIGNMENT_OUTPUT_PATH="${TMP_FILE_PATH}/full_alignment_output"
PHASE_VCF_PATH="${PHASE_OUTPUT_PATH}/phase_vcf"
PHASE_BAM_PATH="${PHASE_OUTPUT_PATH}/phase_bam"
CANDIDATE_BED_PATH="${FULL_ALIGNMENT_OUTPUT_PATH}/candidate_bed"
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export OMP_NUM_THREADS=1

echo "[INFO] Check environment variables"
${PYTHON} ${CLAIR3} CheckEnvs \
    --bam_fn ${BAM_FILE_PATH} \
    --bed_fn ${BED_FILE_PATH} \
    --output_fn_prefix ${OUTPUT_FOLDER} \
    --ref_fn ${REFERENCE_FILE_PATH} \
    --vcf_fn ${VCF_FILE_PATH} \
    --ctg_name ${CONTIGS} \
    --chunk_num ${CHUNK_NUM} \
    --chunk_size ${CHUNK_SIZE} \
    --include_all_ctgs ${INCLUDE_ALL_CTGS} \
    --threads ${THREADS} \
    --python ${PYTHON} \
    --pypy ${PYPY} \
    --samtools ${SAMTOOLS} \
    --whatshap ${WHATSHAP} \
    --parallel ${PARALLEL} \
    --qual ${QUAL} \
    --sampleName ${SAMPLE} \
    --var_pct_full ${PRO} \
    --ref_pct_full ${REF_PRO} \
    --snp_min_af ${SNP_AF} \
    --indel_min_af ${INDEL_AF}
readarray -t CHR < "${OUTPUT_FOLDER}/tmp/CONTIGS"
if [ ${#CHR[@]} -eq 0 ]; then echo "[INFO] Exit in environment checking"; exit 0; fi
# GenarchBench: Tensorflow sometimes uses more than 1 thread, so Clair3 uses
# ${THREADS}*3/4 simultaneous threads as an heuristic. We obtained better scalability
# using all threads.
# THREADS_LOW=$((${THREADS}*3/4))
THREADS_LOW=${THREADS}
if [[ ${THREADS_LOW} < 1 ]]; then THREADS_LOW=1; fi

cd ${OUTPUT_FOLDER}
# Pileup calling
#-----------------------------------------------------------------------------------------------------------------------
# GenarchBnech: --delay = 0 (no need to wait at startup) and 
# --tensorflow_threads=1 (we want that each chunk is processed by a single physical CPU)

export CUDA_VISIBLE_DEVICES=""
echo "[INFO] 1/1 Call variants using pileup model"
\time -v ${PARALLEL} --retries ${RETRIES} -C ' ' --joblog ${LOG_PATH}/parallel_1_call_var_bam_pileup.log -j ${THREADS_LOW} \
"${PYTHON} ${CLAIR3} CallVarBam \
    --chkpnt_fn ${PILEUP_CHECKPOINT_PATH} \
    --bam_fn ${BAM_FILE_PATH} \
    --call_fn ${PILEUP_VCF_PATH}/pileup_{1}_{2}.vcf \
    --sampleName ${SAMPLE} \
    --ref_fn ${REFERENCE_FILE_PATH} \
    --extend_bed ${SPLIT_BED_PATH}/{1} \
    --bed_fn ${BED_FILE_PATH} \
    --vcf_fn ${VCF_FILE_PATH} \
    --ctgName {1} \
    --chunk_id {2} \
    --chunk_num {3} \
    --platform ${PLATFORM} \
    --fast_mode ${FAST_MODE} \
    --snp_min_af ${SNP_AF} \
    --indel_min_af ${INDEL_AF} \
    --call_snp_only ${SNP_ONLY} \
    --gvcf ${GVCF} \
    --enable_long_indel ${ENABLE_LONG_INDEL} \
    --python ${PYTHON} \
    --pypy ${PYPY} \
    --samtools ${SAMTOOLS} \
    --temp_file_dir ${GVCF_TMP_PATH} \
    --delay=0 \
    --tensorflow_threads=1 \
    --pileup" :::: ${OUTPUT_FOLDER}/tmp/CHUNK_LIST |& tee ${LOG_PATH}/1_call_var_bam_pileup.log

# Write to prof_pipe to tell the c++ wrapper to stop profiling.
echo "ee" > "prof_pipe"


