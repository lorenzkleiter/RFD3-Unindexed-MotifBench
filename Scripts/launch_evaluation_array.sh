#!/bin/bash
# Evaluates all 30 problems' converted scaffolds against the standard
# MotifBench pipeline. One array task per problem, %10 concurrency.
#
# NOT SUBMITTED YET -- review before running: sbatch launch_evaluation_array.sh
#
# Must match the SUFFIX used for the corresponding run_rfd3_array.sh /
# convert_all.sh submissions, e.g.: SUFFIX=_rep2 sbatch launch_evaluation_array.sh
# A config<SUFFIX>.txt is generated on the fly (each array task writes the same
# content, harmless) rather than requiring a pre-made static file per replicate.

#SBATCH --job-name=Unindexed_eval_fixed
#SBATCH --partition=barracuda
#SBATCH --reservation=AIPD
#SBATCH --mem=24G
#SBATCH --gres=gpu:1
#SBATCH --time=1-00:00:00
#SBATCH --array=0-29%4
#SBATCH --output=/zfs/s01/z04/home/kleiter/Projects/MotifBench/RFD3/Unindexed_Fixed/benchmark/logs/eval_%A_%a.out

SUFFIX=${SUFFIX:-}

UNINDEXED_FIXED_DIR="/zfs/s01/z04/home/kleiter/Projects/MotifBench/RFD3/Unindexed_Fixed"
EVALUATE_SCRIPT="/work/AIPDlab/software/MotifBench/scripts/evaluate_bbs.sh"
CONFIG_PATH="${UNINDEXED_FIXED_DIR}/config${SUFFIX}.txt"

# Written per array task, so use a per-task temp file + atomic mv -- all 30
# tasks write identical content, but concurrent writes to the same path
# without this could otherwise interleave into a corrupt file.
TMP_CONFIG="$(mktemp "${UNINDEXED_FIXED_DIR}/.config${SUFFIX}.XXXXXX")"
cat > "${TMP_CONFIG}" <<EOF
scaffold_base_dir=${UNINDEXED_FIXED_DIR}/scaffolds_pdb${SUFFIX}
benchmark_dir=/work/AIPDlab/software/MotifBench
foldseek_db_path=/work/AIPDlab/data/foldseek_database/pdb
base_output_dir=${UNINDEXED_FIXED_DIR}/results${SUFFIX}
python_path=/zfs/s01/z04/home/kleiter/.conda/envs/motif_bench/bin/python
EOF
mv "${TMP_CONFIG}" "${CONFIG_PATH}"

MOTIFS=(01_1LDB 02_1ITU 03_2CGA 04_5WN9 05_5ZE9 06_6E6R 07_6E6R 08_7AD5 09_7CG5 10_7WRK 11_3TQB 12_4JHW 13_4JHW 14_5IUS 15_7A8S 16_7BNY 17_7DGW 18_7MQQ 19_7MQQ 20_7UWL 21_1B73 22_1BCF 23_1MPY 24_1QY3 25_2RKX 26_3B5V 27_4XOJ 28_5YUI 29_6CPA 30_7UWL)
motif_name="${MOTIFS[$SLURM_ARRAY_TASK_ID]}"

echo "========================================"
echo " Motif  : ${motif_name}"
echo " Suffix : ${SUFFIX:-(none)}"
echo " Config : ${CONFIG_PATH}"
echo "========================================"

cd /work/AIPDlab/software/MotifBench/scripts
bash "${EVALUATE_SCRIPT}" "${motif_name}" "${CONFIG_PATH}"
EXIT_CODE=$?

if [[ ${EXIT_CODE} -eq 0 ]]; then
    echo "Task finished successfully — $(date)"
else
    echo "Task FAILED with exit code ${EXIT_CODE} — $(date)"
fi
exit ${EXIT_CODE}
