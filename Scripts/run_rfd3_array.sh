#!/bin/bash
# ==============================================================================
# RFD3 SLURM Array Job -- Unindexed benchmark, all 30 problems, PATCHED fork.
# Unlike AutoContigmap/ExpertGuess (explicit contig recipe: segment order +
# gap sizes), this marks the whole motif region as "unindexed" (multi-chain
# structure preserved as in ExpertGuess's motif_pdbs_gly) and lets RFD3 decide
# placement itself. Runs against the patched RFD3 fork
# (github.com/lorenzkleiter/RFD3-Fix, conda env fix_rfd3) instead of stock
# RFD3. One array task per motif.
#
# NOT SUBMITTED YET -- review before running: sbatch run_rfd3_array.sh
#
# For triplicates, set SUFFIX before submitting so each replicate writes to
# its own scaffolds<SUFFIX> dir instead of clobbering the others, e.g.:
#   SUFFIX=_rep1 sbatch run_rfd3_array.sh   (or leave SUFFIX unset for rep1)
#   SUFFIX=_rep2 sbatch run_rfd3_array.sh
#   SUFFIX=_rep3 sbatch run_rfd3_array.sh
# (RFD3's own seed is left unset/null in the input JSONs, so each submission
# already gets an independently random seed regardless of SUFFIX -- SUFFIX
# only controls where the output goes.)
# ==============================================================================

#SBATCH --job-name=Unindexed_RFD3_Fixed
#SBATCH -c 1
#SBATCH --nodes=1
#SBATCH --partition=barracuda
#SBATCH --nodelist=t38cn054,t38cn055
#SBATCH --mem=16G
#SBATCH --gres=gpu:1
#SBATCH --time=10:00:00
#SBATCH --array=0-29%3
#SBATCH --output=/zfs/s01/z04/home/kleiter/Projects/MotifBench/RFD3/Unindexed_Fixed/benchmark/logs/rfd3_%A_%a.out

SUFFIX=${SUFFIX:-}

INPUTS_DIR="/zfs/s01/z04/home/kleiter/Projects/MotifBench/RFD3/Unindexed_Fixed/benchmark/input_jsons"
SCAFFOLD_DIR="/zfs/s01/z04/home/kleiter/Projects/MotifBench/RFD3/Unindexed_Fixed/scaffolds${SUFFIX}"
CKPT_PATH="/work/kleiter/checkpoints/rfd3_latest.ckpt"

N_BATCHES=${N_BATCHES:-100}
DIFFUSION_BATCH_SIZE=1

MOTIFS=(01_1LDB 02_1ITU 03_2CGA 04_5WN9 05_5ZE9 06_6E6R 07_6E6R 08_7AD5 09_7CG5 10_7WRK 11_3TQB 12_4JHW 13_4JHW 14_5IUS 15_7A8S 16_7BNY 17_7DGW 18_7MQQ 19_7MQQ 20_7UWL 21_1B73 22_1BCF 23_1MPY 24_1QY3 25_2RKX 26_3B5V 27_4XOJ 28_5YUI 29_6CPA 30_7UWL)
motif_name="${MOTIFS[$SLURM_ARRAY_TASK_ID]}"

INPUT_JSON="${INPUTS_DIR}/${motif_name}.json"
OUTPUT_DIR="${SCAFFOLD_DIR}/${motif_name}"
mkdir -p "${OUTPUT_DIR}"

echo "========================================"
echo " Unindexed RFD3 Fixed Array Task Starting"
echo " $(date)"
echo "========================================"
echo " Array task : ${SLURM_ARRAY_TASK_ID}"
echo " Motif      : ${motif_name}"
echo " Suffix     : ${SUFFIX:-(none)}"
echo " Input JSON : ${INPUT_JSON}"
echo " Output dir : ${OUTPUT_DIR}"
echo " N_BATCHES  : ${N_BATCHES}"
echo "========================================"

if [[ ! -f "${INPUT_JSON}" ]]; then
    echo "ERROR: Input JSON not found: ${INPUT_JSON}"
    exit 1
fi

conda run -n fix_rfd3 rfd3 design \
    out_dir="${OUTPUT_DIR}" \
    inputs="${INPUT_JSON}" \
    prevalidate_inputs=True \
    n_batches="${N_BATCHES}" \
    diffusion_batch_size="${DIFFUSION_BATCH_SIZE}" \
    ckpt_path="${CKPT_PATH}"

EXIT_CODE=$?

echo "========================================"
if [[ ${EXIT_CODE} -eq 0 ]]; then
    echo " Task finished successfully — $(date)"
else
    echo " Task FAILED with exit code ${EXIT_CODE} — $(date)"
fi
echo "========================================"

exit ${EXIT_CODE}
