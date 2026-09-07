#!/bin/bash
# ==============================================================================
# Glue script: launches all 3 replicates of the Unindexed_Fixed pipeline, each
# as a 3-stage chain (generation -> convert -> evaluate), 9 SLURM jobs total.
# The 3 stages within a replicate are dependent on each other (afterok); the
# 3 replicates themselves are independent of one another and run in parallel.
# No external cross-job dependency -- all three stages are pinned to
# t38cn054/t38cn055 (via --nodelist in each .sh, no --reservation) and just
# scheduled normally against whatever else those two nodes are running.
#
# NOT SUBMITTED YET -- review before running: bash launch_triplicates.sh
# ==============================================================================

set -e

BENCH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATH="/apps/slurm/latest/bin:$PATH"

for SUFFIX in _rep1 _rep2 _rep3; do
    echo "========================================"
    echo " Replicate: ${SUFFIX}"
    echo "========================================"

    JOB1=$(sbatch --parsable --export=ALL,SUFFIX="${SUFFIX}" "${BENCH_DIR}/run_rfd3_array.sh")
    echo "  run_rfd3_array.sh          -> ${JOB1}"

    JOB2=$(sbatch --parsable --export=ALL,SUFFIX="${SUFFIX}" --dependency=afterok:"${JOB1}" "${BENCH_DIR}/convert_all.sh")
    echo "  convert_all.sh             -> ${JOB2}  (depends on ${JOB1})"

    JOB3=$(sbatch --parsable --export=ALL,SUFFIX="${SUFFIX}" --dependency=afterok:"${JOB2}" "${BENCH_DIR}/launch_evaluation_array.sh")
    echo "  launch_evaluation_array.sh -> ${JOB3}  (depends on ${JOB2})"

    echo ""
done

echo "All 9 jobs submitted (3 replicates x 3 chained stages)."
