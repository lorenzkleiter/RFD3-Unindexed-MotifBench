#!/bin/bash
# Converts all 30 problems' raw RFD3 output (CIF -> PDB + scaffold_info.csv) in
# one pass, using convert_rfd3_outputs_unindexed.py's built-in "no motif_name
# given -> convert every subdirectory" mode. Unlike ExpertGuess/AutoContigmap,
# "unindex" mode never records a per-design contig (specification.extra.
# sampled_contig is just the total length, e.g. "200P" -- there's no explicit
# contig recipe to record). Instead it records `diffused_index_map` (input
# motif residue -> output residue), which convert_rfd3_outputs_unindexed.py
# (copied verbatim from the original Unindexed pipeline) uses to reconstruct
# scaffold_info.csv's motif placements. Still reuses ExpertGuess's
# contig_specifications.csv, but only for chain order per problem -- a
# property of the motif itself, unrelated to placement strategy.
#
# NOT SUBMITTED YET -- review before running: sbatch convert_all.sh
#
# Must match the SUFFIX used for the corresponding run_rfd3_array.sh submission,
# e.g.: SUFFIX=_rep2 sbatch convert_all.sh

#SBATCH --job-name=Unindexed_convert_all_fixed
#SBATCH --partition=barracuda
#SBATCH --nodelist=t38cn054,t38cn055
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=/zfs/s01/z04/home/kleiter/Projects/MotifBench/RFD3/Unindexed_Fixed/benchmark/logs/convert_all_%j.out

set -e

SUFFIX=${SUFFIX:-}

EXPERTGUESS_DIR="/zfs/s01/z04/home/kleiter/Projects/MotifBench/RFD3/ExpertGuess"
UNINDEXED_FIXED_DIR="/zfs/s01/z04/home/kleiter/Projects/MotifBench/RFD3/Unindexed_Fixed"
PY="/zfs/s01/z04/home/kleiter/.conda/envs/motif_bench/bin/python"

MOTIFS=(01_1LDB 02_1ITU 03_2CGA 04_5WN9 05_5ZE9 06_6E6R 07_6E6R 08_7AD5 09_7CG5 10_7WRK 11_3TQB 12_4JHW 13_4JHW 14_5IUS 15_7A8S 16_7BNY 17_7DGW 18_7MQQ 19_7MQQ 20_7UWL 21_1B73 22_1BCF 23_1MPY 24_1QY3 25_2RKX 26_3B5V 27_4XOJ 28_5YUI 29_6CPA 30_7UWL)

echo "========================================"
echo " Convert-all (Unindexed_Fixed) Starting"
echo " $(date)"
echo " Suffix: ${SUFFIX:-(none)}"
echo "========================================"

MISSING=0
for m in "${MOTIFS[@]}"; do
    n=$(ls "${UNINDEXED_FIXED_DIR}/scaffolds${SUFFIX}/${m}"/*.cif.gz 2>/dev/null | wc -l)
    echo "${m}: ${n} .cif.gz files"
    if [[ "${n}" -ne 100 ]]; then
        echo "  ERROR: expected 100, found ${n}"
        MISSING=1
    fi
done
if [[ "${MISSING}" -ne 0 ]]; then
    echo "ERROR: one or more motifs missing designs. Aborting."
    exit 1
fi

"${PY}" "${UNINDEXED_FIXED_DIR}/benchmark/convert_rfd3_outputs_unindexed.py" \
    "${UNINDEXED_FIXED_DIR}/scaffolds${SUFFIX}" \
    "${EXPERTGUESS_DIR}/benchmark/contig_specifications.csv" \
    "${UNINDEXED_FIXED_DIR}/scaffolds_pdb${SUFFIX}"
EXIT_CODE=$?

echo "========================================"
if [[ ${EXIT_CODE} -eq 0 ]]; then
    echo " Task finished successfully — $(date)"
else
    echo " Task FAILED with exit code ${EXIT_CODE} — $(date)"
fi
echo "========================================"

exit ${EXIT_CODE}
