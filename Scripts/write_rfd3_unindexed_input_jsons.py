#!/usr/bin/env python3
"""
write_rfd3_unindexed_input_jsons.py (Unindexed_Fixed variant)

Generates one RFD3 "unindexed" input JSON per MotifBench problem, e.g.:

{
  "motif": {
    "dialect": 2,
    "input": ".../ExpertGuess/benchmark/motif_pdbs_gly/21_1B73.pdb",
    "length": "125",
    "unindex": "A1-2,B1,C1-3",
    "select_fixed_atoms": {"A1-2,B1,C1-3": "ALL"},
    "select_unfixed_sequence": "C2",
    "is_non_loopy": true
  }
}

Unlike AutoContigmap/ExpertGuess (which specify an exact contig recipe:
segment order and gap sizes), this just marks the whole motif region as
"unindexed" and lets RFD3 decide placement itself.

Source of truth for everything here is ExpertGuess's already-prepared,
already-validated files -- nothing is re-derived from scratch:
  - motif_pdbs_gly/<problem>.pdb: the UNK->GLY relabeled motif PDBs (multi-chain
    structure preserved as-is, NOT collapsed onto a single chain).
  - input_jsons/<problem>.json: reused directly for `length` and
    `select_unfixed_sequence` (already correctly computed against
    motif_pdbs_gly's own numbering).

`unindex` is derived fresh by reading each chain's residue range straight out
of the motif_pdbs_gly PDB file, so it's correct regardless of how many chains
a given problem's motif spans.

Usage
-----
python write_rfd3_unindexed_input_jsons.py
"""

import json
from pathlib import Path

EXPERTGUESS_DIR = Path("/zfs/s01/z04/home/kleiter/Projects/MotifBench/RFD3/ExpertGuess/benchmark")
GLY_PDB_DIR = EXPERTGUESS_DIR / "motif_pdbs_gly"
EXISTING_JSON_DIR = EXPERTGUESS_DIR / "input_jsons"
OUT_DIR = Path(__file__).resolve().parent / "input_jsons"

MOTIFS = ["01_1LDB", "02_1ITU", "03_2CGA", "04_5WN9", "05_5ZE9", "06_6E6R", "07_6E6R",
          "08_7AD5", "09_7CG5", "10_7WRK", "11_3TQB", "12_4JHW", "13_4JHW", "14_5IUS",
          "15_7A8S", "16_7BNY", "17_7DGW", "18_7MQQ", "19_7MQQ", "20_7UWL", "21_1B73",
          "22_1BCF", "23_1MPY", "24_1QY3", "25_2RKX", "26_3B5V", "27_4XOJ", "28_5YUI",
          "29_6CPA", "30_7UWL"]


def chain_ranges_from_pdb(pdb_path):
    """Read (chain, resnum) pairs from ATOM records, return per-chain contiguous
    ranges in first-seen chain order, e.g. [("A", 1, 2), ("B", 1, 1), ("C", 1, 3)]."""
    seen = {}  # chain -> sorted set of resnums
    order = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            chain = line[21]
            resnum = int(line[22:26])
            if chain not in seen:
                seen[chain] = set()
                order.append(chain)
            seen[chain].add(resnum)

    ranges = []
    for chain in order:
        resnums = sorted(seen[chain])
        first, last = resnums[0], resnums[-1]
        if list(range(first, last + 1)) != resnums:
            raise ValueError(f"{pdb_path}: chain {chain} is not contiguous ({resnums})")
        ranges.append((chain, first, last))
    return ranges


def format_unindex(ranges):
    parts = []
    for chain, first, last in ranges:
        parts.append(f"{chain}{first}" if first == last else f"{chain}{first}-{last}")
    return ",".join(parts)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    for problem in MOTIFS:
        gly_pdb = GLY_PDB_DIR / f"{problem}.pdb"
        existing_json_path = EXISTING_JSON_DIR / f"{problem}.json"

        with open(existing_json_path) as f:
            existing = json.load(f)["motif"]
        length = existing["length"]
        select_unfixed_sequence = existing.get("select_unfixed_sequence")

        ranges = chain_ranges_from_pdb(gly_pdb)
        unindex = format_unindex(ranges)

        motif_spec = {
            "dialect": 2,
            "input": str(gly_pdb.resolve()),
            "length": str(length),
            "unindex": unindex,
            "select_fixed_atoms": {unindex: "ALL"},
            "is_non_loopy": True,
        }
        if select_unfixed_sequence:
            motif_spec["select_unfixed_sequence"] = select_unfixed_sequence

        payload = {"motif": motif_spec}

        out_path = OUT_DIR / f"{problem}.json"
        with open(out_path, "w") as f:
            json.dump(payload, f, indent=2)
        print(f"wrote {out_path}  (unindex={unindex}, unfixed_seq={select_unfixed_sequence})")


if __name__ == "__main__":
    main()
