#!/usr/bin/env python3
"""
convert_rfd3_outputs.py (Unindexed variant)

Same purpose as AutoContigmap's convert_rfd3_outputs.py (RFD3's .cif.gz+.json
-> Scaffold-Lab's expected .pdb + scaffold_info.csv), but Unindexed mode's
JSON metadata doesn't expose a per-design contig (`specification.extra.sampled_contig`
is just the total length, since "unindex" mode doesn't build an explicit contig
recipe). Instead, `diffused_index_map` directly maps each INPUT motif residue
(e.g. "A1".."A24", positions in the input motif_pdbs_gly file) to its OUTPUT
residue number in the generated 100-residue chain -- so the set of output
values IS the motif placement for that design.

BUG HISTORY (fixed 2026-09-04)
------------------------------
The original version of this script (still present, unfixed, in
../../Unindexed/benchmark/convert_rfd3_outputs.py) assigned chain labels to the
reconstructed motif segments POSITIONALLY -- the k-th contiguous run of motif
residues found in the output got chains[k], the k-th chain of the problem's
contig. That silently assumes RFD3 emits the motif segments in contig order.
"unindex" mode makes no such promise: it places each segment freely and permutes
them in ~2/3 of designs (per-problem: 8/79 for 11_3TQB up to 100/100 for
23_1MPY, which never once came out in contig order).

Consequences of the positional assignment:
  * Wrong labels are only *sometimes* fatal. write_motifInfo_from_scaffoldInfo.py
    expands "B:17" as "chain B's first 17 motif residues", so a swap between
    chains of unequal length asks for more residues than the chain has, the
    contig comes up short, and Scaffold-Lab dies with
    "IndexError: boolean index did not match indexed array along dimension 0;
     dimension is 100 but corresponding boolean dimension is 90"
    (this is what killed 12_4JHW in the rep1 evaluation, job 526691_11).
  * Otherwise it is SILENT: the design is evaluated with its motif residues
    paired against the wrong reference residues, i.e. a meaningless motif-RMSD.
    The old Unindexed run never crashed only because it used single-collapsed-
    chain motif PDBs (every segment labeled "A"), so no label could overflow --
    63% of its kept multi-chain designs were scored against the wrong reference
    positions instead.

The fix: take each segment's chain from the diffused_index_map KEYS (the input
residue labels, which say exactly which chain each motif residue came from)
rather than from run order, and validate that each run is one whole chain in
ascending source order -- the only shape write_motifInfo can expand correctly.

Usage
-----
python convert_rfd3_outputs.py <scaffolds_dir> <contig_specifications.csv> <out_dir> [motif_name]
"""
import csv
import gzip
import json
import re
import shutil
import sys
import tempfile
from pathlib import Path

from Bio.PDB import MMCIFParser, PDBIO


def parse_contig_specifications(contig_file):
    chain_order = {}
    with open(contig_file, newline="") as f:
        for row in csv.DictReader(f):
            problem = row["problem"]
            contig = row["contig"]
            chains = [seg[0] for seg in contig.split(";") if seg and seg[0].isalpha()]
            chain_order[problem] = chains
    return chain_order


class SegmentCountMismatch(Exception):
    """Raised when the reconstructed motif segments aren't a shape the
    downstream MotifBench scripts can expand: RFD3 split a segment into extra
    pieces, merged two, interleaved residues from different chains into one
    contiguous run, or laid a chain's residues out non-monotonically.
    check_segment_validity.py hard-crashes the ENTIRE motif's evaluation on any
    single such sample, so these must be excluded during conversion rather than
    passed through."""


def placements_from_index_map(diffused_index_map, total_length, chains):
    """diffused_index_map maps each INPUT motif residue label ("B17") to its
    OUTPUT position ("A67") in the generated chain. Group the output positions
    into contiguous runs (= the placed motif segments), label each run with the
    chain its residues actually came from, and fill the gaps between runs with
    scaffold lengths.

    The chain label comes from the map's keys, never from run order: "unindex"
    mode permutes the segments freely, so the k-th run is NOT the k-th chain of
    the contig (see BUG HISTORY at the top of this file).

    Raises SegmentCountMismatch unless every run is exactly one chain's motif
    residues, complete and in ascending residue order. That is the only layout
    write_motifInfo_from_scaffoldInfo.py can expand correctly: it turns "B:17"
    into chain B's first 17 motif residues in ascending order, so a run that is
    partial, reordered, or mixes chains would be paired against the wrong
    reference residues.
    """
    if len(set(chains)) != len(chains):
        raise ValueError(
            f"expected one motif segment per distinct chain, got {chains} -- "
            "the legacy single-collapsed-chain convention is not supported"
        )

    # (output position, source chain, source residue number) per motif residue
    residues = sorted(
        (int(out_label[1:]), in_label[0], int(in_label[1:]))
        for in_label, out_label in diffused_index_map.items()
    )

    # each chain's motif residues in the order write_motifInfo hands them out
    chain_residues = {}
    for _, chain, resnum in residues:
        chain_residues.setdefault(chain, []).append(resnum)
    for nums in chain_residues.values():
        nums.sort()

    placements = []
    found = []
    prev_end = 0  # 1-indexed; "prev_end" = last position already accounted for
    i = 0
    n = len(residues)
    while i < n:
        j = i + 1
        while j < n and residues[j][0] == residues[j - 1][0] + 1:
            j += 1
        run = residues[i:j]

        gap = run[0][0] - prev_end - 1
        if gap > 0:
            placements.append(str(gap))

        chain = run[0][1]
        if any(r[1] != chain for r in run):
            raise SegmentCountMismatch(
                f"contiguous output run at {run[0][0]}-{run[-1][0]} mixes input "
                f"chains {sorted({r[1] for r in run})}"
            )
        if [r[2] for r in run] != chain_residues[chain]:
            raise SegmentCountMismatch(
                f"run for chain {chain} is not that chain's motif residues in "
                f"ascending order (got {[r[2] for r in run]}, "
                f"expected {chain_residues[chain]})"
            )

        placements.append(f"{chain}:{len(run)}")
        found.append(chain)
        prev_end = run[-1][0]
        i = j

    if sorted(found) != sorted(chains):
        raise SegmentCountMismatch(
            f"found motif segments for chains {found}, expected {chains}"
        )

    trailing_gap = total_length - prev_end
    if trailing_gap > 0:
        placements.append(str(trailing_gap))

    return "/".join(placements)


def convert_one_motif(scaffolds_dir, motif_name, chains, out_dir):
    src_dir = scaffolds_dir / motif_name
    dst_dir = out_dir / motif_name
    dst_dir.mkdir(parents=True, exist_ok=True)

    parser = MMCIFParser(QUIET=True)
    io = PDBIO()
    rows = []

    cif_files = sorted(src_dir.glob(f"{motif_name}_motif_*_model_0.cif.gz"))
    if not cif_files:
        print(f"warning: no .cif.gz files found for {motif_name} in {src_dir}", file=sys.stderr)
        return

    n_skipped = 0
    for cif_gz in cif_files:
        m = re.match(rf"{re.escape(motif_name)}_motif_(\d+)_model_0\.cif\.gz", cif_gz.name)
        sample_idx = int(m.group(1))
        json_path = cif_gz.with_suffix("").with_suffix(".json")

        with open(json_path) as f:
            meta = json.load(f)
        total_length = int(meta["specification"]["length"])
        diffused_index_map = meta["diffused_index_map"]
        try:
            motif_placements = placements_from_index_map(diffused_index_map, total_length, chains)
        except SegmentCountMismatch:
            n_skipped += 1
            continue

        with gzip.open(cif_gz, "rt") as f_in, tempfile.NamedTemporaryFile(
            mode="w", suffix=".cif", delete=False
        ) as f_tmp:
            shutil.copyfileobj(f_in, f_tmp)
            tmp_path = f_tmp.name

        structure = parser.get_structure(motif_name, tmp_path)
        Path(tmp_path).unlink()

        out_pdb = dst_dir / f"{motif_name}_{sample_idx}.pdb"
        io.set_structure(structure)
        io.save(str(out_pdb))

        rows.append([sample_idx, motif_placements])

    rows.sort(key=lambda r: r[0])
    with open(dst_dir / "scaffold_info.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sample_num", "motif_placements"])
        writer.writerows(rows)

    print(f"{motif_name}: converted {len(rows)} designs, skipped {n_skipped} (segment-count mismatch) -> {dst_dir}")


def main():
    if len(sys.argv) not in (4, 5):
        print(f"Usage: {sys.argv[0]} <scaffolds_dir> <contig_specifications.csv> <out_dir> [motif_name]", file=sys.stderr)
        sys.exit(1)

    scaffolds_dir = Path(sys.argv[1])
    contig_file = sys.argv[2]
    out_dir = Path(sys.argv[3])
    only_motif = sys.argv[4] if len(sys.argv) == 5 else None

    chain_order = parse_contig_specifications(contig_file)

    motifs = [only_motif] if only_motif else sorted(d.name for d in scaffolds_dir.iterdir() if d.is_dir())
    for motif_name in motifs:
        chains = chain_order.get(motif_name, [])
        convert_one_motif(scaffolds_dir, motif_name, chains, out_dir)


if __name__ == "__main__":
    main()
