#!/usr/bin/env python3
"""
Append newly called Clippy peak BEDs to RBPeekSamplesheet.tsv.

Run this ON THE HPC, from the RBPeek directory. It discovers the peak files rather than
taking hardcoded names, because the exact Clippy output filenames depend on the parameters
it was run with.

Two categories are skipped by default, both for reasons that already cost us a run:

  TRA2A  already in the panel as HepG2-TRA2A-eCLIP and K562-TRA2A-eCLIP, called from the
         same underlying experiments. Adding a second copy would create a column pair that
         correlates by construction and reads as a reproducible partnership.
  _Mm    mouse. The inference BED is hg38, and mouse shares chromosome NAMES with human but
         not coordinates, so bedtools would report coincidental overlaps as real signal.

Paths are written relative to --xldir, since that is how intersect_inference_bed.py resolves
the `file` column. The script is idempotent: rows already present are reported and skipped,
so it is safe to re-run after calling more peaks.
"""

import argparse
import csv
import os
import re
import sys
from pathlib import Path

DEFAULT_SKIP = "TRA2A,_Mm"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--peaks-dir", type=Path, default=Path("."), help="Directory to scan for *_Peaks.bed (default: cwd)")
    p.add_argument("--samplesheet", type=Path, default=Path("RBPeekSamplesheet.tsv"))
    p.add_argument("--xldir", type=Path, default=Path("../CLIP"), help="The -x value the run uses; paths are made relative to it")
    p.add_argument("--skip", default=DEFAULT_SKIP, help=f"Comma-separated substrings to skip (default: {DEFAULT_SKIP})")
    p.add_argument("--glob", default="*_Peaks.bed", help="Filename pattern to scan for (default: *_Peaks.bed)")
    p.add_argument("--recursive", action="store_true", help="Also scan subdirectories (e.g. clippy_peaks/)")
    p.add_argument("--dry-run", action="store_true", help="Show what would be added and exit without writing")
    return p.parse_args()


def group_from_filename(name: str) -> str:
    """
    Clippy names outputs <prefix>_rollmean<N>_..._Peaks.bed, so the prefix is the group.
    Falls back to stripping the known suffixes when the name does not match.
    """
    stem = re.sub(r"\.bed$", "", name)
    for marker in ("_rollmean", "_roll", "_Summits", "_Peaks"):
        i = stem.find(marker)
        if i > 0:
            return stem[:i]
    return stem


def main():
    args = parse_args()
    if not args.samplesheet.is_file():
        sys.exit(f"samplesheet not found: {args.samplesheet} (run this from the RBPeek directory)")

    skips = [s.strip() for s in args.skip.split(",") if s.strip()]
    it = args.peaks_dir.rglob(args.glob) if args.recursive else args.peaks_dir.glob(args.glob)
    found = sorted(p for p in it if p.is_file())
    if not found:
        sys.exit(
            f"no files matching {args.glob!r} under {args.peaks_dir}"
            + ("" if args.recursive else " (try --recursive if peaks are in a subdirectory)")
        )

    with open(args.samplesheet, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    existing_files = {r["file"] for r in rows}
    existing_groups = {r["group"] for r in rows}

    to_add, skipped = [], []
    for path in found:
        name = path.name
        hit = next((s for s in skips if s in name), None)
        if hit:
            skipped.append((name, f"matches skip pattern {hit!r}"))
            continue
        rel = os.path.relpath(path.resolve(), args.xldir.resolve())
        group = group_from_filename(name)
        if rel in existing_files:
            skipped.append((name, "path already in samplesheet"))
            continue
        if group in existing_groups:
            skipped.append((name, f"group {group!r} already in samplesheet"))
            continue
        if not path.stat().st_size:
            skipped.append((name, "file is empty"))
            continue
        to_add.append((rel, group))
        existing_groups.add(group)
        existing_files.add(rel)

    print(f"scanned {len(found)} file(s) under {args.peaks_dir}")
    if skipped:
        print("\nskipped:")
        for name, why in skipped:
            print(f"  {name}\n      {why}")
    if not to_add:
        print("\nnothing to add - samplesheet already up to date.")
        return

    print(f"\nto add ({len(to_add)}):")
    for rel, group in to_add:
        print(f"  {group:28} <- {rel}")

    if args.dry_run:
        print("\n--dry-run set, samplesheet not modified.")
        return

    with open(args.samplesheet, "a", newline="") as fh:
        for rel, group in to_add:
            fh.write(f"{rel}\t{group}\n")
    total = len(rows) + len(to_add)
    print(f"\nappended to {args.samplesheet}; now {total} panel entries (was {len(rows)}).")
    print("Check with: git diff RBPeekSamplesheet.tsv")


if __name__ == "__main__":
    main()
