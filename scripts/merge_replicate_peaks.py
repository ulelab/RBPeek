#!/usr/bin/env python3
"""
Merge replicate peak BEDs by summing scores at IDENTICAL intervals.

Implements the pipeline as specified:

    cat *.bed | sort -k1,1 -k2,2n -k3,3n -k6,6 \\
      | bedtools groupby -g 1,2,3,6 -c 5 -o sum

then reformats groupby's output (chrom, start, end, strand, sum) back to BED6 by inserting
a name column and returning strand to column 6.

Two deliberate departures from the pasted command:

  * the stray `cat merged_Prpf8.bed` between the sort and the groupby is dropped - it
    discards stdin and reads a file instead, so the sorted stream never reaches groupby.
  * `-k3,3n` is added to the sort. groupby only merges ADJACENT rows sharing the group key,
    and the key includes `end` (field 3). Without sorting on it, two rows with the same
    chrom/start but different ends can interleave and split a group that should be one.

What this does and does not do
------------------------------
Grouping is on EXACT chrom/start/end/strand, so only byte-identical intervals combine. That
is the right operation for 1 nt crosslink sites. For peaks it collapses very little: on the
four THRAP3 replicates, exact matching found 2+ replicate support at 4.5% of intervals where
overlap-based merging found 36%. Use build_inference_bed_from_peaks.py if you want
overlap-based reproducibility instead; this script is for pooling replicate signal.

Chromosome naming is left exactly as found. intersect_inference_bed.py harmonises panel
naming against the inference BED at run time.
"""

import argparse
import os
import subprocess
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--inputs", nargs="+", required=True, help="Replicate peak BEDs to merge")
    p.add_argument("-l", "--label", required=True, help="Group label, used as the BED name column")
    p.add_argument("-o", "--out", required=True, help="Output BED6")
    return p.parse_args()


def require(tool):
    from shutil import which
    if which(tool) is None:
        sys.exit("%s not found on PATH. Try: conda activate rbpeek" % tool)


def main():
    args = parse_args()
    require("bedtools")
    require("sort")
    missing = [f for f in args.inputs if not os.path.isfile(f)]
    if missing:
        sys.exit("not found: %s" % ", ".join(missing))

    n_in = 0
    for f in args.inputs:
        with open(f) as fh:
            n_in += sum(1 for l in fh if l.strip())

    d = os.path.dirname(args.out)
    if d:
        os.makedirs(d, exist_ok=True)

    cmd = (
        "cat " + " ".join('"%s"' % f for f in args.inputs) +
        " | sort -k1,1 -k2,2n -k3,3n -k6,6"
        " | bedtools groupby -g 1,2,3,6 -c 5 -o sum"
        " | awk -v OFS='\\t' -v L=\"%s\" '{print $1,$2,$3,L\"_\"NR,$5,$4}'" % args.label +
        ' > "%s"' % args.out
    )
    r = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, universal_newlines=True)
    if r.returncode != 0:
        sys.exit("merge failed:\n%s" % r.stderr[:600])

    with open(args.out) as fh:
        n_out = sum(1 for l in fh if l.strip())
    bad = 0
    with open(args.out) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) != 6 or c[5] not in ("+", "-"):
                bad += 1
    print("  %-22s %d replicate(s), %7d rows in -> %7d rows out  (%.1f%% collapsed)%s"
          % (args.label, len(args.inputs), n_in, n_out,
             100.0 * (1 - n_out / float(n_in)) if n_in else 0.0,
             "" if not bad else "   MALFORMED ROWS: %d" % bad))


if __name__ == "__main__":
    main()
