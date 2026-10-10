#!/usr/bin/env python3
from __future__ import annotations

import argparse
import pathlib

import pysam
from Bio import SeqIO


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Filter BAM alignments to references present in a FASTA file.")
    parser.add_argument("--bam", type=pathlib.Path, required=True)
    parser.add_argument("--fasta", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args(argv)

    keep_refs = {record.id for record in SeqIO.parse(args.fasta, "fasta")}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    kept = 0
    total = 0
    with pysam.AlignmentFile(args.bam) as src, pysam.AlignmentFile(args.output, "wb", header=src.header) as dst:
        for read in src.fetch(until_eof=True):
            total += 1
            if read.reference_name in keep_refs:
                dst.write(read)
                kept += 1
    print(f"kept {kept} of {total} alignments in {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
