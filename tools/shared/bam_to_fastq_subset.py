#!/usr/bin/env python3
from __future__ import annotations

import argparse
import pathlib

import pysam


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Write a FASTQ from primary mapped BAM reads, optionally read-limited.")
    parser.add_argument("--bam", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--max-reads", type=int, default=250000, help="Maximum reads to write; use 0 for all reads.")
    args = parser.parse_args(argv)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with pysam.AlignmentFile(args.bam) as bam, args.output.open("w") as out:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            seq = read.get_forward_sequence()
            quals = read.get_forward_qualities()
            if not seq or quals is None:
                continue
            qual_text = "".join(chr(q + 33) for q in quals)
            out.write(f"@{read.query_name}\n{seq}\n+\n{qual_text}\n")
            written += 1
            if args.max_reads > 0 and written >= args.max_reads:
                break
    print(f"wrote {written} reads to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
