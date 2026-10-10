#!/usr/bin/env python3
from __future__ import annotations

import argparse
import pathlib

import pandas as pd


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Convert riboWaltz length offsets to B1 predictions.")
    parser.add_argument("--test", type=pathlib.Path, required=True)
    parser.add_argument("--offset-dir", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--offset-min", type=int, default=8)
    parser.add_argument("--offset-max", type=int, default=20)
    args = parser.parse_args(argv)

    offset_frames = []
    for path in sorted(args.offset_dir.glob("*.ribowaltz_offsets.tsv")):
        frame = pd.read_csv(path, sep="\t")
        if frame.empty:
            continue
        offset_frames.append(
            frame[["sample", "length", "corrected_offset_from_5"]].rename(
                columns={
                    "sample": "sample_id",
                    "length": "read_length",
                    "corrected_offset_from_5": "candidate_offset",
                }
            )
        )
    if not offset_frames:
        predictions = pd.DataFrame(columns=["read_id", "candidate_offset", "p_offset"])
    else:
        offsets = pd.concat(offset_frames, ignore_index=True)
        offsets["candidate_offset"] = pd.to_numeric(offsets["candidate_offset"], errors="coerce")
        offsets = offsets.dropna(subset=["candidate_offset"]).copy()
        offsets["candidate_offset"] = offsets["candidate_offset"].round().astype(int)
        offsets = offsets[
            offsets["candidate_offset"].between(args.offset_min, args.offset_max)
        ].copy()

        test = pd.read_parquet(args.test, columns=["read_id", "sample_id", "read_length"])
        reads = test.drop_duplicates("read_id")
        predictions = reads.merge(offsets, on=["sample_id", "read_length"], how="inner")
        predictions = predictions[["read_id", "candidate_offset"]].drop_duplicates()
        predictions["p_offset"] = 1.0

    args.output.parent.mkdir(parents=True, exist_ok=True)
    predictions.to_parquet(args.output, index=False)
    print(f"wrote {len(predictions)} riboWaltz predictions to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
