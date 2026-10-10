#!/usr/bin/env python3
from __future__ import annotations

import argparse
import pathlib
import re

import pandas as pd


OFFSET_COLUMNS = ["sample_id", "read_length", "five_prime_frame", "asite_offset", "psite_offset"]
PREDICTION_COLUMNS = ["read_id", "candidate_offset", "p_offset"]
OFFSET_RE = re.compile(r"^[+-]?\d+$")


def parse_native_offset(raw: object, ambiguous_policy: str = "drop") -> int | None:
    value = str(raw).strip()
    if not value or value.upper() == "NA":
        return None
    if "/" in value:
        if ambiguous_policy == "drop":
            return None
        raise ValueError(f"ambiguous Asite IP offset: {value}")
    if not OFFSET_RE.match(value):
        if ambiguous_policy == "drop":
            return None
        raise ValueError(f"unrecognized Asite IP offset: {value}")
    return int(value)


def parse_results_ip_algorithm(
    path: pathlib.Path,
    sample_id: str,
    ambiguous_policy: str = "drop",
) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    in_offset_table = False
    saw_offset_header = False

    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if not in_offset_table:
                if fields[:4] == ["Frag size", "Frame_0", "Frame_1", "Frame_2"]:
                    in_offset_table = True
                    saw_offset_header = True
                continue

            if not line.strip():
                break
            if len(fields) < 4:
                break
            try:
                read_length = int(fields[0])
            except ValueError:
                break

            for frame, raw_offset in enumerate(fields[1:4]):
                asite_offset = parse_native_offset(raw_offset, ambiguous_policy=ambiguous_policy)
                if asite_offset is None:
                    continue
                records.append(
                    {
                        "sample_id": sample_id,
                        "read_length": read_length,
                        "five_prime_frame": frame,
                        "asite_offset": asite_offset,
                        "psite_offset": asite_offset - 3,
                    }
                )

    if not saw_offset_header:
        raise ValueError(f"{path} does not contain a 'Frag size/Frame_0/Frame_1/Frame_2' offset table")

    return pd.DataFrame.from_records(records, columns=OFFSET_COLUMNS)


def load_offset_tables(
    results: pathlib.Path | None = None,
    sample_id: str | None = None,
    baseline_root: pathlib.Path | None = None,
    ambiguous_policy: str = "drop",
) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []

    if results is not None:
        if sample_id is None:
            raise ValueError("--sample-id is required when --results is used")
        frames.append(parse_results_ip_algorithm(results, sample_id=sample_id, ambiguous_policy=ambiguous_policy))

    if baseline_root is not None:
        for path in sorted(baseline_root.glob("*/output/Results_IP_algorithm.tab")):
            inferred_sample_id = path.parent.parent.name
            frames.append(
                parse_results_ip_algorithm(
                    path,
                    sample_id=inferred_sample_id,
                    ambiguous_policy=ambiguous_policy,
                )
            )

    if not frames:
        return pd.DataFrame(columns=OFFSET_COLUMNS)
    return pd.concat(frames, ignore_index=True)


def convert_offsets_to_predictions(
    test: pd.DataFrame,
    offsets: pd.DataFrame,
    offset_min: int = 8,
    offset_max: int = 20,
) -> pd.DataFrame:
    if offsets.empty or test.empty:
        return pd.DataFrame(columns=PREDICTION_COLUMNS)

    offset_candidates = offsets.copy()
    offset_candidates["candidate_offset"] = offset_candidates["psite_offset"].astype(int)
    offset_candidates = offset_candidates[
        offset_candidates["candidate_offset"].between(offset_min, offset_max)
    ].copy()
    if offset_candidates.empty:
        return pd.DataFrame(columns=PREDICTION_COLUMNS)

    join_columns = ["sample_id", "read_length", "five_prime_frame", "candidate_offset"]
    test_candidates = test[join_columns + ["read_id"]].copy()
    for column in ["read_length", "five_prime_frame", "candidate_offset"]:
        test_candidates[column] = pd.to_numeric(test_candidates[column], errors="raise").astype(int)
        offset_candidates[column] = pd.to_numeric(offset_candidates[column], errors="raise").astype(int)

    predictions = test_candidates.merge(
        offset_candidates[join_columns].drop_duplicates(),
        on=join_columns,
        how="inner",
        sort=False,
    )
    if predictions.empty:
        return pd.DataFrame(columns=PREDICTION_COLUMNS)

    predictions = predictions[["read_id", "candidate_offset"]].drop_duplicates()
    predictions["p_offset"] = 1.0
    return predictions[PREDICTION_COLUMNS].reset_index(drop=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Convert Asite IP A-site offsets to P-site candidate predictions.")
    parser.add_argument("--test", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    source_group = parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument("--results", type=pathlib.Path)
    source_group.add_argument("--baseline-root", type=pathlib.Path)
    parser.add_argument("--sample-id")
    parser.add_argument("--offset-min", type=int, default=8)
    parser.add_argument("--offset-max", type=int, default=20)
    parser.add_argument("--ambiguous-policy", choices=("drop", "error"), default="drop")
    args = parser.parse_args(argv)

    if args.results is not None and args.sample_id is None:
        parser.error("--sample-id is required when --results is used")

    offsets = load_offset_tables(
        results=args.results,
        sample_id=args.sample_id,
        baseline_root=args.baseline_root,
        ambiguous_policy=args.ambiguous_policy,
    )
    test = pd.read_parquet(
        args.test,
        columns=["read_id", "sample_id", "read_length", "five_prime_frame", "candidate_offset"],
    )
    predictions = convert_offsets_to_predictions(
        test,
        offsets,
        offset_min=args.offset_min,
        offset_max=args.offset_max,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    predictions.to_parquet(args.output, index=False)
    print(f"wrote {len(predictions)} Asite IP predictions to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
