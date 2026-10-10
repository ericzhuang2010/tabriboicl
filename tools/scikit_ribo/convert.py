#!/usr/bin/env python3
from __future__ import annotations

import argparse
import pathlib
from collections import defaultdict

import pandas as pd


def normalize_scores(scores: list[float], pseudocount: float) -> list[float]:
    adjusted = [max(float(score), 0.0) + pseudocount for score in scores]
    total = sum(adjusted)
    if total <= 0:
        return [1.0 / len(scores)] * len(scores)
    return [score / total for score in adjusted]


def target_positions(sample_frame: pd.DataFrame, window_nt: int) -> dict[str, set[int]]:
    targets: dict[str, set[int]] = defaultdict(set)
    for row in sample_frame[["contig", "candidate_ref_pos"]].drop_duplicates().itertuples(index=False):
        chrom = str(getattr(row, "contig"))
        pos = int(getattr(row, "candidate_ref_pos"))
        if window_nt <= 0:
            targets[chrom].add(pos)
        else:
            targets[chrom].update(range(pos - window_nt, pos + window_nt + 1))
    return targets


def profile_scores(
    path: pathlib.Path,
    site: str,
    targets: dict[str, set[int]] | None = None,
    chunksize: int = 1_000_000,
) -> dict[tuple[str, int], float]:
    scores: dict[tuple[str, int], float] = defaultdict(float)
    if not path.exists() or path.stat().st_size == 0:
        return {}

    required = ["chrom", "start", "end", "gene_strand", "ribosome_count"]
    header = pd.read_csv(path, sep="\t", nrows=0)
    missing = sorted(set(required) - set(header.columns))
    if missing:
        raise ValueError(f"{path} is missing required scikit-ribo profile columns: {missing}")

    target_chroms = set(targets) if targets is not None else None
    for frame in pd.read_csv(path, sep="\t", usecols=required, chunksize=chunksize):
        frame = frame[frame["ribosome_count"].astype(float) > 0]
        if frame.empty:
            continue
        for row in frame.itertuples(index=False):
            count = float(getattr(row, "ribosome_count"))
            chrom = str(getattr(row, "chrom"))
            if target_chroms is not None and chrom not in target_chroms:
                continue
            start = int(getattr(row, "start"))
            end = int(getattr(row, "end"))
            strand = str(getattr(row, "gene_strand"))
            coord = start if strand == "+" else end - 1
            if site == "psite":
                coord = coord - 3 if strand == "+" else coord + 3
            if targets is not None and coord not in targets.get(chrom, set()):
                continue
            scores[(chrom, coord)] += count
    return dict(scores)

def score_candidate(
    scores: dict[tuple[str, int], float],
    contig: str,
    pos: int,
    window_nt: int,
) -> float:
    if window_nt <= 0:
        return float(scores.get((contig, pos), 0.0))
    return max(float(scores.get((contig, nearby), 0.0)) for nearby in range(pos - window_nt, pos + window_nt + 1))


def predictions_for_test(args: argparse.Namespace, test: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    profile_cache: dict[str, dict[tuple[str, int], float]] = {}

    for sample_id, sample_frame in test.groupby("sample_id", sort=False):
        sample = str(sample_id)
        profile_path = args.baseline_root / sample / "riboseq_input.txt"
        targets = target_positions(sample_frame, args.window_nt)
        scores = profile_cache.setdefault(
            sample, profile_scores(profile_path, args.site, targets=targets, chunksize=args.profile_chunksize)
        )
        for read_id, read_frame in sample_frame.groupby("read_id", sort=False):
            candidate_rows = list(read_frame.itertuples(index=False))
            raw_scores = [
                score_candidate(
                    scores,
                    contig=str(getattr(row, "contig")),
                    pos=int(getattr(row, "candidate_ref_pos")),
                    window_nt=args.window_nt,
                )
                for row in candidate_rows
            ]
            probs = normalize_scores(raw_scores, args.pseudocount)
            rows.extend(
                {
                    "read_id": str(read_id),
                    "candidate_offset": int(getattr(row, "candidate_offset")),
                    "p_offset": prob,
                }
                for row, prob in zip(candidate_rows, probs)
            )

    if not rows:
        return pd.DataFrame(columns=["read_id", "candidate_offset", "p_offset"])
    return pd.DataFrame.from_records(rows, columns=["read_id", "candidate_offset", "p_offset"]).drop_duplicates(
        ["read_id", "candidate_offset"]
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Convert scikit-ribo A-site profiles into B1 candidate predictions.")
    parser.add_argument("--test", type=pathlib.Path, required=True)
    parser.add_argument("--baseline-root", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--site", choices=("psite", "asite"), default="psite")
    parser.add_argument("--window-nt", type=int, default=0)
    parser.add_argument("--pseudocount", type=float, default=1e-9)
    parser.add_argument("--profile-chunksize", type=int, default=1_000_000)
    args = parser.parse_args(argv)

    test = pd.read_parquet(
        args.test,
        columns=["read_id", "sample_id", "contig", "candidate_ref_pos", "candidate_offset"],
    )
    predictions = predictions_for_test(args, test)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    predictions.to_parquet(args.output, index=False)
    print(f"wrote {len(predictions)} scikit-ribo predictions to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
