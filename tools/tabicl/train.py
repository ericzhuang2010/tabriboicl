from __future__ import annotations

import argparse
import pathlib
import time

import pandas as pd

from common.feature_sets import select_feature_columns
from common.subsets import filter_read_groups, load_read_subset
from tools.shared.learned_candidate_models import (
    estimate_offset_prior,
    sample_read_groups,
    score_candidate_probabilities,
)

def prepare_features(frame: pd.DataFrame, feature_cols: list[str]) -> pd.DataFrame:
    out = frame[feature_cols].copy()
    for column in out.columns:
        if out[column].dtype == "object":
            out[column] = out[column].fillna("missing").astype(str)
        else:
            out[column] = out[column].fillna(-1)
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--train", type=pathlib.Path, required=True)
    parser.add_argument("--test", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--top-output", type=pathlib.Path)
    parser.add_argument("--n-estimators", type=int, default=4)
    parser.add_argument("--max-train-reads", type=int)
    parser.add_argument("--max-test-reads", type=int)
    parser.add_argument("--target-read-ids", type=pathlib.Path)
    parser.add_argument("--seed", type=int, default=20260517)
    parser.add_argument("--feature-set", choices=("full", "no_boundary", "minimal"), default="no_boundary")
    parser.add_argument("--tabicl-score-power", type=float, default=1.5)
    parser.add_argument("--offset-prior-power", type=float, default=1.0)
    parser.add_argument("--offset-prior-smoothing", type=float, default=1.0)
    parser.add_argument("--device", default="auto", help="TabICL device: auto, cpu, cuda, cuda:0, ...")
    parser.add_argument("--kv-cache", choices=("false", "true", "kv", "repr"), default="kv")
    args = parser.parse_args(argv)

    from tabicl import TabICLClassifier
    import torch

    start_time = time.perf_counter()
    if args.device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    else:
        device = args.device
    print(f"TabICL device: {device}")

    if args.kv_cache == "false":
        kv_cache_arg: bool | str = False
    elif args.kv_cache == "true":
        kv_cache_arg = True
    else:
        kv_cache_arg = args.kv_cache
    print(f"TabICL kv_cache: {args.kv_cache}")

    train = sample_read_groups(pd.read_parquet(args.train), args.max_train_reads, args.seed)
    test = pd.read_parquet(args.test)
    if args.target_read_ids:
        test = filter_read_groups(test, load_read_subset(args.target_read_ids))
    else:
        test = sample_read_groups(test, args.max_test_reads, args.seed + 1)
    if test.empty:
        raise ValueError("selected test read subset is empty")
    feature_cols = select_feature_columns(set(train.columns), set(test.columns), args.feature_set)
    if not feature_cols:
        raise ValueError("no configured feature columns found in train/test data")

    clf = TabICLClassifier(kv_cache=kv_cache_arg, n_estimators=args.n_estimators, device=device)
    clf.fit(prepare_features(train, feature_cols), train["label"].astype(int))

    out = test.copy()
    out["p_raw"] = clf.predict_proba(prepare_features(out, feature_cols))[:, 1]
    out["offset_prior"] = estimate_offset_prior(train, out, args.offset_prior_smoothing)
    out["p_offset"] = score_candidate_probabilities(
        out,
        raw_power=args.tabicl_score_power,
        prior_power=args.offset_prior_power,
    )

    pred = (
        out.sort_values(["read_id", "p_offset"], ascending=[True, False])
        .groupby("read_id", as_index=False)
        .head(1)
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out[["read_id", "candidate_offset", "p_raw", "offset_prior", "p_offset"]].to_parquet(args.output, index=False)
    if args.top_output:
        args.top_output.parent.mkdir(parents=True, exist_ok=True)
        pred[["read_id", "candidate_offset", "p_offset"]].to_parquet(args.top_output, index=False)
    elapsed = time.perf_counter() - start_time
    print(f"TabICL finished in {elapsed:.2f} seconds")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
