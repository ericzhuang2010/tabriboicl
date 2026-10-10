from __future__ import annotations

import argparse
import pathlib

import pandas as pd

from common.evaluation import evaluate_candidate_predictions


def mode_int(series: pd.Series) -> int:
    counts = series.astype(int).value_counts()
    return int(counts[counts.eq(counts.max())].index.min())


def choose_available_offsets(test: pd.DataFrame, target_by_read: dict[str, int]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for read_id, group in test.groupby("read_id", sort=False):
        candidates = sorted(group["candidate_offset"].astype(int).unique())
        # the default value candidates[0] should not be needed
        target = int(target_by_read.get(str(read_id), candidates[0]))
        if target in candidates:
            chosen = target
        else:
            chosen = min(candidates, key=lambda value: (abs(value - target), value))
        for candidate in candidates:
            rows.append(
                {
                    "read_id": read_id,
                    "candidate_offset": int(candidate),
                    "p_offset": 1.0 if int(candidate) == chosen else 0.0,
                }
            )
    return pd.DataFrame(rows)

def run_baselines(train: pd.DataFrame, test: pd.DataFrame, global_offset: int | None) -> dict[str, pd.DataFrame]:
    # Look only at training rows where label == 1
    # For the intended labeled benchmark data, each usable read_id should
    # have exactly one row where label=1
    positives = train[train["label"].eq(1)].copy()
    if positives.empty:
        raise ValueError("training split has no positive candidate rows")

   # Count which candidate_offset value appears most often.
    estimated_global = mode_int(positives["candidate_offset"])
    fixed_offset = estimated_global if global_offset is None else int(global_offset)
    # gets all the read ids in test
    read_ids = test["read_id"].astype(str).unique()
    # map each read_id to the fixed_offset
    fixed_targets = {read_id: fixed_offset for read_id in read_ids}

    # dict of read_length: most_common_true_offset_for_that_length
    length_offsets = positives.groupby("read_length")["candidate_offset"].agg(mode_int).to_dict()
    
    # if the same read_id appears with two different read_length values, only the 
    # first occurrence is kept.
    # In the intended data model, a read_id should have exactly one read_length. 
    # The same read appears multiple times only because it has multiple candidate_offset 
    # rows, not because its length changes.
    test_read_lengths = test[["read_id", "read_length"]].drop_duplicates("read_id")

    # a dictionary mapping each (test read id, read_length) to the offset predicted
    # for its read length:
    length_targets = {
        str(row.read_id): int(length_offsets.get(row.read_length, fixed_offset))
        for row in test_read_lengths.itertuples(index=False)
    }

    return {
        "fixed_offset": choose_available_offsets(test, fixed_targets),
        "read_length_offset": choose_available_offsets(test, length_targets),
    }


# For the intended labeled benchmark data, each usable read_id should
# have exactly one row where label=1
# each read_id should have the same read_length, but candidate_offset can be different
def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run deterministic p-site offset baselines.")
    parser.add_argument("--train", type=pathlib.Path, required=True)
    parser.add_argument("--test", type=pathlib.Path, required=True)
    parser.add_argument("--outdir", type=pathlib.Path, required=True)
    parser.add_argument(
        "--global-offset",
        default="auto",
        help="Use this fixed offset, or 'auto' to estimate the modal positive offset from training labels.",
    )
    parser.add_argument(
        "--models",
        nargs="+",
        default=["fixed_offset", "read_length_offset"],
        choices=["fixed_offset", "read_length_offset"],
        help="Deterministic offset baseline models to run.",
    )
    parser.add_argument("--split", default="")
    args = parser.parse_args(argv)

    global_offset = None if args.global_offset == "auto" else int(args.global_offset)
    train = pd.read_parquet(args.train)
    test = pd.read_parquet(args.test)
    args.outdir.mkdir(parents=True, exist_ok=True)

    baseline_predictions = run_baselines(train, test, global_offset)
    for model_name in args.models:
        predictions = baseline_predictions[model_name]
        model_dir = args.outdir / model_name
        model_dir.mkdir(parents=True, exist_ok=True)
        pred_path = model_dir / "candidate_predictions.parquet"
        read_path = model_dir / "read_predictions.parquet"
        predictions.to_parquet(pred_path, index=False)
        metrics, pred_read = evaluate_candidate_predictions(test, predictions)
        metrics.update({"model": model_name, "split": args.split, "predictions_path": str(pred_path)})
        pred_read.to_parquet(read_path, index=False)
        print(pd.DataFrame([metrics]).to_string(index=False))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
