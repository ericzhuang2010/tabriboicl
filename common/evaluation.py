from __future__ import annotations

import argparse
from collections.abc import Iterable
import pathlib

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from common.metrics import binary_metrics, offset_metrics
from common.subsets import load_read_subset


def table_sep(path: pathlib.Path) -> str:
    return "," if path.suffix == ".csv" else "	"


def normalize_prediction_columns(predictions: pd.DataFrame) -> pd.DataFrame:
    out = predictions.copy()
    if "pred_offset" in out.columns and "candidate_offset" not in out.columns:
        out = out.rename(columns={"pred_offset": "candidate_offset"})
    required = {"read_id", "candidate_offset"}
    missing = sorted(required - set(out.columns))
    if missing:
        raise ValueError(f"prediction table is missing required columns: {missing}")
    if "p_offset" not in out.columns:
        out["p_offset"] = 1.0
    return out[["read_id", "candidate_offset", "p_offset"]].copy()


def _prediction_columns(path: pathlib.Path) -> list[str]:
    names = set(pq.ParquetFile(path).schema_arrow.names)
    columns = ["read_id"]
    if "candidate_offset" in names:
        columns.append("candidate_offset")
    elif "pred_offset" in names:
        columns.append("pred_offset")
    if "p_offset" in names:
        columns.append("p_offset")
    return columns


def _empty_predictions() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "read_id": pd.Series(dtype="object"),
            "candidate_offset": pd.Series(dtype="int64"),
            "p_offset": pd.Series(dtype="float64"),
        }
    )


def _parquet_row_count(path: pathlib.Path) -> int:
    return int(pq.ParquetFile(path).metadata.num_rows)


def _read_prediction_parquet_for_read_ids(
    path: pathlib.Path,
    read_ids: Iterable[object],
    *,
    batch_size: int,
) -> pd.DataFrame:
    target_ids = _unique_read_ids(read_ids)
    if not target_ids:
        return _empty_predictions()

    value_set = pa.array(target_ids, type=pa.string())
    parquet = pq.ParquetFile(path)
    columns = _prediction_columns(path)
    frames: list[pd.DataFrame] = []
    for batch in parquet.iter_batches(batch_size=batch_size, columns=columns):
        read_id_column = batch.column(batch.schema.get_field_index("read_id"))
        if not pa.types.is_string(read_id_column.type):
            read_id_column = pc.cast(read_id_column, pa.string())
        mask = pc.is_in(read_id_column, value_set=value_set)
        if not bool(pc.any(mask).as_py()):
            continue
        frame = batch.filter(mask).to_pandas()
        frame["read_id"] = frame["read_id"].astype(str)
        frames.append(frame)

    if not frames:
        return _empty_predictions()
    return normalize_prediction_columns(pd.concat(frames, ignore_index=True, sort=False))


def read_prediction_table(
    path: pathlib.Path,
    *,
    filter_read_ids: Iterable[object] | None = None,
    stream_threshold: int = 5_000_000,
    batch_size: int = 100_000,
) -> pd.DataFrame:
    path = pathlib.Path(path)
    if filter_read_ids is not None and path.suffix == ".parquet" and _parquet_row_count(path) > stream_threshold:
        return _read_prediction_parquet_for_read_ids(path, filter_read_ids, batch_size=batch_size)

    if path.suffix == ".parquet":
        predictions = pd.read_parquet(path, columns=_prediction_columns(path))
    else:
        predictions = pd.read_csv(path, sep=table_sep(path))
    predictions = normalize_prediction_columns(predictions)
    if filter_read_ids is not None:
        target_ids = set(_unique_read_ids(filter_read_ids))
        predictions = predictions[predictions["read_id"].astype(str).isin(target_ids)].copy()
        predictions["read_id"] = predictions["read_id"].astype(str)
    return predictions


def _require_test_columns(test: pd.DataFrame) -> None:
    required = {"read_id", "candidate_offset", "label", "candidate_frame"}
    missing = sorted(required - set(test.columns))
    if missing:
        raise ValueError(f"test table is missing required columns: {missing}")


def _unique_read_ids(read_ids: Iterable[object]) -> list[str]:
    return list(dict.fromkeys(str(read_id) for read_id in read_ids))


def _true_rows(test: pd.DataFrame) -> pd.DataFrame:
    return (
        test[test["label"].eq(1)][["read_id", "candidate_offset", "candidate_frame"]]
        .rename(
            columns={
                "candidate_offset": "true_offset",
                "candidate_frame": "true_frame",
            }
        )
        .drop_duplicates("read_id")
    )


def _candidate_predictions(test: pd.DataFrame, predictions: pd.DataFrame) -> pd.DataFrame:
    pred_candidates = test[["read_id", "candidate_offset", "label", "candidate_frame"]].merge(
        predictions,
        on=["read_id", "candidate_offset"],
        how="left",
    )
    pred_candidates["p_offset"] = pred_candidates["p_offset"].fillna(0.0).astype(float)
    return pred_candidates


def _top_predicted_rows(pred_candidates: pd.DataFrame) -> pd.DataFrame:
    return (
        pred_candidates.sort_values(["read_id", "p_offset", "candidate_offset"], ascending=[True, False, True])
        .groupby("read_id", as_index=False)
        .head(1)[["read_id", "candidate_offset", "candidate_frame", "p_offset"]]
        .rename(
            columns={
                "candidate_offset": "pred_offset",
                "candidate_frame": "pred_frame",
                "p_offset": "p_pred",
            }
        )
    )


def _q_true(pred_candidates: pd.DataFrame) -> pd.DataFrame:
    return (
        pred_candidates[pred_candidates["label"].eq(1)][["read_id", "p_offset"]]
        .rename(columns={"p_offset": "q_true"})
        .drop_duplicates("read_id")
    )


def evaluate_candidate_predictions(
    test: pd.DataFrame,
    predictions: pd.DataFrame,
) -> tuple[dict[str, float], pd.DataFrame]:
    _require_test_columns(test)

    predictions = normalize_prediction_columns(predictions)
    # This evaluation only scores reads with at least one prediction.
    predicted_read_ids = set(predictions["read_id"].astype(str).unique())
    test = test[test["read_id"].astype(str).isin(predicted_read_ids)].copy()
    if test.empty:
        raise ValueError("no test rows match the predicted read IDs")

    pred_candidates = _candidate_predictions(test, predictions)
    true_rows = _true_rows(test)
    pred_rows = _top_predicted_rows(pred_candidates)
    q_true = _q_true(pred_candidates)
    pred_read = pred_rows.merge(true_rows, on="read_id", how="inner").merge(
        q_true,
        on="read_id",
        how="left",
    )
    pred_read["q_true"] = pred_read["q_true"].fillna(0.0)

    metrics = offset_metrics(pred_read)
    candidate_metrics = binary_metrics(pred_candidates["label"], pred_candidates["p_offset"])
    metrics.update({f"candidate_{name}": value for name, value in candidate_metrics.items()})
    metrics["reads"] = float(pred_read["read_id"].nunique())
    metrics["candidate_rows"] = float(len(pred_candidates))
    return metrics, pred_read


def strict_offset_metrics(pred_read: pd.DataFrame) -> dict[str, float]:
    required = {"true_offset", "pred_offset", "true_frame", "pred_frame", "q_true", "has_prediction"}
    missing = sorted(required - set(pred_read.columns))
    if missing:
        raise ValueError(f"missing required columns: {missing}")

    has_prediction = pred_read["has_prediction"].astype(bool)
    # Missing predictions stay in the denominator and count as incorrect.
    exact = (has_prediction & pred_read["true_offset"].eq(pred_read["pred_offset"])).mean()
    within1 = (has_prediction & (pred_read["true_offset"] - pred_read["pred_offset"]).abs().le(1)).mean()
    frame_acc = (has_prediction & pred_read["true_frame"].eq(pred_read["pred_frame"])).mean()
    nll = -np.log(np.clip(pred_read["q_true"].astype(float), 1e-12, 1.0)).mean()
    return {
        "exact_offset_accuracy": float(exact),
        "within_1nt_accuracy": float(within1),
        "frame_accuracy": float(frame_acc),
        "read_level_nll": float(nll),
    }


def evaluate_targeted_candidate_predictions(
    test: pd.DataFrame,
    predictions: pd.DataFrame,
    target_read_ids: Iterable[object],
    missing_policy: str = "strict",
    include_candidate_metrics: bool = True,
) -> tuple[dict[str, float | str], pd.DataFrame]:
    if missing_policy != "strict":
        raise ValueError(f"unsupported missing_policy for targeted evaluation: {missing_policy!r}")
    _require_test_columns(test)

    target_ids = _unique_read_ids(target_read_ids)
    if not target_ids:
        raise ValueError("target_read_ids is empty")
    target_set = set(target_ids)

    test = test.copy()
    test["read_id"] = test["read_id"].astype(str)
    predictions = normalize_prediction_columns(predictions)
    predictions["read_id"] = predictions["read_id"].astype(str)
    predictions = predictions.groupby(["read_id", "candidate_offset"], as_index=False)["p_offset"].max()

    target_test = test[test["read_id"].isin(target_set)].copy()
    present_targets = set(target_test["read_id"].unique())
    missing_from_test = sorted(target_set - present_targets)
    if missing_from_test:
        raise ValueError(f"target read IDs are absent from the test table: {len(missing_from_test)}")
    if target_test.empty:
        raise ValueError("no test rows match the target read IDs")

    true_rows = _true_rows(target_test)
    reads_without_truth = sorted(target_set - set(true_rows["read_id"].astype(str)))
    if reads_without_truth:
        raise ValueError(f"target read IDs are missing positive true-offset rows: {len(reads_without_truth)}")

    target_candidate_keys = target_test[["read_id", "candidate_offset"]].drop_duplicates()
    matched_predictions = predictions.merge(target_candidate_keys, on=["read_id", "candidate_offset"], how="inner")
    predicted_target_ids = set(matched_predictions["read_id"].astype(str).unique())

    pred_candidates = _candidate_predictions(target_test, predictions)
    pred_rows = _top_predicted_rows(pred_candidates)
    q_true = _q_true(pred_candidates)

    target_order = pd.DataFrame({"read_id": target_ids})
    pred_read = (
        target_order.merge(true_rows, on="read_id", how="left")
        .merge(pred_rows, on="read_id", how="left")
        .merge(q_true, on="read_id", how="left")
    )
    pred_read["has_prediction"] = pred_read["read_id"].isin(predicted_target_ids)
    missing_mask = ~pred_read["has_prediction"]
    pred_read.loc[missing_mask, ["pred_offset", "pred_frame", "p_pred"]] = np.nan
    pred_read["q_true"] = pred_read["q_true"].fillna(0.0)

    target_reads = len(target_ids)
    predicted_target_reads = len(predicted_target_ids)
    metrics: dict[str, float | str] = strict_offset_metrics(pred_read)
    metrics.update(
        {
            "target_reads": float(target_reads),
            "predicted_target_reads": float(predicted_target_reads),
            "missing_target_reads": float(target_reads - predicted_target_reads),
            "coverage": float(predicted_target_reads / target_reads),
            "reads": float(target_reads),
            "candidate_rows": float(len(pred_candidates)),
        }
    )
    if include_candidate_metrics:
        candidate_metrics = binary_metrics(pred_candidates["label"], pred_candidates["p_offset"])
        metrics.update({f"candidate_{name}": value for name, value in candidate_metrics.items()})
        metrics["candidate_metric_scope"] = "common_test_subset_candidate_strict"
    return metrics, pred_read


def _filter_to_target_reads(
    test: pd.DataFrame,
    predictions: pd.DataFrame,
    target_read_ids: Iterable[object],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    target = set(_unique_read_ids(target_read_ids))
    test_out = test[test["read_id"].astype(str).isin(target)].copy()
    pred_out = predictions[predictions["read_id"].astype(str).isin(target)].copy()
    return test_out, pred_out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Evaluate P-site offset predictions for each read.")
    parser.add_argument("--test", type=pathlib.Path, required=True)
    parser.add_argument("--predictions", type=pathlib.Path, required=True)
    parser.add_argument("--metrics-output", type=pathlib.Path)
    parser.add_argument("--read-output", type=pathlib.Path)
    parser.add_argument("--target-read-ids", type=pathlib.Path)
    parser.add_argument("--missing-policy", choices=("conditional", "strict"), default="conditional")
    parser.add_argument("--prediction-filter-mode", choices=("auto", "always", "never"), default="auto")
    parser.add_argument("--prediction-stream-threshold", type=int, default=5_000_000)
    parser.add_argument("--prediction-batch-size", type=int, default=100_000)
    parser.add_argument("--model", default="")
    parser.add_argument("--split", default="")
    args = parser.parse_args(argv)

    test = pd.read_parquet(args.test)
    target_read_ids = load_read_subset(args.target_read_ids) if args.target_read_ids else None
    if args.prediction_filter_mode == "never":
        prediction_filter_ids = None
        stream_threshold = args.prediction_stream_threshold
    else:
        prediction_filter_ids = target_read_ids if target_read_ids is not None else test["read_id"]
        stream_threshold = 0 if args.prediction_filter_mode == "always" else args.prediction_stream_threshold
    predictions = read_prediction_table(
        args.predictions,
        filter_read_ids=prediction_filter_ids,
        stream_threshold=stream_threshold,
        batch_size=args.prediction_batch_size,
    )
    if args.missing_policy == "strict":
        if target_read_ids is None:
            raise ValueError("--missing-policy strict requires --target-read-ids")
        metrics, pred_read = evaluate_targeted_candidate_predictions(
            test,
            predictions,
            target_read_ids,
            missing_policy="strict",
        )
        metrics["metric_scope"] = "common_test_subset"
    else:
        if target_read_ids is not None:
            test, predictions = _filter_to_target_reads(test, predictions, target_read_ids)
        metrics, pred_read = evaluate_candidate_predictions(test, predictions)
        metrics["metric_scope"] = "common_test_subset_predicted_only" if target_read_ids is not None else "all_predicted_reads"
    if args.model:
        metrics["model"] = args.model
    if args.split:
        metrics["split"] = args.split

    metrics_frame = pd.DataFrame([metrics])
    if args.metrics_output:
        args.metrics_output.parent.mkdir(parents=True, exist_ok=True)
        metrics_frame.to_csv(args.metrics_output, sep=table_sep(args.metrics_output), index=False)
    if args.read_output:
        args.read_output.parent.mkdir(parents=True, exist_ok=True)
        pred_read.to_parquet(args.read_output, index=False)
    print(metrics_frame.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
