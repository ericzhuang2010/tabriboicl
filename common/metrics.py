from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    log_loss,
    roc_auc_score,
)


# Metric glossary:
#
# read-level cares about the final chosen answer per read.
# candidate-level cares about the quality/calibration/ranking of all candidate scores.
#
# offset_metrics(pred_read_df) returns read-level offset metrics:
# - exact_offset_accuracy: fraction of reads where pred_offset exactly equals
#   true_offset.
# - within_1nt_accuracy: fraction of reads where pred_offset is within one
#   nucleotide of true_offset.
# - frame_accuracy: fraction of reads where pred_frame equals true_frame.
#   This can be correct even when the exact offset is wrong by a multiple of
#   three nucleotides.
# - read_level_nll: negative log likelihood of the true read-level offset,
#   computed from q_true (the probability assigned to the true candidate).
#
# binary_metrics(y_true, p) returns candidate-level binary classification
# metrics before naming prefixes are applied:
# - accuracy: fraction of candidate rows whose thresholded prediction
#   (p >= 0.5) matches the binary label.
# - balanced_accuracy: mean of sensitivity and specificity after thresholding
#   at p >= 0.5; useful when positive candidate rows are rare.
# - auprc: area under the precision-recall curve using raw probabilities.
# - brier: mean squared error between the binary labels and probabilities.
# - nll: binary negative log likelihood / cross entropy.
# - ece: expected calibration error; lower means predicted probabilities better
#   match observed positive rates across probability bins.
# - auroc: area under the ROC curve using raw probabilities. This is NaN when
#   the labels contain only one class.
#
# evaluation.evaluate_candidate_predictions() prefixes binary_metrics outputs
# with "candidate_": for example, auprc becomes candidate_auprc. It also adds
# reads, the number of evaluated read IDs, and candidate_rows, the number of
# candidate rows included in the candidate-level evaluation.


def expected_calibration_error(y_true, p, n_bins: int = 15) -> float:
    y_true = np.asarray(y_true)
    p = np.asarray(p)
    bins = np.linspace(0.0, 1.0, n_bins + 1)
    ece = 0.0
    for lo, hi in zip(bins[:-1], bins[1:]):
        mask = (p >= lo) & (p < hi)
        if mask.sum() == 0:
            continue
        conf = p[mask].mean()
        acc = y_true[mask].mean()
        ece += mask.mean() * abs(acc - conf)
    return float(ece)


def binary_metrics(y_true, p) -> dict[str, float]:
    y_true = np.asarray(y_true)
    p = np.asarray(p)
    pred = (p >= 0.5).astype(int)
    out = {
        "accuracy": accuracy_score(y_true, pred),
        "balanced_accuracy": balanced_accuracy_score(y_true, pred),
        "auprc": average_precision_score(y_true, p),
        "brier": brier_score_loss(y_true, p),
        "nll": log_loss(y_true, np.vstack([1 - p, p]).T, labels=[0, 1]),
        "ece": expected_calibration_error(y_true, p),
    }
    try:
        out["auroc"] = roc_auc_score(y_true, p)
    except ValueError:
        out["auroc"] = np.nan
    return {key: float(value) for key, value in out.items()}


def offset_metrics(pred_read_df: pd.DataFrame) -> dict[str, float]:
    required = {"true_offset", "pred_offset", "true_frame", "pred_frame", "q_true"}
    missing = sorted(required - set(pred_read_df.columns))
    if missing:
        raise ValueError(f"missing required columns: {missing}")

    exact = (pred_read_df["true_offset"] == pred_read_df["pred_offset"]).mean()
    within1 = (pred_read_df["true_offset"] - pred_read_df["pred_offset"]).abs().le(1).mean()
    frame_acc = (pred_read_df["true_frame"] == pred_read_df["pred_frame"]).mean()
    nll = -np.log(np.clip(pred_read_df["q_true"], 1e-12, 1.0)).mean()
    return {
        "exact_offset_accuracy": float(exact),
        "within_1nt_accuracy": float(within1),
        "frame_accuracy": float(frame_acc),
        "read_level_nll": float(nll),
    }


def entropy_by_read(candidate_df: pd.DataFrame) -> pd.DataFrame:
    required = {"read_id", "p_offset"}
    missing = sorted(required - set(candidate_df.columns))
    if missing:
        raise ValueError(f"missing required columns: {missing}")

    out = candidate_df[["read_id", "p_offset"]].copy()
    q = out["p_offset"].clip(1e-12, 1.0)
    out["entropy_term"] = -q * np.log(q)
    return (
        out.groupby("read_id", as_index=False)["entropy_term"]
        .sum()
        .rename(columns={"entropy_term": "offset_entropy"})
    )
