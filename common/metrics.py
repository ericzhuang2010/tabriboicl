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


def expected_calibration_error(y_true, p, n_bins: int = 15) -> float:
    """Compare the mean score and positive rate in each probability bin."""
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
    """Score candidate rows. The evaluator adds the 'candidate_' prefix."""
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
        # AUROC needs both positive and negative labels.
        out["auroc"] = np.nan
    return {key: float(value) for key, value in out.items()}


def offset_metrics(pred_read_df: pd.DataFrame) -> dict[str, float]:
    """Score one chosen offset per read; q_true is the true offset's probability."""
    required = {"true_offset", "pred_offset", "true_frame", "pred_frame", "q_true"}
    missing = sorted(required - set(pred_read_df.columns))
    if missing:
        raise ValueError(f"missing required columns: {missing}")

    exact = (pred_read_df["true_offset"] == pred_read_df["pred_offset"]).mean()
    within1 = (pred_read_df["true_offset"] - pred_read_df["pred_offset"]).abs().le(1).mean()
    # An offset off by a multiple of three can still get the frame right.
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
