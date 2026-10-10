#!/usr/bin/env python3
from __future__ import annotations

import os
import random


class _NoOpFacetGrid:
    def __init__(self, *args, **kwargs) -> None:
        pass

    def map(self, *args, **kwargs):
        return self

    def set(self, *args, **kwargs):
        return self

    def savefig(self, *args, **kwargs):
        return self


def _truthy_env(name: str, default: str = "1") -> bool:
    value = os.environ.get(name, default).strip().lower()
    return value not in {"0", "false", "no", "off"}


def _disable_train_qc_plots() -> None:
    import seaborn as sns

    sns.FacetGrid = _NoOpFacetGrid
    sns.histplot = lambda *args, **kwargs: None


def _align_sklearn_prediction_features() -> None:
    """Align independently encoded PSite chunks to fitted GBT columns."""
    from sklearn.ensemble import GradientBoostingClassifier

    if getattr(GradientBoostingClassifier.predict, "_tabribo_feature_aligned", False):
        return
    original_predict = GradientBoostingClassifier.predict

    def predict_with_aligned_features(self, X, *args, **kwargs):
        feature_names = getattr(self, "feature_names_in_", None)
        if feature_names is not None and hasattr(X, "reindex"):
            X = X.reindex(columns=list(feature_names), fill_value=False)
        return original_predict(self, X, *args, **kwargs)

    predict_with_aligned_features._tabribo_feature_aligned = True
    GradientBoostingClassifier.predict = predict_with_aligned_features


def main() -> int:
    random.seed(int(os.environ.get("TABRIBO_PSITE_RANDOM_SEED", "20260517")))
    if _truthy_env("TABRIBO_PSITE_DISABLE_QC_PLOTS"):
        _disable_train_qc_plots()
    _align_sklearn_prediction_features()

    from psite.main import psite

    return psite()


if __name__ == "__main__":
    raise SystemExit(main())
