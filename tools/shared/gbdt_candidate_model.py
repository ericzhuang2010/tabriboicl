from __future__ import annotations

import argparse
import pathlib

import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from common.feature_sets import select_feature_columns
from common.subsets import filter_read_groups, load_read_subset
from tools.shared.learned_candidate_models import normalize_candidate_probabilities, sample_read_groups


def make_one_hot_encoder() -> OneHotEncoder:
    try:
        return OneHotEncoder(handle_unknown="ignore", sparse_output=True)
    except TypeError:
        return OneHotEncoder(handle_unknown="ignore", sparse=True)


def make_model(model_name: str, seed: int):
    if model_name == "xgboost":
        from xgboost import XGBClassifier

        return XGBClassifier(
            n_estimators=250,
            max_depth=4,
            learning_rate=0.05,
            subsample=0.8,
            colsample_bytree=0.8,
            objective="binary:logistic",
            eval_metric="logloss",
            tree_method="hist",
            random_state=seed,
            n_jobs=4,
        )
    if model_name == "lightgbm":
        from lightgbm import LGBMClassifier

        return LGBMClassifier(
            n_estimators=250,
            learning_rate=0.05,
            num_leaves=31,
            subsample=0.8,
            colsample_bytree=0.8,
            random_state=seed,
            n_jobs=4,
            verbose=-1,
        )
    if model_name == "logreg":
        return LogisticRegression(max_iter=1000, n_jobs=4, class_weight="balanced")
    raise ValueError(f"unsupported model: {model_name}")


def is_categorical_dtype(dtype: object) -> bool:
    dtype_name = str(dtype)
    return dtype_name in {"object", "str", "string", "category"} or dtype_name.startswith("string")


def prepare_pipeline(train: pd.DataFrame, feature_cols: list[str], model_name: str, seed: int) -> Pipeline:
    categorical_cols = [
        column
        for column in feature_cols
        if is_categorical_dtype(train[column].dtype)
    ]
    numeric_cols = [column for column in feature_cols if column not in categorical_cols]
    preprocessor = ColumnTransformer(
        transformers=[
            ("cat", make_one_hot_encoder(), categorical_cols),
            ("num", StandardScaler(with_mean=False), numeric_cols),
        ],
        sparse_threshold=1.0,
    )
    return Pipeline(
        steps=[
            ("preprocess", preprocessor),
            ("model", make_model(model_name, seed)),
        ]
    )


def clean_features(frame: pd.DataFrame, feature_cols: list[str]) -> pd.DataFrame:
    out = frame[feature_cols].copy()
    for column in out.columns:
        if is_categorical_dtype(out[column].dtype):
            out[column] = out[column].fillna("missing").astype(str)
        else:
            out[column] = out[column].fillna(-1)
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Train p-site offset prediction XGBoost/LightGBM/logistic comparators.")
    parser.add_argument("--train", type=pathlib.Path, required=True)
    parser.add_argument("--test", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--top-output", type=pathlib.Path)
    parser.add_argument("--model", choices=("xgboost", "lightgbm", "logreg"), required=True)
    parser.add_argument("--feature-set", choices=("full", "no_boundary", "minimal"), default="no_boundary")
    parser.add_argument("--max-train-reads", type=int)
    parser.add_argument("--max-test-reads", type=int)
    parser.add_argument("--target-read-ids", type=pathlib.Path)
    parser.add_argument("--seed", type=int, default=20260517)
    args = parser.parse_args(argv)

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
        raise ValueError("no feature columns are available for this feature set")

    pipeline = prepare_pipeline(train, feature_cols, args.model, args.seed)
    pipeline.fit(clean_features(train, feature_cols), train["label"].astype(int))

    out = test.copy()
    out["p_raw"] = pipeline.predict_proba(clean_features(out, feature_cols))[:, 1]
    out["p_offset"] = normalize_candidate_probabilities(out)
    pred = (
        out.sort_values(["read_id", "p_offset"], ascending=[True, False])
        .groupby("read_id", as_index=False)
        .head(1)
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out[["read_id", "candidate_offset", "p_raw", "p_offset"]].to_parquet(args.output, index=False)
    if args.top_output:
        args.top_output.parent.mkdir(parents=True, exist_ok=True)
        pred[["read_id", "candidate_offset", "p_offset"]].to_parquet(args.top_output, index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
