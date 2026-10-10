from __future__ import annotations

import hashlib

import numpy as np
import pandas as pd

OFFSET_PRIOR_HIERARCHY: tuple[tuple[str, ...], ...] = (
    ("sample_id", "read_length"),
    ("study_id", "read_length"),
    ("organism", "nuclease", "library_method", "read_length"),
    ("organism", "nuclease", "read_length"),
    ("organism", "read_length"),
    ("read_length",),
)


def stable_score(value: str, seed: int) -> int:
    digest = hashlib.blake2b(f"{seed}:{value}".encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, byteorder="big", signed=False)


def sample_read_groups(frame: pd.DataFrame, max_reads: int | None, seed: int) -> pd.DataFrame:
    if max_reads is None:
        return frame
    read_ids = sorted(frame["read_id"].astype(str).unique(), key=lambda value: stable_score(value, seed))
    keep = set(read_ids[:max_reads])
    return frame[frame["read_id"].astype(str).isin(keep)].copy()


def normalize_candidate_probabilities(test: pd.DataFrame, raw_col: str = "p_raw") -> pd.Series:
    return test.groupby("read_id")[raw_col].transform(lambda s: s / max(float(s.sum()), 1e-12))


def offset_values(train: pd.DataFrame, test: pd.DataFrame) -> list[int]:
    values = pd.concat(
        [
            train["candidate_offset"].astype(int),
            test["candidate_offset"].astype(int),
        ],
        ignore_index=True,
    )
    return sorted(values.dropna().unique().tolist())


def offset_prior_table(
    positives: pd.DataFrame,
    key_cols: list[str],
    offsets: list[int],
    smoothing: float,
) -> pd.DataFrame:
    if key_cols:
        groups = positives[key_cols].drop_duplicates()
        grid = groups.merge(pd.DataFrame({"candidate_offset": offsets}), how="cross")
        count_keys = key_cols + ["candidate_offset"]
    else:
        grid = pd.DataFrame({"candidate_offset": offsets})
        count_keys = ["candidate_offset"]

    counts = positives.groupby(count_keys, dropna=False).size().rename("count").reset_index()
    table = grid.merge(counts, on=count_keys, how="left")
    table["count"] = table["count"].fillna(0.0).astype(float) + smoothing
    if key_cols:
        table["offset_prior"] = table["count"] / table.groupby(key_cols, dropna=False)["count"].transform("sum")
    else:
        table["offset_prior"] = table["count"] / table["count"].sum()
    return table[count_keys + ["offset_prior"]]


def estimate_offset_prior(
    train: pd.DataFrame,
    test: pd.DataFrame,
    smoothing: float,
) -> pd.Series:
    if smoothing <= 0:
        raise ValueError("--offset-prior-smoothing must be positive")
    required = {"candidate_offset", "label"}
    missing = sorted(required - set(train.columns))
    if missing:
        raise ValueError(f"training data is missing columns needed for offset priors: {missing}")
    if "candidate_offset" not in test.columns:
        raise ValueError("test data is missing candidate_offset needed for offset priors")

    positives = train[train["label"].eq(1)].copy()
    if positives.empty:
        return pd.Series(1.0, index=test.index, dtype=float)

    offsets = offset_values(train, test)
    if not offsets:
        return pd.Series(1.0, index=test.index, dtype=float)

    prior = pd.Series(np.nan, index=test.index, dtype=float)
    work = test.reset_index(drop=True).copy()
    work["_row_index"] = np.arange(len(work))
    for key_cols_tuple in OFFSET_PRIOR_HIERARCHY:
        key_cols = [
            column
            for column in key_cols_tuple
            if column in positives.columns and column in work.columns
        ]
        if len(key_cols) != len(key_cols_tuple):
            continue
        table = offset_prior_table(positives, key_cols, offsets, smoothing)
        merged = work[["_row_index", *key_cols, "candidate_offset"]].merge(
            table,
            on=[*key_cols, "candidate_offset"],
            how="left",
            sort=False,
        )
        merged = merged.sort_values("_row_index")
        values = pd.Series(merged["offset_prior"].to_numpy(), index=test.index, dtype=float)
        prior = prior.fillna(values)

    global_table = offset_prior_table(positives, [], offsets, smoothing)
    global_prior = test["candidate_offset"].astype(int).map(
        dict(zip(global_table["candidate_offset"].astype(int), global_table["offset_prior"]))
    )
    return prior.fillna(global_prior).fillna(1.0 / len(offsets)).astype(float)


def score_candidate_probabilities(
    test: pd.DataFrame,
    raw_col: str = "p_raw",
    prior_col: str = "offset_prior",
    raw_power: float = 1.0,
    prior_power: float = 0.0,
) -> pd.Series:
    raw = test[raw_col].astype(float).clip(1e-12, 1.0)
    score = raw.pow(raw_power)
    if prior_power:
        score = score * test[prior_col].astype(float).clip(1e-12, 1.0).pow(prior_power)
    return test.assign(p_score=score).groupby("read_id")["p_score"].transform(
        lambda s: s / max(float(s.sum()), 1e-12)
    )
