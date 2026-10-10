from __future__ import annotations

import pandas as pd


def apply_fixed_offset(
    reads: pd.DataFrame,
    global_offset: int = 12,
    five_prime_col: str = "five_prime_pos",
) -> pd.DataFrame:
    if five_prime_col not in reads.columns:
        raise ValueError(f"missing required column: {five_prime_col}")
    out = reads.copy()
    out["p_site_pos_fixed"] = out[five_prime_col] + global_offset
    out["a_site_pos_fixed"] = out["p_site_pos_fixed"] + 3
    return out


def estimate_offset_by_length(
    train_landmark_reads: pd.DataFrame,
    read_length_col: str = "read_length",
    silver_offset_col: str = "silver_offset",
) -> dict[int, int]:
    required = {read_length_col, silver_offset_col}
    missing = sorted(required - set(train_landmark_reads.columns))
    if missing:
        raise ValueError(f"missing required columns: {missing}")

    return (
        train_landmark_reads.groupby(read_length_col)[silver_offset_col]
        .agg(lambda s: int(s.value_counts().idxmax()))
        .to_dict()
    )


def apply_read_length_offsets(
    reads: pd.DataFrame,
    offset_by_length: dict[int, int],
    global_offset: int = 12,
    read_length_col: str = "read_length",
    five_prime_col: str = "five_prime_pos",
) -> pd.DataFrame:
    required = {read_length_col, five_prime_col}
    missing = sorted(required - set(reads.columns))
    if missing:
        raise ValueError(f"missing required columns: {missing}")

    out = reads.copy()
    out["offset_len"] = out[read_length_col].map(offset_by_length).fillna(global_offset).astype(int)
    out["p_site_pos_len"] = out[five_prime_col] + out["offset_len"]
    out["a_site_pos_len"] = out["p_site_pos_len"] + 3
    return out

