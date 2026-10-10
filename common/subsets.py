from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable, Mapping
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

READ_ID_COLUMN = "read_id"
READ_ID_HASH_ALGORITHM = "blake2b_64"


def stable_score(value: str, seed: int) -> int:
    digest = hashlib.blake2b(f"{seed}:{value}".encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, byteorder="big", signed=False)


def normalize_max_reads(max_reads: int | str | None) -> int | None:
    if max_reads is None:
        return None
    if isinstance(max_reads, str):
        text = max_reads.strip().lower()
        if text in {"", "all", "full", "none", "0"}:
            return None
        if not text.isdigit():
            raise ValueError(f"max_reads must be an integer or 'all', got {max_reads!r}")
        max_reads = int(text)
    if max_reads < 0:
        raise ValueError(f"max_reads must be non-negative, got {max_reads}")
    if max_reads == 0:
        return None
    return int(max_reads)


def _read_ids_from_iterable(read_ids: Iterable[object]) -> set[str]:
    return {str(read_id) for read_id in read_ids}


def select_read_ids(
    frame: pd.DataFrame,
    max_reads: int | str | None,
    seed: int,
    read_id_column: str = READ_ID_COLUMN,
) -> list[str]:
    if read_id_column not in frame.columns:
        raise ValueError(f"frame is missing required column: {read_id_column}")

    limit = normalize_max_reads(max_reads)
    read_ids = sorted(
        frame[read_id_column].astype(str).unique(),
        key=lambda value: (stable_score(value, seed), value),
    )
    if limit is None:
        return list(read_ids)
    return list(read_ids[:limit])


def filter_read_groups(
    frame: pd.DataFrame,
    read_ids: Iterable[object],
    read_id_column: str = READ_ID_COLUMN,
) -> pd.DataFrame:
    if read_id_column not in frame.columns:
        raise ValueError(f"frame is missing required column: {read_id_column}")
    keep = _read_ids_from_iterable(read_ids)
    return frame[frame[read_id_column].astype(str).isin(keep)].copy()


def read_subset_metadata(
    *,
    test_path: str | Path,
    test_reads_requested: int | str | None,
    test_reads_selected: int,
    selection_seed: int,
    selection_method: str = "stable_hash",
    created_by: str = "common.subsets",
    extra: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    path = Path(test_path)
    metadata: dict[str, Any] = {
        "test_path": str(path),
        "test_path_resolved": str(path.resolve()),
        "test_reads_requested": test_reads_requested,
        "test_reads_selected": int(test_reads_selected),
        "selection_seed": int(selection_seed),
        "selection_method": selection_method,
        "read_id_hash_algorithm": READ_ID_HASH_ALGORITHM,
        "created_by": created_by,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    if extra:
        metadata.update(dict(extra))
    return metadata


def write_read_subset(
    path: str | Path,
    read_ids: Iterable[object],
    metadata: Mapping[str, Any] | None = None,
) -> None:
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ids = [str(read_id) for read_id in read_ids]
    pd.DataFrame({READ_ID_COLUMN: ids}).to_csv(out_path, sep="\t", index=False)
    if metadata is not None:
        metadata_path = out_path.with_name(out_path.stem.replace("_reads", "") + "_metadata.json")
        metadata_path.write_text(json.dumps(dict(metadata), indent=2, sort_keys=True) + "\n")


def load_read_subset(path: str | Path, read_id_column: str = READ_ID_COLUMN) -> list[str]:
    frame = pd.read_csv(path, sep="\t")
    if read_id_column not in frame.columns:
        raise ValueError(f"subset file is missing required column: {read_id_column}")
    return frame[read_id_column].astype(str).tolist()


def assert_prediction_read_set_matches_target(
    predictions: pd.DataFrame,
    target_read_ids: Iterable[object],
    read_id_column: str = READ_ID_COLUMN,
) -> None:
    if read_id_column not in predictions.columns:
        raise ValueError(f"predictions are missing required column: {read_id_column}")
    predicted = set(predictions[read_id_column].astype(str).unique())
    target = _read_ids_from_iterable(target_read_ids)
    missing = sorted(target - predicted)
    extra = sorted(predicted - target)
    if missing or extra:
        parts = []
        if missing:
            parts.append(f"missing target reads: {len(missing)}")
        if extra:
            parts.append(f"extra non-target reads: {len(extra)}")
        raise ValueError("prediction read IDs do not match target read IDs; " + ", ".join(parts))
