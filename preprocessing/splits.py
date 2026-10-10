from __future__ import annotations

import argparse
import hashlib
import pathlib
import re

import pandas as pd


DEFAULT_HOLDOUT_COLUMNS = ("study_id", "nuclease", "organism")


def slugify(value: object) -> str:
    text = str(value).strip().lower()
    text = re.sub(r"[^a-z0-9._-]+", "_", text)
    return text.strip("_") or "missing"


def stable_score(value: str, seed: int) -> int:
    digest = hashlib.blake2b(f"{seed}:{value}".encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, byteorder="big", signed=False)


def validate_candidate_table(candidates: pd.DataFrame) -> pd.DataFrame:
    required = {"read_id", "label", "candidate_offset"}
    missing = sorted(required - set(candidates.columns))
    if missing:
        raise ValueError(f"candidate table is missing required columns: {missing}")

    labeled = candidates[candidates["label"].notna()].copy()
    labeled["label"] = labeled["label"].astype("int8")
    positive_counts = labeled.groupby("read_id")["label"].sum()
    valid_read_ids = positive_counts[positive_counts.eq(1)].index
    if len(valid_read_ids) == 0:
        raise ValueError("no read groups with exactly one positive candidate were found")
    return labeled[labeled["read_id"].isin(valid_read_ids)].copy()


def split_counts(frame: pd.DataFrame, split_name: str, split_type: str) -> dict[str, object]:
    return {
        "split_name": split_name,
        "split_type": split_type,
        "rows": len(frame),
        "reads": frame["read_id"].nunique(),
        "positive_rows": int(frame["label"].sum()) if "label" in frame.columns else "",
    }


def write_split(
    candidates: pd.DataFrame,
    reads: pd.DataFrame | None,
    train_read_ids: set[str],
    test_read_ids: set[str],
    split_dir: pathlib.Path,
    split_name: str,
    split_type: str,
    holdout_column: str = "",
    holdout_value: str = "",
) -> dict[str, object]:
    split_dir.mkdir(parents=True, exist_ok=True)
    train = candidates[candidates["read_id"].isin(train_read_ids)].copy()
    test = candidates[candidates["read_id"].isin(test_read_ids)].copy()
    if train.empty or test.empty:
        raise ValueError(f"empty train/test split for {split_name}")

    train_path = split_dir / "train.parquet"
    test_path = split_dir / "test.parquet"
    train.to_parquet(train_path, index=False)
    test.to_parquet(test_path, index=False)

    train_reads_path = ""
    test_reads_path = ""
    if reads is not None:
        train_reads = reads[reads["read_id"].isin(train_read_ids)].copy()
        test_reads = reads[reads["read_id"].isin(test_read_ids)].copy()
        train_reads_path = split_dir / "train_reads.parquet"
        test_reads_path = split_dir / "test_reads.parquet"
        train_reads.to_parquet(train_reads_path, index=False)
        test_reads.to_parquet(test_reads_path, index=False)

    return {
        "split_name": split_name,
        "split_type": split_type,
        "holdout_column": holdout_column,
        "holdout_value": holdout_value,
        "train_rows": len(train),
        "test_rows": len(test),
        "train_reads": train["read_id"].nunique(),
        "test_reads": test["read_id"].nunique(),
        "train_positive_rows": int(train["label"].sum()),
        "test_positive_rows": int(test["label"].sum()),
        "train_path": str(train_path),
        "test_path": str(test_path),
        "train_reads_path": str(train_reads_path),
        "test_reads_path": str(test_reads_path),
    }


def deterministic_group_split(
    candidates: pd.DataFrame,
    group_col: str,
    test_frac: float,
    seed: int,
) -> tuple[set[str], set[str]]:
    train_read_ids: set[str] = set()
    test_read_ids: set[str] = set()
    grouped_reads = candidates[["read_id", group_col]].drop_duplicates()

    for _group_value, group in grouped_reads.groupby(group_col, dropna=False):
        read_ids = sorted(group["read_id"].astype(str).unique())
        ranked = sorted(read_ids, key=lambda read_id: stable_score(read_id, seed))
        test_n = max(1, int(round(len(ranked) * test_frac))) if len(ranked) > 1 else 1
        test_set = set(ranked[:test_n])
        train_set = set(ranked[test_n:])
        if not train_set and len(test_set) > 1:
            moved = sorted(test_set)[0]
            test_set.remove(moved)
            train_set.add(moved)
        train_read_ids.update(train_set)
        test_read_ids.update(test_set)

    return train_read_ids, test_read_ids


def create_splits(
    candidates: pd.DataFrame,
    reads: pd.DataFrame | None,
    outdir: pathlib.Path,
    holdout_columns: tuple[str, ...],
    test_frac: float,
    seed: int,
) -> pd.DataFrame:
    manifests: list[dict[str, object]] = []

    train_ids, test_ids = deterministic_group_split(
        candidates,
        group_col="sample_id" if "sample_id" in candidates.columns else "study_id",
        test_frac=test_frac,
        seed=seed,
    )
    manifests.append(
        write_split(
            candidates,
            reads,
            train_ids,
            test_ids,
            outdir / "within_sample",
            split_name="within_sample",
            split_type="deterministic_grouped_within_sample",
        )
    )

    all_read_ids = set(candidates["read_id"].astype(str).unique())
    for column in holdout_columns:
        if column not in candidates.columns:
            continue
        for value in sorted(candidates[column].dropna().astype(str).unique()):
            test_read_ids = set(
                candidates.loc[candidates[column].astype(str).eq(value), "read_id"]
                .astype(str)
                .unique()
            )
            train_read_ids = all_read_ids - test_read_ids
            if not train_read_ids or not test_read_ids:
                continue
            split_type = f"leave_one_{column}"
            split_name = f"{split_type}={slugify(value)}"
            manifests.append(
                write_split(
                    candidates,
                    reads,
                    train_read_ids,
                    test_read_ids,
                    outdir / split_type / slugify(value),
                    split_name=split_name,
                    split_type=split_type,
                    holdout_column=column,
                    holdout_value=value,
                )
            )

    return pd.DataFrame(manifests)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Create p-site offset prediction train/test Parquet splits.")
    parser.add_argument(
        "--candidates",
        type=pathlib.Path,
        default=pathlib.Path(
            "data/datasets/p_site_offset_prediction/landmark_riboseq/tables/candidate_offsets.parquet"
        ),
    )
    parser.add_argument(
        "--reads",
        type=pathlib.Path,
        default=pathlib.Path("data/datasets/p_site_offset_prediction/landmark_riboseq/tables/reads.parquet"),
    )
    parser.add_argument(
        "--outdir",
        type=pathlib.Path,
        default=pathlib.Path("data/datasets/p_site_offset_prediction/landmark_riboseq/splits"),
    )
    parser.add_argument("--test-frac", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=20260517)
    parser.add_argument("--holdout-columns", nargs="*", default=list(DEFAULT_HOLDOUT_COLUMNS))
    args = parser.parse_args(argv)

    if not 0.0 < args.test_frac < 1.0:
        raise ValueError("--test-frac must be between 0 and 1")

    candidates = validate_candidate_table(pd.read_parquet(args.candidates))
    reads = None
    if args.reads.exists():
        reads = pd.read_parquet(args.reads)
        reads = reads[reads["read_id"].isin(candidates["read_id"].unique())].copy()

    args.outdir.mkdir(parents=True, exist_ok=True)
    manifest = create_splits(
        candidates=candidates,
        reads=reads,
        outdir=args.outdir,
        holdout_columns=tuple(args.holdout_columns),
        test_frac=args.test_frac,
        seed=args.seed,
    )
    manifest_path = args.outdir / "split_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)
    print(f"wrote {manifest_path} with {len(manifest)} splits")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
