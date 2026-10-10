#!/usr/bin/env python3
from __future__ import annotations

import argparse

import csv
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common.p_site_offset_prediction.evaluation import (
    evaluate_candidate_predictions,
    evaluate_targeted_candidate_predictions,
    read_prediction_table,
)
from common.p_site_offset_prediction.subsets import (
    assert_prediction_read_set_matches_target,
    load_read_subset,
    read_subset_metadata,
    select_read_ids,
    write_read_subset,
)

LEARNED_MODEL_FAMILIES = ("tabicl", "xgboost", "lightgbm")
DIRECT_EXTERNAL_TOOLS = ("psite", "ribowaltz", "asite_ip")
PROFILE_EXTERNAL_TOOLS = ("ribodeblur", "riboba", "scikit_ribo")
OFFSET_BASELINES = ("fixed_offset", "read_length_offset")
DEFAULT_CONFIG = Path("configs/p_site_offset_prediction_landmark_riboseq.tsv")


@dataclass(frozen=True)
class ExperimentConfig:
    experiment: str
    enabled: bool
    algorithms: tuple[str, ...]
    train_specs: tuple[str, ...]
    compare_train_spec: str


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)
    print(f"wrote {path}")


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def split_list(value: object, default: tuple[str, ...] = ()) -> tuple[str, ...]:
    text = str(value or "").replace(",", " ").replace(";", " ")
    items = tuple(item.strip() for item in text.split() if item.strip())
    return items or default


def expand_algorithm_aliases(algorithms: tuple[str, ...]) -> tuple[str, ...]:
    out: list[str] = []
    for algorithm in algorithms:
        if algorithm == "all":
            out.extend([*OFFSET_BASELINES, *LEARNED_MODEL_FAMILIES, *DIRECT_EXTERNAL_TOOLS, *PROFILE_EXTERNAL_TOOLS])
        elif algorithm == "offset_baselines":
            out.extend(OFFSET_BASELINES)
        elif algorithm == "learned":
            out.extend(LEARNED_MODEL_FAMILIES)
        elif algorithm == "external":
            out.extend([*DIRECT_EXTERNAL_TOOLS, *PROFILE_EXTERNAL_TOOLS])
        else:
            out.append(algorithm)
    return tuple(dict.fromkeys(out))


def read_eval_config(path: Path) -> list[ExperimentConfig]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        configs = []
        for row in reader:
            experiment = str(row.get("experiment", "")).strip()
            if not experiment or experiment.startswith("#"):
                continue
            algorithms = expand_algorithm_aliases(split_list(row.get("algorithms")))
            train_specs = tuple(normalize_train_spec(spec) for spec in split_list(row.get("train_specs"), ("all",)))
            compare_train_spec = normalize_train_spec(
                str(row.get("compare_train_spec") or train_specs[-1] or "all")
            )
            configs.append(
                ExperimentConfig(
                    experiment=experiment,
                    enabled=truthy(row.get("enabled", "1")),
                    algorithms=algorithms,
                    train_specs=train_specs,
                    compare_train_spec=compare_train_spec,
                )
            )
    return configs


def config_for_experiment(
    configs: list[ExperimentConfig],
    experiment: str,
) -> ExperimentConfig | None:
    aliases = experiment_aliases(experiment)
    for config in reversed(configs):
        if config.experiment in aliases:
            return config if config.enabled else None

    for config in reversed(configs):
        if not config.experiment.endswith("/*"):
            continue
        prefix = config.experiment[:-1]
        if any(alias.startswith(prefix) for alias in aliases):
            return config if config.enabled else None
    return None


def experiment_aliases(experiment: str) -> tuple[str, ...]:
    aliases = [experiment]
    if experiment.startswith("leave_one_study_id/"):
        aliases.append("leave_one_study/" + experiment.split("/", 1)[1])
    return tuple(dict.fromkeys(aliases))


def configs_for_experiment(
    configs: list[ExperimentConfig],
    experiment: str,
) -> list[ExperimentConfig]:
    aliases = experiment_aliases(experiment)
    exact = [
        config
        for config in configs
        if config.enabled and config.experiment in aliases
    ]
    if exact:
        return exact

    matches = []
    for config in configs:
        if not config.enabled or not config.experiment.endswith("/*"):
            continue
        prefix = config.experiment[:-1]
        if any(alias.startswith(prefix) for alias in aliases):
            matches.append(config)
    return matches


def normalize_train_spec(spec: str) -> str:
    spec = str(spec).strip().lower()
    if spec in {"all", "full", "none", ""}:
        return "all"
    if spec.startswith("k") and spec[1:].isdigit():
        return spec[1:]
    return spec


def model_suffix(train_spec: str) -> str:
    spec = normalize_train_spec(train_spec)
    return "all" if spec == "all" else f"k{spec}"


def model_dir_name(model_family: str, train_spec: str) -> str:
    return f"{model_family}_{model_suffix(train_spec)}"


def learned_model_names(train_specs: str) -> set[str]:
    specs = [normalize_train_spec(spec) for spec in split_list(train_specs, ("all",))]
    return {
        model_dir_name(model_family, spec)
        for model_family in LEARNED_MODEL_FAMILIES
        for spec in specs
    }


def configured_algorithm_names(configs: list[ExperimentConfig]) -> set[str]:
    names: set[str] = set()
    for config in configs:
        names.update(algorithm for algorithm in config.algorithms if algorithm in OFFSET_BASELINES)
        names.update(algorithm for algorithm in config.algorithms if algorithm in DIRECT_EXTERNAL_TOOLS)
        names.update(algorithm for algorithm in config.algorithms if algorithm in PROFILE_EXTERNAL_TOOLS)
        for model_family in LEARNED_MODEL_FAMILIES:
            if model_family not in config.algorithms:
                continue
            for train_spec in config.train_specs:
                names.add(model_dir_name(model_family, train_spec))
    return names


def configured_external_algorithms(configs: list[ExperimentConfig]) -> tuple[str, ...]:
    algorithms = []
    for config in configs:
        for algorithm in config.algorithms:
            if algorithm in (*DIRECT_EXTERNAL_TOOLS, *PROFILE_EXTERNAL_TOOLS) and algorithm not in algorithms:
                algorithms.append(algorithm)
    return tuple(algorithms)


def subset_compare_train_spec(configs: list[ExperimentConfig]) -> str | None:
    for config in configs:
        if any(algorithm in config.algorithms for algorithm in OFFSET_BASELINES):
            return config.compare_train_spec
    if configs:
        return configs[0].compare_train_spec
    return None


def infer_train_spec(model: str) -> str | None:
    if model.endswith("_all"):
        return "all"
    if "_k" not in model:
        return None
    _family, k = model.rsplit("_k", 1)
    if not k.isdigit():
        return None
    return k


def attach_train_read_fields(row: dict[str, object], train_reads_all: int | None = None) -> None:
    train_spec = infer_train_spec(str(row.get("model", row.get("algorithm", ""))))
    if train_spec == "all":
        row["train_read_limit"] = "all"
        if train_reads_all is not None:
            row["train_reads"] = train_reads_all
    elif train_spec is not None:
        row["train_read_limit"] = int(train_spec)
        row["train_reads"] = int(train_spec)


def train_read_count(train_path: Path | None) -> int | None:
    if train_path is None or not train_path.exists():
        return None
    return int(pd.read_parquet(train_path, columns=["read_id"])["read_id"].astype(str).nunique())


def is_learned_algorithm(algorithm: str) -> bool:
    return any(algorithm.startswith(f"{family}_") for family in LEARNED_MODEL_FAMILIES)


def ensure_common_test_subset(
    *,
    category_dir: Path,
    test_path: Path,
    requested_reads: str,
    seed: int,
) -> list[str]:
    subset_path = category_dir / "common_test_subset_reads.tsv"
    if subset_path.exists():
        return load_read_subset(subset_path)

    test_reads = pd.read_parquet(test_path, columns=["read_id"])
    read_ids = select_read_ids(test_reads, requested_reads, seed)
    metadata = read_subset_metadata(
        test_path=test_path,
        test_reads_requested=requested_reads,
        test_reads_selected=len(read_ids),
        selection_seed=seed,
        created_by="scripts/summarize_p_site_offset_prediction_evals.py",
        extra={
            "COMMON_TEST_READS": requested_reads,
            "COMMON_TEST_SEED": seed,
            "TEST_READS": os.environ.get("TEST_READS"),
            "test_table_rows": int(len(test_reads)),
            "test_table_unique_reads": int(test_reads["read_id"].astype(str).nunique()),
        },
    )
    write_read_subset(subset_path, read_ids, metadata)
    print(f"wrote {subset_path}")
    return read_ids


def metric_row(
    *,
    test: pd.DataFrame,
    algorithm_dir: Path,
    split: str,
    override: dict[str, object] | None = None,
    train_reads_all: int | None = None,
    target_read_ids: list[str] | None = None,
) -> dict[str, object] | None:
    algorithm = algorithm_dir.name
    pred_path = algorithm_dir / "candidate_predictions.parquet"
    if not pred_path.exists():
        return None

    if override is None:
        filter_ids = target_read_ids if target_read_ids is not None else test["read_id"]
        predictions = read_prediction_table(pred_path, filter_read_ids=filter_ids)
        if predictions.empty:
            print(f"skipping {algorithm}; empty predictions at {pred_path}")
            return None
        if target_read_ids is not None and is_learned_algorithm(algorithm):
            assert_prediction_read_set_matches_target(predictions, target_read_ids)
        metrics, _pred_read = evaluate_candidate_predictions(test, predictions)
        metrics.update(
            {
                "model": algorithm,
                "split": split,
                "predictions_path": str(pred_path),
                "metric_scope": "all_predicted_reads",
            }
        )
        metric_source = "candidate_predictions.parquet"
        _pred_read.to_parquet(algorithm_dir / "read_predictions.parquet", index=False)
    else:
        metrics = dict(override)
        metric_source = "candidate_predictions.parquet:model_subset"

    row: dict[str, object] = {
        "algorithm": algorithm,
        "metric_source": metric_source,
        **metrics,
    }
    attach_train_read_fields(row, train_reads_all)
    return row


def missing_common_metric_row(
    *,
    algorithm: str,
    pred_path: Path,
    split: str,
    target_reads: int,
    target_candidate_rows: int,
    train_reads_all: int | None,
    status: str,
    missing_reason: str,
) -> dict[str, object]:
    row: dict[str, object] = {
        "algorithm": algorithm,
        "metric_source": "candidate_predictions.parquet:common_test_subset",
        "exact_offset_accuracy": float("nan"),
        "within_1nt_accuracy": float("nan"),
        "frame_accuracy": float("nan"),
        "read_level_nll": float("nan"),
        "target_reads": float(target_reads),
        "predicted_target_reads": 0.0,
        "missing_target_reads": float(target_reads),
        "coverage": 0.0 if target_reads else float("nan"),
        "reads": float(target_reads),
        "candidate_rows": float(target_candidate_rows),
        "model": algorithm,
        "split": split,
        "predictions_path": str(pred_path),
        "metric_scope": "common_test_subset",
        "status": status,
        "missing_reason": missing_reason,
    }
    attach_train_read_fields(row, train_reads_all)
    return row


def common_metric_row(
    *,
    test: pd.DataFrame,
    algorithm: str,
    algorithm_dir: Path,
    split: str,
    target_read_ids: list[str],
    train_reads_all: int | None = None,
) -> dict[str, object]:
    pred_path = algorithm_dir / "candidate_predictions.parquet"
    target_set = set(target_read_ids)
    target_candidate_rows = int(test["read_id"].astype(str).isin(target_set).sum())
    if not pred_path.exists():
        return missing_common_metric_row(
            algorithm=algorithm,
            pred_path=pred_path,
            split=split,
            target_reads=len(target_read_ids),
            target_candidate_rows=target_candidate_rows,
            train_reads_all=train_reads_all,
            status="missing_predictions",
            missing_reason="candidate_predictions.parquet not found",
        )

    predictions = read_prediction_table(pred_path, filter_read_ids=target_read_ids)
    if predictions.empty:
        return missing_common_metric_row(
            algorithm=algorithm,
            pred_path=pred_path,
            split=split,
            target_reads=len(target_read_ids),
            target_candidate_rows=target_candidate_rows,
            train_reads_all=train_reads_all,
            status="empty_predictions",
            missing_reason="candidate_predictions.parquet is empty",
        )

    if is_learned_algorithm(algorithm):
        assert_prediction_read_set_matches_target(predictions, target_read_ids)
    metrics, pred_read = evaluate_targeted_candidate_predictions(
        test,
        predictions,
        target_read_ids,
        missing_policy="strict",
    )
    metrics.update(
        {
            "model": algorithm,
            "split": split,
            "predictions_path": str(pred_path),
            "metric_scope": "common_test_subset",
            "status": "ok",
            "missing_reason": "",
        }
    )
    pred_read.to_parquet(algorithm_dir / "read_predictions_common_test_subset.parquet", index=False)
    row: dict[str, object] = {
        "algorithm": algorithm,
        "metric_source": "candidate_predictions.parquet:common_test_subset",
        **metrics,
    }
    attach_train_read_fields(row, train_reads_all)
    return row


def collect_category_common_metrics(
    *,
    category_dir: Path,
    test_path: Path,
    split: str,
    target_read_ids: list[str],
    train_path: Path | None = None,
    allowed_algorithms: set[str] | None = None,
) -> pd.DataFrame:
    test = pd.read_parquet(test_path)
    all_train_reads = train_read_count(train_path)
    if allowed_algorithms is None:
        algorithms = sorted(path.name for path in category_dir.iterdir() if path.is_dir())
    else:
        algorithms = sorted(allowed_algorithms)
    rows = [
        common_metric_row(
            test=test,
            algorithm=algorithm,
            algorithm_dir=category_dir / algorithm,
            split=split,
            target_read_ids=target_read_ids,
            train_reads_all=all_train_reads,
        )
        for algorithm in algorithms
    ]
    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows)
    sort_cols = [col for col in ["train_reads", "algorithm"] if col in out.columns]
    if sort_cols:
        out = out.sort_values(sort_cols, na_position="last")
    return out


def write_category_common_summary(
    *,
    category_dir: Path,
    test_path: Path,
    split: str,
    target_read_ids: list[str],
    train_path: Path | None = None,
    allowed_algorithms: set[str] | None = None,
) -> pd.DataFrame:
    summary = collect_category_common_metrics(
        category_dir=category_dir,
        test_path=test_path,
        split=split,
        target_read_ids=target_read_ids,
        train_path=train_path,
        allowed_algorithms=allowed_algorithms,
    )
    if not summary.empty:
        write_csv(summary, category_dir / "summary.csv")
    return summary


def collect_category_metrics(
    *,
    category_dir: Path,
    test_path: Path,
    split: str,
    subset_metrics: dict[str, dict[str, object]] | None = None,
    train_path: Path | None = None,
    allowed_algorithms: set[str] | None = None,
    target_read_ids: list[str] | None = None,
) -> pd.DataFrame:
    test = pd.read_parquet(test_path)
    all_train_reads = train_read_count(train_path)
    rows: list[dict[str, object]] = []
    subset_metrics = subset_metrics or {}
    for algorithm_dir in sorted(path for path in category_dir.iterdir() if path.is_dir()):
        if allowed_algorithms is not None and algorithm_dir.name not in allowed_algorithms:
            continue
        row = metric_row(
            test=test,
            algorithm_dir=algorithm_dir,
            split=split,
            override=subset_metrics.get(algorithm_dir.name),
            train_reads_all=all_train_reads,
            target_read_ids=target_read_ids,
        )
        if row is None:
            continue
        rows.append(row)
    if not rows:
        return pd.DataFrame()

    out = pd.DataFrame(rows)
    sort_cols = [col for col in ["train_reads", "algorithm"] if col in out.columns]
    if sort_cols:
        out = out.sort_values(sort_cols, na_position="last")
    return out


def write_category_summary(
    *,
    category_dir: Path,
    test_path: Path,
    split: str,
    subset_metrics: dict[str, dict[str, object]] | None = None,
    train_path: Path | None = None,
    allowed_algorithms: set[str] | None = None,
    target_read_ids: list[str] | None = None,
) -> pd.DataFrame:
    summary = collect_category_metrics(
        category_dir=category_dir,
        test_path=test_path,
        split=split,
        subset_metrics=subset_metrics,
        train_path=train_path,
        allowed_algorithms=allowed_algorithms,
        target_read_ids=target_read_ids,
    )
    if not summary.empty:
        write_csv(summary, category_dir / "summary_diagnostic.csv")
    return summary


def clear_materialized_predictions(category_dir: Path, algorithm: str) -> None:
    algorithm_dir = category_dir / algorithm
    for filename in ("candidate_predictions.parquet", "read_predictions.parquet"):
        path = algorithm_dir / filename
        if path.exists():
            path.unlink()


def materialize_psite_predictions(
    *,
    baseline_root: Path,
    sample_table: Path,
    category_dir: Path,
    test_path: Path,
) -> None:
    bam_dir = baseline_root / "psite" / "tagged_bams"
    if not bam_dir.exists():
        clear_materialized_predictions(category_dir, "psite")
        return

    target = category_dir / "psite" / "candidate_predictions.parquet"
    target.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            sys.executable,
            "tools/psite/convert.py",
            "--sample-table",
            str(sample_table),
            "--bam-dir",
            str(bam_dir),
            "--test",
            str(test_path),
            "--output",
            str(target),
        ],
        check=True,
    )


def materialize_ribowaltz_predictions(
    *,
    baseline_root: Path,
    category_dir: Path,
    test_path: Path,
    offset_min: int = 8,
    offset_max: int = 20,
) -> None:
    offset_dir = baseline_root / "ribowaltz" / "offsets"
    offset_paths = sorted(offset_dir.glob("*.ribowaltz_offsets.tsv"))
    if not offset_paths:
        clear_materialized_predictions(category_dir, "ribowaltz")
        return

    offset_frames = []
    for path in offset_paths:
        frame = pd.read_csv(path, sep="\t")
        if frame.empty:
            continue
        offset_frames.append(
            frame[["sample", "length", "corrected_offset_from_5"]].rename(
                columns={
                    "sample": "sample_id",
                    "length": "read_length",
                    "corrected_offset_from_5": "candidate_offset",
                }
            )
        )
    if not offset_frames:
        clear_materialized_predictions(category_dir, "ribowaltz")
        return

    offsets = pd.concat(offset_frames, ignore_index=True)
    offsets["candidate_offset"] = pd.to_numeric(offsets["candidate_offset"], errors="coerce")
    offsets = offsets.dropna(subset=["candidate_offset"]).copy()
    offsets["candidate_offset"] = offsets["candidate_offset"].round().astype(int)
    offsets = offsets[offsets["candidate_offset"].between(offset_min, offset_max)].copy()
    if offsets.empty:
        clear_materialized_predictions(category_dir, "ribowaltz")
        return

    test = pd.read_parquet(test_path, columns=["read_id", "sample_id", "read_length"])
    reads = test.drop_duplicates("read_id")
    predictions = reads.merge(offsets, on=["sample_id", "read_length"], how="inner")
    predictions = predictions[["read_id", "candidate_offset"]].drop_duplicates()
    if predictions.empty:
        print(f"skipping ribowaltz materialization for {category_dir}; no matching reads")
        clear_materialized_predictions(category_dir, "ribowaltz")
        return
    predictions["p_offset"] = 1.0

    target = category_dir / "ribowaltz" / "candidate_predictions.parquet"
    target.parent.mkdir(parents=True, exist_ok=True)
    predictions.to_parquet(target, index=False)


def materialize_asite_ip_predictions(
    *,
    baseline_root: Path,
    category_dir: Path,
    test_path: Path,
) -> None:
    from tools.asite_ip.convert import (
        PREDICTION_COLUMNS,
        convert_offsets_to_predictions,
        load_offset_tables,
    )

    method_root = baseline_root / "asite_ip"
    target = category_dir / "asite_ip" / "candidate_predictions.parquet"
    result_paths = sorted(method_root.glob("*/output/Results_IP_algorithm.tab"))
    if not result_paths:
        print(f"skipping asite_ip materialization for {category_dir}; missing native Results_IP_algorithm.tab files")
        clear_materialized_predictions(category_dir, "asite_ip")
        return

    offsets = load_offset_tables(baseline_root=method_root)
    if offsets.empty:
        predictions = pd.DataFrame(columns=PREDICTION_COLUMNS)
    else:
        test = pd.read_parquet(
            test_path,
            columns=["read_id", "sample_id", "read_length", "five_prime_frame", "candidate_offset"],
        )
        predictions = convert_offsets_to_predictions(test, offsets)

    target.parent.mkdir(parents=True, exist_ok=True)
    predictions.to_parquet(target, index=False)
    if predictions.empty:
        print(f"materialized empty asite_ip predictions for {category_dir}")


def materialize_profile_baseline_predictions(
    *,
    method: str,
    baseline_root: Path,
    annotation_dir: Path,
    sample_table: Path,
    category_dir: Path,
    test_path: Path,
) -> None:
    if method not in PROFILE_EXTERNAL_TOOLS:
        raise ValueError(f"unsupported profile baseline: {method}")

    method_root = baseline_root / method
    target = category_dir / method / "candidate_predictions.parquet"
    if not method_root.exists():
        print(f"skipping {method} materialization for {category_dir}; missing {method_root}")
        clear_materialized_predictions(category_dir, method)
        return

    target.parent.mkdir(parents=True, exist_ok=True)
    if method == "scikit_ribo":
        subprocess.run(
            [
                sys.executable,
                "tools/scikit_ribo/convert.py",
                "--test",
                str(test_path),
                "--baseline-root",
                str(method_root),
                "--output",
                str(target),
            ],
            check=True,
        )
    else:
        subprocess.run(
            [
                sys.executable,
                "tools/shared/convert_profile_predictions.py",
                "--method",
                method,
                "--test",
                str(test_path),
                "--sample-table",
                str(sample_table),
                "--baseline-root",
                str(method_root),
                "--annotation-dir",
                str(annotation_dir),
                "--output",
                str(target),
            ],
            check=True,
        )


def materialize_leave_external_predictions(
    *,
    baseline_root: Path,
    annotation_dir: Path,
    sample_table: Path,
    category_dir: Path,
    test_path: Path,
    algorithms: tuple[str, ...],
) -> None:
    if "psite" in algorithms:
        materialize_psite_predictions(
            baseline_root=baseline_root,
            sample_table=sample_table,
            category_dir=category_dir,
            test_path=test_path,
        )
    if "ribodeblur" in algorithms:
        materialize_profile_baseline_predictions(
            method="ribodeblur",
            baseline_root=baseline_root,
            annotation_dir=annotation_dir,
            sample_table=sample_table,
            category_dir=category_dir,
            test_path=test_path,
        )
    if "riboba" in algorithms:
        materialize_profile_baseline_predictions(
            method="riboba",
            baseline_root=baseline_root,
            annotation_dir=annotation_dir,
            sample_table=sample_table,
            category_dir=category_dir,
            test_path=test_path,
        )
    if "scikit_ribo" in algorithms:
        materialize_profile_baseline_predictions(
            method="scikit_ribo",
            baseline_root=baseline_root,
            annotation_dir=annotation_dir,
            sample_table=sample_table,
            category_dir=category_dir,
            test_path=test_path,
        )
    if "ribowaltz" in algorithms:
        materialize_ribowaltz_predictions(
            baseline_root=baseline_root,
            category_dir=category_dir,
            test_path=test_path,
        )
    if "asite_ip" in algorithms:
        materialize_asite_ip_predictions(
            baseline_root=baseline_root,
            category_dir=category_dir,
            test_path=test_path,
        )


def evaluate_subset_baselines(
    *,
    test_path: Path,
    category_dir: Path,
    model_read_path: Path,
    split: str,
) -> dict[str, dict[str, object]]:
    if not model_read_path.exists():
        print(f"skipping subset baselines; missing {model_read_path}")
        return {}

    read_ids = set(pd.read_parquet(model_read_path, columns=["read_id"])["read_id"].astype(str))
    test = pd.read_parquet(test_path)
    test = test[test["read_id"].astype(str).isin(read_ids)].copy()
    rows: dict[str, dict[str, object]] = {}
    for model_name in ["fixed_offset", "read_length_offset"]:
        model_dir = category_dir / model_name
        pred_path = model_dir / "candidate_predictions.parquet"
        if not pred_path.exists():
            print(f"skipping {model_name} subset baseline; missing {pred_path}")
            continue

        predictions = read_prediction_table(pred_path, filter_read_ids=read_ids)
        metrics, pred_read = evaluate_candidate_predictions(test, predictions)
        metrics.update(
            {
                "model": model_name,
                "split": split,
                "predictions_path": str(pred_path),
                "subset_source": str(model_read_path),
                "metric_scope": "model_subset",
            }
        )
        pred_read.to_parquet(model_dir / "read_predictions_model_subset.parquet", index=False)
        rows[model_name] = metrics

    return rows


def write_leave_one_summary(result_root: Path, split_family: str, frames: list[pd.DataFrame]) -> None:
    if frames:
        write_csv(pd.concat(frames, ignore_index=True, sort=False), result_root / split_family / "summary_diagnostic.csv")


def write_leave_one_common_summary(result_root: Path, split_family: str, frames: list[pd.DataFrame]) -> None:
    if frames:
        write_csv(
            pd.concat(frames, ignore_index=True, sort=False),
            result_root / split_family / "summary.csv",
        )


def main() -> int:
    argparse.ArgumentParser(
        description="Summarize configured predictions on the strict common target reads.",
        epilog="Paths and subsets are configured through environment variables; see docs/reproduction.md.",
    ).parse_args()
    root = Path(os.environ.get("ROOT", "data/datasets/p_site_offset_prediction/landmark_riboseq"))
    result_root = Path(os.environ.get("RESULT_ROOT", "results/p_site_offset_prediction/landmark_riboseq"))
    baseline_root = Path(os.environ.get("NATIVE_OUTPUT_ROOT", os.environ.get("BASELINE_ROOT", str(result_root / "native_outputs"))))
    baseline_input_root = Path(os.environ.get("BASELINE_INPUT_ROOT", str(root / "views")))
    annotation_dir = Path(os.environ.get("ANNOT_DIR", str(baseline_input_root / "annotations")))
    sample_table = Path(os.environ.get("SAMPLE_TABLE", "data/datasets/p_site_offset_prediction/landmark_riboseq/samples.tsv"))
    eval_config = Path(os.environ.get("EVAL_CONFIG", str(DEFAULT_CONFIG)))
    configs = read_eval_config(eval_config)
    common_test_reads = str(os.environ.get("COMMON_TEST_READS", os.environ.get("TEST_READS", "5000")))
    common_test_seed = int(os.environ.get("COMMON_TEST_SEED", "20260518"))

    within_sample_dir = result_root / "within_sample"
    within_sample_train = root / "splits" / "within_sample" / "train.parquet"
    within_sample_test = root / "splits" / "within_sample" / "test.parquet"
    within_sample_configs = configs_for_experiment(configs, "within_sample")
    if within_sample_configs:
        within_sample_target_read_ids = ensure_common_test_subset(
            category_dir=within_sample_dir,
            test_path=within_sample_test,
            requested_reads=common_test_reads,
            seed=common_test_seed,
        )
        allowed_within_algorithms = configured_algorithm_names(within_sample_configs)
        write_category_common_summary(
            category_dir=within_sample_dir,
            test_path=within_sample_test,
            split="p_site_offset_prediction_landmark_riboseq_within_sample",
            target_read_ids=within_sample_target_read_ids,
            train_path=within_sample_train,
            allowed_algorithms=allowed_within_algorithms,
        )

    for split_family_dir in sorted((root / "splits").glob("leave_one_*")):
        if not split_family_dir.is_dir():
            continue
        split_family = split_family_dir.name
        leave_common_frames: list[pd.DataFrame] = []
        for heldout_dir in sorted(split_family_dir.glob("*")):
            if not heldout_dir.is_dir():
                continue
            heldout = heldout_dir.name
            experiment = f"{split_family}/{heldout}"
            heldout_configs = configs_for_experiment(configs, experiment)
            if not heldout_configs:
                continue
            category_dir = result_root / split_family / heldout
            heldout_train = heldout_dir / "train.parquet"
            heldout_test = heldout_dir / "test.parquet"
            split_name = f"p_site_offset_prediction_landmark_riboseq_{split_family}_{heldout}"
            heldout_target_read_ids = ensure_common_test_subset(
                category_dir=category_dir,
                test_path=heldout_test,
                requested_reads=common_test_reads,
                seed=common_test_seed,
            )
            materialize_leave_external_predictions(
                baseline_root=baseline_root,
                annotation_dir=annotation_dir,
                sample_table=sample_table,
                category_dir=category_dir,
                test_path=heldout_test,
                algorithms=configured_external_algorithms(heldout_configs),
            )
            allowed_heldout_algorithms = configured_algorithm_names(heldout_configs)
            common_summary = write_category_common_summary(
                category_dir=category_dir,
                test_path=heldout_test,
                split=split_name,
                target_read_ids=heldout_target_read_ids,
                train_path=heldout_train,
                allowed_algorithms=allowed_heldout_algorithms,
            )
            if not common_summary.empty:
                common_frame = common_summary.copy()
                common_frame.insert(0, "heldout", heldout)
                common_frame.insert(0, "split_family", split_family)
                leave_common_frames.append(common_frame)

        write_leave_one_common_summary(result_root, split_family, leave_common_frames)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
