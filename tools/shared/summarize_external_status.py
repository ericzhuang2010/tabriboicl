#!/usr/bin/env python3
"""Summarize p-site offset prediction external-tool within-sample outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    sep = "," if path.suffix == ".csv" else "\t"
    return pd.read_csv(path, sep=sep)


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)


def write_tool_statuses(root: Path, metrics_dir: Path) -> None:
    for tool in ("psite", "ribowaltz", "ribodeblur", "scikit_ribo", "asite_ip"):
        status = read_table(root / tool / f"{tool}_status.tsv")
        if not status.empty:
            write_csv(status, metrics_dir / tool / "status.csv")

    riboba_status = read_table(root / "riboba" / "riboba_status.tsv")
    if not riboba_status.empty:
        write_csv(riboba_status, metrics_dir / "riboba" / "status.csv")

    riboba_pipeline = read_table(root / "riboba" / "riboba_pipeline_status.tsv")
    if not riboba_pipeline.empty:
        write_csv(riboba_pipeline, metrics_dir / "riboba" / "pipeline_status.csv")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--baseline-root",
        type=Path,
        default=Path("results/p_site_offset_prediction/landmark_riboseq/native_outputs"),
    )
    parser.add_argument(
        "--metrics-dir",
        type=Path,
        default=Path("results/p_site_offset_prediction/landmark_riboseq/within_sample"),
    )
    args = parser.parse_args()

    args.metrics_dir.mkdir(parents=True, exist_ok=True)
    write_tool_statuses(args.baseline_root, args.metrics_dir)


if __name__ == "__main__":
    main()
