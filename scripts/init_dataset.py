#!/usr/bin/env python3
"""Set up the data folders and copy the starting configs."""
from __future__ import annotations

import argparse
from pathlib import Path
import shutil

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def initialize(data_root: Path, samples: Path) -> Path:
    data_root = data_root.expanduser().resolve()
    dataset = data_root / "datasets/p_site_offset_prediction/landmark_riboseq"
    dataset.mkdir(parents=True, exist_ok=True)
    for source, target in (
        (samples, dataset / "samples.tsv"),
        (PROJECT_ROOT / "configs/dataset.yaml", dataset / "dataset.yaml"),
    ):
        if not target.exists():
            shutil.copy2(source, target)
    for directory in ("tables", "splits", "views", "references", "alignments"):
        (dataset / directory).mkdir(exist_ok=True)
    for name, target in (
        ("references/dmel", data_root / "references/dmel"),
        ("references/scer", data_root / "references/scer"),
        ("alignments/genome_bams", data_root / "intermediates/alignments/genome_bams"),
    ):
        target.mkdir(parents=True, exist_ok=True)
        link = dataset / name
        if not link.exists() and not link.is_symlink():
            link.symlink_to(target, target_is_directory=True)
    return dataset


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=PROJECT_ROOT / "data")
    parser.add_argument("--samples", type=Path, default=PROJECT_ROOT / "configs/samples.tsv")
    args = parser.parse_args()
    print(initialize(args.data_root, args.samples))


if __name__ == "__main__":
    main()
