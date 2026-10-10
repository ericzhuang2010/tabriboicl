#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import pathlib
import shutil
import sys
from collections import Counter, defaultdict

import pandas as pd
import pysam

PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from tools.shared.convert_profile_predictions import reference_pos_for_candidate, read_query_length


def read_sample_table(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_annotation(path: pathlib.Path, ref_key: str, mix_start: int) -> tuple[list[dict[str, object]], int]:
    frame = pd.read_csv(path, sep="\t")
    rows: list[dict[str, object]] = []
    cursor = mix_start
    for idx, row in enumerate(frame.itertuples(index=False), start=1):
        tx_len = int(row.l_tr)
        rows.append(
            {
                "ref_key": ref_key,
                "tx_id": f"{ref_key}:{idx}",
                "tx_name": str(row.transcript),
                "gene_id": str(row.transcript),
                "tx_len": tx_len,
                "cds_len": int(row.l_cds),
                "utr5_len": int(row.l_utr5),
                "utr3_len": int(row.l_utr3),
                "mix_start": cursor,
                "mix_cds_start": cursor + int(row.l_utr5),
                "mix_cds_end": cursor + int(row.l_utr5) + int(row.l_cds),
                "mix_end": cursor + tx_len,
            }
        )
        cursor += tx_len + 1
    return rows, cursor


def primary_mapped(read: pysam.AlignedSegment) -> bool:
    return not read.is_unmapped and not read.is_secondary and not read.is_supplementary


def ensure_link_or_copy(source: pathlib.Path, target: pathlib.Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() or target.is_symlink():
        target.unlink()
    rel_source = os.path.relpath(source, target.parent)
    try:
        target.symlink_to(rel_source)
    except OSError:
        shutil.copy2(source, target)


def choose_offsets(
    bam_path: pathlib.Path,
    tx_info: dict[str, dict[str, object]],
    offset_min: int,
    offset_max: int,
    default_offset: int | None = None,
) -> dict[int, int]:
    frame_counts: dict[int, Counter[int]] = defaultdict(Counter)
    with pysam.AlignmentFile(bam_path) as bam:
        for read in bam.fetch(until_eof=True):
            if not primary_mapped(read):
                continue
            info = tx_info.get(str(read.reference_name))
            if info is None:
                continue
            qwidth = read_query_length(read)
            if qwidth <= 0:
                continue
            utr5_len = int(info["utr5_len"])
            cds_len = int(info["cds_len"])
            for offset in range(offset_min, offset_max + 1):
                local_pos = reference_pos_for_candidate(read, offset)
                if local_pos is None:
                    continue
                cds_index = int(local_pos) - utr5_len
                if cds_index < 0 or cds_index >= cds_len:
                    continue
                frame_counts[qwidth][offset] += 1 if cds_index % 3 == 0 else 0
    offsets: dict[int, int] = {}
    if default_offset is None:
        default_offset = (offset_min + offset_max) // 2
    for qwidth, counts in frame_counts.items():
        if not counts:
            offsets[qwidth] = default_offset
            continue
        offsets[qwidth] = min(
            range(offset_min, offset_max + 1),
            key=lambda offset: (-counts[offset], abs(offset - default_offset), offset),
        )
    return offsets


def build_profile(
    *,
    sample: str,
    bam_path: pathlib.Path,
    tx_info: dict[str, dict[str, object]],
    offsets: dict[int, int],
    output_path: pathlib.Path,
    offset_min: int,
    offset_max: int,
    default_offset: int | None = None,
) -> tuple[int, int]:
    profile: Counter[tuple[int, int]] = Counter()
    reads_seen = 0
    if default_offset is None:
        default_offset = (offset_min + offset_max) // 2
    with pysam.AlignmentFile(bam_path) as bam:
        for read in bam.fetch(until_eof=True):
            if not primary_mapped(read):
                continue
            info = tx_info.get(str(read.reference_name))
            if info is None:
                continue
            qwidth = read_query_length(read)
            if qwidth <= 0:
                continue
            offset = offsets.get(qwidth, default_offset)
            local_pos = reference_pos_for_candidate(read, offset)
            if local_pos is None:
                continue
            utr5_len = int(info["utr5_len"])
            cds_len = int(info["cds_len"])
            cds_index = int(local_pos) - utr5_len
            if cds_index < 0 or cds_index >= cds_len:
                continue
            global_pos = int(info["mix_start"]) + int(local_pos)
            profile[(global_pos, qwidth)] += 1
            reads_seen += 1

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["pos", "qwidth", "weight"], delimiter="\t")
        writer.writeheader()
        for (pos, qwidth), weight in sorted(profile.items()):
            writer.writerow({"pos": pos, "qwidth": qwidth, "weight": weight})
    return reads_seen, len(profile)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build RiboBA-compatible profile fallback artifacts.")
    parser.add_argument("--sample-table", type=pathlib.Path, required=True)
    parser.add_argument("--annotation-dir", type=pathlib.Path, required=True)
    parser.add_argument("--bam-dir", type=pathlib.Path, required=True)
    parser.add_argument("--outdir", type=pathlib.Path, required=True)
    parser.add_argument("--offset-min", type=int, default=8)
    parser.add_argument("--offset-max", type=int, default=20)
    parser.add_argument("--offset-mode", choices=("fixed", "frame"), default="fixed")
    parser.add_argument("--fixed-offset", type=int, default=13)
    args = parser.parse_args(argv)

    samples = read_sample_table(args.sample_table)
    ref_keys = sorted({row["ref_key"] for row in samples if row.get("ref_key")})

    tx_rows: list[dict[str, object]] = []
    cursor = 0
    for ref_key in ref_keys:
        rows, cursor = read_annotation(
            args.annotation_dir / f"{ref_key}.ribowaltz_annotation.tsv",
            ref_key,
            cursor,
        )
        tx_rows.extend(rows)
    tx_frame = pd.DataFrame(tx_rows)
    args.outdir.mkdir(parents=True, exist_ok=True)
    tx_frame.to_csv(args.outdir / "riboba_tx_map.tsv", sep="\t", index=False)
    tx_info = {str(row["tx_name"]): row for row in tx_rows}

    bam_outdir = args.outdir / "prepare_bam" / "bam"
    downstream = args.outdir / "downstream"
    status_path = args.outdir / "riboba_status.tsv"
    with status_path.open("w", newline="") as status_handle:
        status = csv.DictWriter(status_handle, fieldnames=["sample_id", "status", "message"], delimiter="\t")
        status.writeheader()
        for row in samples:
            sample = row["run"]
            source_bam = args.bam_dir / f"{sample}.transcriptome.bam"
            source_bai = args.bam_dir / f"{sample}.transcriptome.bam.bai"
            if not source_bam.exists():
                status.writerow({"sample_id": sample, "status": "failed", "message": f"missing {source_bam}"})
                continue
            target_bam = bam_outdir / f"{sample}_cds_txsorted.bam"
            target_bai = bam_outdir / f"{sample}_cds_txsorted.bam.bai"
            ensure_link_or_copy(source_bam, target_bam)
            if source_bai.exists():
                ensure_link_or_copy(source_bai, target_bai)

            if args.offset_mode == "fixed":
                offsets = {}
                default_offset = args.fixed_offset
            else:
                offsets = choose_offsets(source_bam, tx_info, args.offset_min, args.offset_max)
                default_offset = None
            reads_seen, profile_rows = build_profile(
                sample=sample,
                bam_path=source_bam,
                tx_info=tx_info,
                offsets=offsets,
                output_path=downstream / sample / f"{sample}.riboba_corrected_profile.tsv",
                offset_min=args.offset_min,
                offset_max=args.offset_max,
                default_offset=default_offset,
            )
            status.writerow(
                {
                    "sample_id": sample,
                    "status": "ok" if profile_rows else "failed",
                    "message": f"profile_rows={profile_rows}; reads={reads_seen}; read_lengths={len(offsets)}; offset_mode={args.offset_mode}; fixed_offset={args.fixed_offset if args.offset_mode == 'fixed' else 'NA'}",
                }
            )
    pipeline_status = args.outdir / "riboba_pipeline_status.tsv"
    pipeline_status.write_text(f"stage\tstatus\tmessage\npipeline\tok\tprofile fallback artifacts generated; offset_mode={args.offset_mode}; fixed_offset={args.fixed_offset}\n")
    print(f"wrote {args.outdir / 'riboba_tx_map.tsv'}")
    print(f"wrote profiles under {downstream}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
