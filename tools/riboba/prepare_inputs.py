#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import pathlib
import random
from dataclasses import dataclass

import pandas as pd
import pysam


SUPPORTED_RNASES = {
    "rnase_i": ("rnase-i", "reported"),
    "rnase-i": ("rnase-i", "reported"),
    "mnase": ("mnase", "reported"),
    "p1": ("p1", "reported"),
    # OTTR is a library method, not a nuclease. RiboBA has no OTTR or unknown
    # model, so use RNase-I and record that assumption.
    "ottr": ("rnase-i", "assumed_rnase_i"),
    "unknown": ("rnase-i", "assumed_rnase_i"),
    "": ("rnase-i", "assumed_rnase_i"),
}


@dataclass(frozen=True)
class FastqRecord:
    name: str
    sequence: str
    quality: str


@dataclass(frozen=True)
class ExtractionSummary:
    eligible_background_reads: int
    background_reads_written: int
    required_reads_requested: int
    required_reads_written: int
    total_reads_written: int


def fastq_record(read: pysam.AlignedSegment) -> FastqRecord | None:
    sequence = read.get_forward_sequence()
    qualities = read.get_forward_qualities()
    if not sequence or qualities is None:
        return None
    quality = "".join(chr(int(value) + 33) for value in qualities)
    return FastqRecord(str(read.query_name), sequence, quality)


def write_record(handle, record: FastqRecord) -> None:
    handle.write(f"@{record.name}\n{record.sequence}\n+\n{record.quality}\n")


def extract_fastq(
    *,
    bam_path: pathlib.Path,
    output_path: pathlib.Path,
    required_names: set[str],
    max_background_reads: int,
    seed: int,
) -> ExtractionSummary:
    """Keep required reads; use a seeded background sample when capped."""

    if max_background_reads < 0:
        raise ValueError("max_background_reads must be non-negative")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    rng = random.Random(seed)
    required_records: dict[str, FastqRecord] = {}
    background_records: list[FastqRecord] = []
    background_indexes: dict[str, int] = {}
    eligible_background_reads = 0

    with pysam.AlignmentFile(bam_path) as bam:
        if max_background_reads == 0:
            with output_path.open("w") as output:
                total = 0
                for read in bam.fetch(until_eof=True):
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue
                    if read.is_paired and read.is_read2:
                        continue
                    name = str(read.query_name)
                    record = fastq_record(read)
                    if record is None:
                        continue
                    write_record(output, record)
                    total += 1
                    if name in required_names:
                        required_records[name] = record
                    else:
                        eligible_background_reads += 1
            missing = required_names - set(required_records)
            if missing:
                preview = ", ".join(sorted(missing)[:5])
                raise ValueError(
                    f"{bam_path}: {len(missing)} required read(s) were absent; "
                    f"examples: {preview}"
                )
            return ExtractionSummary(
                eligible_background_reads=eligible_background_reads,
                background_reads_written=eligible_background_reads,
                required_reads_requested=len(required_names),
                required_reads_written=len(required_records),
                total_reads_written=total,
            )

        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.is_paired and read.is_read2:
                continue
            if not read.query_sequence or read.query_qualities is None:
                continue
            name = str(read.query_name)

            if name in required_names:
                if name in required_records:
                    continue
                record = fastq_record(read)
                if record is not None:
                    required_records[name] = record
                continue
            if name in background_indexes:
                continue

            eligible_background_reads += 1
            if len(background_records) < max_background_reads:
                record = fastq_record(read)
                if record is None:
                    continue
                background_indexes[name] = len(background_records)
                background_records.append(record)
                continue
            replacement = rng.randrange(eligible_background_reads)
            if replacement < max_background_reads:
                record = fastq_record(read)
                if record is None:
                    continue
                del background_indexes[background_records[replacement].name]
                background_records[replacement] = record
                background_indexes[name] = replacement

    missing = required_names - set(required_records)
    if missing:
        preview = ", ".join(sorted(missing)[:5])
        raise ValueError(
            f"{bam_path}: {len(missing)} required read(s) were absent; "
            f"examples: {preview}"
        )

    with output_path.open("w") as output:
        for record in background_records:
            write_record(output, record)
        for name in sorted(required_records):
            write_record(output, required_records[name])

    return ExtractionSummary(
        eligible_background_reads=eligible_background_reads,
        background_reads_written=len(background_records),
        required_reads_requested=len(required_names),
        required_reads_written=len(required_records),
        total_reads_written=len(background_records) + len(required_records),
    )


def read_sample_table(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"run", "ref_key"}
    missing = required - set(rows[0] if rows else ())
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    return rows


def read_test_metadata(
    path: pathlib.Path,
) -> tuple[dict[str, set[str]], dict[str, str]]:
    frame = pd.read_parquet(path, columns=["read_id", "sample_id", "nuclease"])
    unique = frame[["read_id", "sample_id", "nuclease"]].drop_duplicates()

    required_by_sample: dict[str, set[str]] = {}
    nuclease_by_sample: dict[str, str] = {}
    for sample_id, group in unique.groupby("sample_id", sort=False):
        sample = str(sample_id)
        prefix = f"{sample}|"
        names = {
            str(read_id)[len(prefix) :]
            for read_id in group["read_id"].astype(str)
            if str(read_id).startswith(prefix)
        }
        if len(names) != group["read_id"].astype(str).nunique():
            raise ValueError(f"{sample}: some read IDs do not begin with {prefix!r}")
        required_by_sample[sample] = names

        values = {
            str(value).strip().lower()
            for value in group["nuclease"].fillna("").tolist()
        }
        if len(values) != 1:
            raise ValueError(f"{sample}: conflicting nuclease values: {sorted(values)}")
        nuclease_by_sample[sample] = next(iter(values))
    return required_by_sample, nuclease_by_sample


def sample_seed(seed: int, sample_id: str) -> int:
    digest = hashlib.sha256(f"{seed}:{sample_id}".encode()).digest()
    return int.from_bytes(digest[:8], "little")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare per-organism FASTQs and RNase tables for native RiboBA. "
            "Every benchmark test read is retained, plus a reservoir-sampled "
            "background used to estimate RiboBA's bias model."
        )
    )
    parser.add_argument("--sample-table", type=pathlib.Path, required=True)
    parser.add_argument("--test", type=pathlib.Path, required=True)
    parser.add_argument("--genome-bam-dir", type=pathlib.Path, required=True)
    parser.add_argument("--outdir", type=pathlib.Path, required=True)
    parser.add_argument(
        "--max-background-reads",
        type=int,
        default=250_000,
        help="Background reads per sample; use 0 for every primary mapped read.",
    )
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--refs", nargs="*", default=None)
    args = parser.parse_args(argv)

    sample_rows = read_sample_table(args.sample_table)
    required_by_sample, nuclease_by_sample = read_test_metadata(args.test)
    selected_refs = set(args.refs) if args.refs else None

    args.outdir.mkdir(parents=True, exist_ok=True)
    status_rows: list[dict[str, object]] = []
    rnase_rows: dict[str, list[dict[str, str]]] = {}

    for row in sample_rows:
        sample = str(row["run"]).strip()
        ref_key = str(row["ref_key"]).strip()
        if selected_refs is not None and ref_key not in selected_refs:
            continue
        required_names = required_by_sample.get(sample)
        if required_names is None:
            continue

        raw_nuclease = nuclease_by_sample.get(sample, "")
        if raw_nuclease not in SUPPORTED_RNASES:
            raise ValueError(
                f"{sample}: unsupported nuclease metadata {raw_nuclease!r}; "
                "expected RNase I, MNase, P1, OTTR, or unknown"
            )
        rnase, rnase_source = SUPPORTED_RNASES[raw_nuclease]
        rnase_rows.setdefault(ref_key, []).append(
            {
                "sample_id": sample,
                "rnase": rnase,
                "rnase_source": rnase_source,
                "reported_nuclease": raw_nuclease or "unknown",
            }
        )

        bam_path = (
            args.genome_bam_dir
            / f"{sample}.Aligned.sortedByCoord.out.bam"
        )
        if not bam_path.is_file():
            raise FileNotFoundError(bam_path)
        output_path = args.outdir / "fastq" / ref_key / f"{sample}.fastq"
        summary = extract_fastq(
            bam_path=bam_path,
            output_path=output_path,
            required_names=required_names,
            max_background_reads=args.max_background_reads,
            seed=sample_seed(args.seed, sample),
        )
        status_rows.append(
            {
                "sample_id": sample,
                "ref_key": ref_key,
                "rnase": rnase,
                "rnase_source": rnase_source,
                **summary.__dict__,
            }
        )
        print(
            f"{sample}: wrote {summary.total_reads_written} reads "
            f"({summary.required_reads_written} required test reads)"
        )

    for ref_key, rows in rnase_rows.items():
        path = args.outdir / "refs" / ref_key / "rnase.tsv"
        path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(path, sep="\t", index=False)

    status = pd.DataFrame(status_rows)
    status.to_csv(args.outdir / "riboba_input_status.tsv", sep="\t", index=False)
    print(f"wrote {args.outdir / 'riboba_input_status.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
