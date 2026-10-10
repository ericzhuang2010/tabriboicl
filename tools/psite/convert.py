#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import pathlib
import sqlite3

import pyarrow as pa
import pyarrow.parquet as pq
import pysam


SCHEMA = pa.schema(
    [
        ("read_id", pa.string()),
        ("candidate_offset", pa.int64()),
        ("p_offset", pa.float64()),
    ]
)


def read_sample_table(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def write_parquet_batch(writer: pq.ParquetWriter, rows: list[tuple[str, int, float]]) -> None:
    table = pa.Table.from_arrays(
        [
            pa.array([row[0] for row in rows], type=pa.string()),
            pa.array([row[1] for row in rows], type=pa.int64()),
            pa.array([row[2] for row in rows], type=pa.float64()),
        ],
        schema=SCHEMA,
    )
    writer.write_table(table)


def flush_inserts(conn: sqlite3.Connection, rows: list[tuple[str, int, float]]) -> None:
    conn.executemany(
        "INSERT OR IGNORE INTO predictions(read_id, candidate_offset, p_offset) VALUES (?, ?, ?)",
        rows,
    )
    conn.commit()
    rows.clear()


def open_prediction_db(path: pathlib.Path) -> sqlite3.Connection:
    if path.exists():
        path.unlink()
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    conn.execute("PRAGMA temp_store=FILE")
    conn.execute("PRAGMA cache_size=-65536")
    conn.execute(
        """
        CREATE TABLE predictions (
            read_id TEXT NOT NULL,
            candidate_offset INTEGER NOT NULL,
            p_offset REAL NOT NULL,
            PRIMARY KEY (read_id, candidate_offset)
        ) WITHOUT ROWID
        """
    )
    return conn


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Convert PSite PS-tagged BAMs to B1 candidate predictions.")
    parser.add_argument("--sample-table", type=pathlib.Path, required=True)
    parser.add_argument("--bam-dir", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--test", type=pathlib.Path)
    parser.add_argument("--offset-min", type=int, default=8)
    parser.add_argument("--offset-max", type=int, default=20)
    parser.add_argument("--batch-size", type=int, default=100_000)
    args = parser.parse_args(argv)

    target_read_ids: set[str] | None = None
    target_sample_ids: set[str] | None = None
    if args.test is not None:
        target_read_ids = set(
            pq.read_table(args.test, columns=["read_id"])["read_id"].to_pylist()
        )
        target_sample_ids = {
            str(read_id).split("|", 1)[0]
            for read_id in target_read_ids
        }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    tmp_output = args.output.with_name(f"{args.output.name}.tmp")
    db_path = args.output.with_name(f"{args.output.name}.sqlite.tmp")
    if tmp_output.exists():
        tmp_output.unlink()

    conn = open_prediction_db(db_path)
    insert_rows: list[tuple[str, int, float]] = []
    try:
        for sample in read_sample_table(args.sample_table):
            sample_id = sample["run"]
            if target_sample_ids is not None and sample_id not in target_sample_ids:
                continue
            bam_path = args.bam_dir / f"{sample_id}.psite.tagged.bam"
            if not bam_path.exists():
                continue
            with pysam.AlignmentFile(bam_path) as bam:
                for read in bam.fetch(until_eof=True):
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue
                    if not read.has_tag("PS"):
                        continue
                    offset = int(read.get_tag("PS"))
                    if offset < args.offset_min or offset > args.offset_max:
                        continue
                    read_id = f"{sample_id}|{read.query_name}"
                    if target_read_ids is not None and read_id not in target_read_ids:
                        continue
                    insert_rows.append((read_id, offset, 1.0))
                    if len(insert_rows) >= args.batch_size:
                        flush_inserts(conn, insert_rows)
        if insert_rows:
            flush_inserts(conn, insert_rows)

        writer = pq.ParquetWriter(tmp_output, SCHEMA)
        row_count = 0
        try:
            cursor = conn.execute("SELECT read_id, candidate_offset, p_offset FROM predictions")
            while True:
                parquet_rows = cursor.fetchmany(args.batch_size)
                if not parquet_rows:
                    break
                write_parquet_batch(writer, parquet_rows)
                row_count += len(parquet_rows)
        finally:
            writer.close()
    finally:
        conn.close()

    tmp_output.replace(args.output)
    db_path.unlink(missing_ok=True)
    scope = "benchmark test reads" if target_read_ids is not None else "all tagged BAM reads"
    print(
        f"wrote {row_count} PSite predictions for {scope} to {args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
