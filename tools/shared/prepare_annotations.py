#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import pathlib
import re
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq


ATTR_RE = re.compile(r'([A-Za-z0-9_]+) "([^"]*)"')


def parse_attrs(raw: str) -> dict[str, str]:
    return {match.group(1): match.group(2) for match in ATTR_RE.finditer(raw)}


def parse_gtf(path: pathlib.Path) -> dict[str, dict]:
    transcripts: dict[str, dict] = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            contig, _source, feature, start, end, _score, strand, _phase, attrs_raw = fields
            attrs = parse_attrs(attrs_raw)
            transcript_id = attrs.get("transcript_id")
            if not transcript_id:
                continue
            row = transcripts.setdefault(
                transcript_id,
                {
                    "transcript": transcript_id,
                    "gene_id": attrs.get("gene_id", transcript_id),
                    "gene_name": attrs.get("gene_name", ""),
                    "gene_biotype": attrs.get("gene_biotype", ""),
                    "transcript_biotype": attrs.get("transcript_biotype", ""),
                    "chrom": contig,
                    "strand": strand,
                    "exons": [],
                    "cds": [],
                },
            )
            for key in ("gene_name", "gene_biotype", "transcript_biotype"):
                if attrs.get(key):
                    row[key] = attrs[key]
            if feature == "exon":
                row["exons"].append((int(start) - 1, int(end)))
            elif feature == "CDS":
                row["cds"].append((int(start) - 1, int(end)))
    return transcripts


def transcript_cds_bounds(row: dict) -> tuple[int, int] | None:
    exons = row["exons"]
    cds_segments = row["cds"]
    if not exons or not cds_segments:
        return None
    if row["strand"] == "+":
        ordered_exons = sorted(exons)
    else:
        ordered_exons = sorted(exons, reverse=True)

    cds_coords: list[tuple[int, int]] = []
    offset = 0
    for exon_start, exon_end in ordered_exons:
        exon_len = exon_end - exon_start
        for cds_start, cds_end in cds_segments:
            ov_start = max(exon_start, cds_start)
            ov_end = min(exon_end, cds_end)
            if ov_start >= ov_end:
                continue
            if row["strand"] == "+":
                tx_start = offset + (ov_start - exon_start)
                tx_end = offset + (ov_end - exon_start)
            else:
                tx_start = offset + (exon_end - ov_end)
                tx_end = offset + (exon_end - ov_start)
            cds_coords.append((tx_start, tx_end))
        offset += exon_len

    if not cds_coords:
        return None
    return min(start for start, _end in cds_coords), max(end for _start, end in cds_coords)


def build_annotation_rows(gtf_path: pathlib.Path, fasta_path: pathlib.Path) -> list[dict]:
    fasta_lengths = {record.id: len(record.seq) for record in SeqIO.parse(fasta_path, "fasta")}
    transcripts = parse_gtf(gtf_path)
    rows: list[dict] = []
    for transcript_id, row in sorted(transcripts.items()):
        tx_len = fasta_lengths.get(transcript_id)
        if tx_len is None:
            continue
        bounds = transcript_cds_bounds(row)
        if bounds is None:
            continue
        cds_start, cds_end = bounds
        if cds_start < 0 or cds_end > tx_len or cds_end <= cds_start:
            continue
        rows.append(
            {
                "transcript": transcript_id,
                "l_tr": tx_len,
                "l_utr5": cds_start,
                "l_cds": cds_end - cds_start,
                "l_utr3": tx_len - cds_end,
                "gene_id": row["gene_id"],
                "gene_name": row["gene_name"],
                "gene_biotype": row["gene_biotype"] or "protein_coding",
                "transcript_biotype": row["transcript_biotype"] or "protein_coding",
                "chrom": row["chrom"],
                "strand": row["strand"],
                "_genomic_start": min(start for start, _end in row["exons"]),
                "_genomic_end": max(end for _start, end in row["exons"]),
            }
        )
    return rows


def build_terminal_flanked_transcripts(
    rows: list[dict],
    transcript_fasta_path: pathlib.Path,
    genome_fasta_path: pathlib.Path,
    minimum_terminal_flank: int,
) -> tuple[list[dict], dict[str, str]]:
    """Add genomic terminal sequence when an annotated transcript lacks UTR context.

    Some compact annotations, notably the Ensembl yeast cDNA/GTF pair, begin at
    the start codon and end immediately after the stop codon. Start/stop
    metaprofile methods then have no reads that can span either boundary. This
    function preserves the annotated/spliced cDNA and adds only the missing
    amount of genomic sequence at its two ends.
    """
    def read_sequences(path: pathlib.Path) -> dict[str, str]:
        with path.open() as handle:
            return {
                record.id: str(record.seq).upper()
                for record in SeqIO.parse(handle, "fasta")
            }

    transcript_sequences = read_sequences(transcript_fasta_path)
    genome_sequences = read_sequences(genome_fasta_path)

    flanked_rows: list[dict] = []
    flanked_sequences: dict[str, str] = {}
    for row in rows:
        transcript_id = str(row["transcript"])
        transcript_sequence = transcript_sequences.get(transcript_id)
        if transcript_sequence is None:
            continue
        if len(transcript_sequence) != int(row["l_tr"]):
            raise ValueError(
                f"{transcript_id}: cDNA length {len(transcript_sequence)} "
                f"does not match annotation length {row['l_tr']}"
            )

        chrom = str(row["chrom"])
        genome_sequence = genome_sequences.get(chrom)
        if genome_sequence is None:
            raise ValueError(f"{transcript_id}: chromosome {chrom!r} is absent from genome FASTA")

        five_needed = max(0, minimum_terminal_flank - int(row["l_utr5"]))
        three_needed = max(0, minimum_terminal_flank - int(row["l_utr3"]))
        genomic_start = int(row["_genomic_start"])
        genomic_end = int(row["_genomic_end"])

        if row["strand"] == "+":
            five_sequence = genome_sequence[max(0, genomic_start - five_needed) : genomic_start]
            three_sequence = genome_sequence[
                genomic_end : min(len(genome_sequence), genomic_end + three_needed)
            ]
        elif row["strand"] == "-":
            five_sequence = str(
                Seq(
                    genome_sequence[
                        genomic_end : min(len(genome_sequence), genomic_end + five_needed)
                    ]
                ).reverse_complement()
            )
            three_sequence = str(
                Seq(
                    genome_sequence[max(0, genomic_start - three_needed) : genomic_start]
                ).reverse_complement()
            )
        else:
            raise ValueError(f"{transcript_id}: unsupported strand {row['strand']!r}")

        flanked_row = dict(row)
        flanked_row["l_utr5"] = int(row["l_utr5"]) + len(five_sequence)
        flanked_row["l_utr3"] = int(row["l_utr3"]) + len(three_sequence)
        flanked_row["l_tr"] = len(five_sequence) + len(transcript_sequence) + len(three_sequence)
        flanked_rows.append(flanked_row)
        flanked_sequences[transcript_id] = five_sequence + transcript_sequence + three_sequence

    return flanked_rows, flanked_sequences


def write_fasta(path: pathlib.Path, rows: list[dict], sequences: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for row in rows:
            transcript_id = str(row["transcript"])
            sequence = sequences[transcript_id]
            handle.write(f">{transcript_id}\n")
            for idx in range(0, len(sequence), 60):
                handle.write(sequence[idx : idx + 60] + "\n")


def write_terminal_flanked_gtf(
    source_path: pathlib.Path,
    output_path: pathlib.Path,
    original_rows: list[dict],
    flanked_rows: list[dict],
) -> None:
    """Extend terminal exon boundaries to match the generated flanked FASTA."""

    flanked_by_transcript = {
        str(row["transcript"]): row
        for row in flanked_rows
    }
    extension_by_transcript: dict[str, dict[str, int | str]] = {}
    gene_bounds: dict[str, tuple[int, int]] = {}

    for original in original_rows:
        transcript_id = str(original["transcript"])
        flanked = flanked_by_transcript.get(transcript_id)
        if flanked is None:
            continue

        five_added = int(flanked["l_utr5"]) - int(original["l_utr5"])
        three_added = int(flanked["l_utr3"]) - int(original["l_utr3"])
        if original["strand"] == "+":
            low_added, high_added = five_added, three_added
        else:
            low_added, high_added = three_added, five_added

        genomic_start = int(original["_genomic_start"])
        genomic_end = int(original["_genomic_end"])
        new_start = genomic_start - low_added
        new_end = genomic_end + high_added
        extension_by_transcript[transcript_id] = {
            "gene_id": str(original["gene_id"]),
            "genomic_start": genomic_start,
            "genomic_end": genomic_end,
            "low_added": low_added,
            "high_added": high_added,
        }

        gene_id = str(original["gene_id"])
        previous = gene_bounds.get(gene_id)
        if previous is None:
            gene_bounds[gene_id] = (new_start, new_end)
        else:
            gene_bounds[gene_id] = (
                min(previous[0], new_start),
                max(previous[1], new_end),
            )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with source_path.open() as source, output_path.open("w") as output:
        for raw_line in source:
            if not raw_line.strip() or raw_line.startswith("#"):
                output.write(raw_line)
                continue

            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 9:
                output.write(raw_line)
                continue

            feature = fields[2]
            attrs = parse_attrs(fields[8])
            transcript_id = attrs.get("transcript_id")
            gene_id = attrs.get("gene_id")
            start0 = int(fields[3]) - 1
            end = int(fields[4])

            extension = (
                extension_by_transcript.get(transcript_id)
                if transcript_id is not None
                else None
            )
            if feature in {"exon", "transcript"} and extension is not None:
                if start0 == int(extension["genomic_start"]):
                    start0 -= int(extension["low_added"])
                if end == int(extension["genomic_end"]):
                    end += int(extension["high_added"])
            elif feature == "gene" and gene_id in gene_bounds:
                gene_start, gene_end = gene_bounds[gene_id]
                start0 = min(start0, gene_start)
                end = max(end, gene_end)

            fields[3] = str(start0 + 1)
            fields[4] = str(end)
            output.write("\t".join(fields) + "\n")


def write_tsv(path: pathlib.Path, rows: list[dict], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row[column] for column in columns})


def write_asite_ip_start_index(path: pathlib.Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for row in rows:
            writer.writerow([row["transcript"], int(row["l_utr5"]), row["l_cds"]])


def write_ribodeblur_fasta(path: pathlib.Path, rows: list[dict], fasta_path: pathlib.Path) -> None:
    row_by_tx = {row["transcript"]: row for row in rows}
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for record in SeqIO.parse(fasta_path, "fasta"):
            row = row_by_tx.get(record.id)
            if row is None:
                continue
            cds_start = int(row["l_utr5"])
            cds_end = cds_start + int(row["l_cds"])
            handle.write(f">{record.id} |CDS:{cds_start}-{cds_end}|length:{row['l_tr']}\n")
            seq = str(record.seq)
            for idx in range(0, len(seq), 60):
                handle.write(seq[idx : idx + 60] + "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Prepare transcript annotation files for external p-site tool runs."
    )
    parser.add_argument(
        "--outdir",
        type=pathlib.Path,
        default=pathlib.Path("data/datasets/p_site_offset_prediction/landmark_riboseq/views/annotations"),
    )
    parser.add_argument(
        "--refs",
        nargs="+",
        default=["dmel", "scer"],
        choices=["dmel", "scer"],
    )
    parser.add_argument(
        "--transcript-outdir",
        type=pathlib.Path,
        default=pathlib.Path(
            "data/datasets/p_site_offset_prediction/landmark_riboseq/views/external_transcripts"
        ),
    )
    parser.add_argument(
        "--minimum-terminal-flank",
        type=int,
        default=50,
        help="Minimum annotated or added sequence on each side of the CDS.",
    )
    args = parser.parse_args(argv)
    if args.minimum_terminal_flank < 0:
        parser.error("--minimum-terminal-flank must be non-negative")

    for ref_key in args.refs:
        ref_dir = pathlib.Path("data/references") / ref_key
        original_rows = build_annotation_rows(ref_dir / "annotation.gtf", ref_dir / "cdna.fa")
        rows, transcript_sequences = build_terminal_flanked_transcripts(
            original_rows,
            ref_dir / "cdna.fa",
            ref_dir / "genome.fa",
            args.minimum_terminal_flank,
        )
        write_fasta(
            args.transcript_outdir / f"{ref_key}.fa",
            rows,
            transcript_sequences,
        )
        prefix = args.outdir / ref_key
        write_terminal_flanked_gtf(
            ref_dir / "annotation.gtf",
            prefix.with_suffix(".riboba_flanked.gtf"),
            original_rows,
            rows,
        )
        write_tsv(
            prefix.with_suffix(".ribowaltz_annotation.tsv"),
            rows,
            ["transcript", "l_tr", "l_utr5", "l_cds", "l_utr3"],
        )
        psite_rows = [
            {
                "gene_id": row["gene_id"],
                "tx_name": row["transcript"],
                "tx_len": row["l_tr"],
                "utr5_len": row["l_utr5"],
                "cds_len": row["l_cds"],
                "utr3_len": row["l_utr3"],
                "gene_biotype": row["gene_biotype"],
                "transcript_biotype": row["transcript_biotype"],
                "chrom": row["chrom"],
                "strand": row["strand"],
            }
            for row in rows
        ]
        write_tsv(
            prefix.with_suffix(".psite_txinfo.tsv"),
            psite_rows,
            [
                "gene_id",
                "tx_name",
                "tx_len",
                "utr5_len",
                "cds_len",
                "utr3_len",
                "gene_biotype",
                "transcript_biotype",
                "chrom",
                "strand",
            ],
        )
        write_asite_ip_start_index(
            prefix.with_suffix(".asite_ip_start_index.tsv"),
            rows,
        )
        write_ribodeblur_fasta(
            prefix.with_suffix(".ribodeblur.fa"),
            rows,
            args.transcript_outdir / f"{ref_key}.fa",
        )
        print(
            f"{ref_key}: wrote {len(rows)} protein-coding transcript annotations "
            f"with >= {args.minimum_terminal_flank}-nt terminal CDS context"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
