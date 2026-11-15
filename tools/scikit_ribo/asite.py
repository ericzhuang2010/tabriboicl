#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import importlib
import os
import pathlib
import re
import subprocess
import sys
import types
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV
from sklearn.model_selection import cross_val_score


class IdentitySelector:
    """Keep all features when RFECV can't run on the training set."""

    def __init__(self, n_features: int):
        self.n_features_ = n_features

    def transform(self, x: Any) -> Any:
        return x


def install_glmnet_shim() -> None:
    try:
        importlib.import_module("glmnet_py")
        return
    except ModuleNotFoundError:
        pass

    glmnet_py = types.ModuleType("glmnet_py")
    dataprocess = types.ModuleType("glmnet_py.dataprocess")
    glmnet_py.dataprocess = dataprocess
    sys.modules["glmnet_py"] = glmnet_py
    sys.modules["glmnet_py.dataprocess"] = dataprocess


def import_scikit_ribo(src_dir: pathlib.Path) -> tuple[Any, Any, Any, Any]:
    src_dir = src_dir.resolve()
    package_dir = src_dir / "scikit_ribo"
    if not package_dir.exists():
        raise FileNotFoundError(f"missing scikit-ribo source directory: {package_dir}")

    install_glmnet_shim()
    for path in (package_dir, src_dir):
        text = str(path)
        if text not in sys.path:
            sys.path.insert(0, text)

    gtf_preprocess = importlib.import_module("gtf_preprocess")
    bam_process = importlib.import_module("bam_process")
    asite_predict = importlib.import_module("asite_predict")
    return (
        gtf_preprocess.GtfPreProcess,
        bam_process.BamProcess,
        asite_predict.PredictAsite,
        asite_predict.VisualizeAsite,
    )


def required_index_paths(index_dir: pathlib.Path, prefix: str) -> list[pathlib.Path]:
    return [
        index_dir / f"{prefix}.start.bed",
        index_dir / f"{prefix}.cds.bed",
        index_dir / f"{prefix}.pos_ranges.txt",
        index_dir / f"{prefix}.nt_table.txt",
        index_dir / f"{prefix}.codons.df",
    ]


ATTR_RE = re.compile(r'([A-Za-z0-9_]+) "([^"]*)"')
REVCOMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def parse_gtf_attrs(raw: str) -> dict[str, str]:
    return {match.group(1): match.group(2) for match in ATTR_RE.finditer(raw)}


def reverse_complement(seq: str) -> str:
    return seq.translate(REVCOMP)[::-1].upper()


def read_cds_records(gtf_path: pathlib.Path) -> tuple[list[dict[str, Any]], dict[str, list[tuple[str, int, int, str]]]]:
    records: dict[str, dict[str, Any]] = {}
    starts: dict[str, list[tuple[str, int, int, str]]] = {}
    with gtf_path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            contig, _source, feature, start_s, end_s, _score, strand, _phase, attrs_raw = fields
            attrs = parse_gtf_attrs(attrs_raw)
            key = attrs.get("transcript_id") or attrs.get("gene_id")
            if not key:
                continue
            start = int(start_s) - 1
            end = int(end_s)
            gene_id = attrs.get("gene_id", key)
            if feature == "CDS":
                rec = records.setdefault(
                    key,
                    {
                        "name": key,
                        "gene_id": gene_id,
                        "contig": contig,
                        "strand": strand,
                        "segments": [],
                    },
                )
                if rec["contig"] != contig or rec["strand"] != strand:
                    continue
                rec["segments"].append((start, end))
            elif feature == "start_codon":
                starts.setdefault(key, []).append((contig, start, end, strand))
    out = [rec for rec in records.values() if rec["segments"]]
    out.sort(key=lambda rec: (str(rec["contig"]), min(start for start, _end in rec["segments"]), str(rec["name"])))
    return out, starts


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        prev_start, prev_end = merged[-1]
        if start <= prev_end:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))
    return merged


def write_lightweight_cds_bed(records: list[dict[str, Any]], path: pathlib.Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for rec in records:
            segments = sorted(rec["segments"])
            chrom_start = min(start for start, _end in segments)
            chrom_end = max(end for _start, end in segments)
            block_sizes = ",".join(str(end - start) for start, end in segments)
            block_starts = ",".join(str(start - chrom_start) for start, _end in segments)
            writer.writerow(
                [
                    rec["contig"],
                    chrom_start,
                    chrom_end,
                    rec["name"],
                    0,
                    rec["strand"],
                    chrom_start,
                    chrom_end,
                    0,
                    len(segments),
                    block_sizes,
                    block_starts,
                ]
            )


def write_lightweight_start_bed(
    records: list[dict[str, Any]],
    starts: dict[str, list[tuple[str, int, int, str]]],
    path: pathlib.Path,
) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for rec in records:
            explicit = starts.get(str(rec["name"]), [])
            if explicit:
                for contig, start, end, strand in explicit:
                    writer.writerow([contig, start, end, rec["name"], 0, strand])
                continue
            segments = sorted(rec["segments"])
            if rec["strand"] == "+":
                start, _end = segments[0]
                writer.writerow([rec["contig"], start, start + 3, rec["name"], 0, rec["strand"]])
            else:
                _start, end = segments[-1]
                writer.writerow([rec["contig"], end - 3, end, rec["name"], 0, rec["strand"]])


def position_ranges_for_record(rec: dict[str, Any], chrom_len: int) -> list[tuple[int, int, int]]:
    segments = sorted(rec["segments"])
    ranges: list[tuple[int, int, int]] = []
    if rec["strand"] == "+":
        phase = 0
        distance = 0
        for start, end in segments:
            ranges.append((start, end, phase))
            distance += end - start
            phase = distance % 3
    else:
        distance = 0
        for start, end in segments:
            distance += end - start
            phase = (3 - (distance % 3)) % 3
            ranges.append((start, end, phase))
    if not ranges:
        return []
    return [
        (max(0, ranges[0][0] - 24), ranges[0][0], 0),
        *ranges,
        (ranges[-1][1], min(chrom_len, ranges[-1][1] + 24), 0),
    ]


def write_lightweight_pos_ranges(records: list[dict[str, Any]], chrom_lengths: dict[str, int], path: pathlib.Path) -> None:
    with path.open("w") as handle:
        handle.write("#gene\tchr\tstrand\tpos_ranges\n")
        for rec in records:
            chrom_len = chrom_lengths.get(str(rec["contig"]), 0)
            ranges = position_ranges_for_record(rec, chrom_len)
            if not ranges:
                continue
            ranges_text = "|".join(f"{start},{end},{phase}" for start, end, phase in ranges if start < end)
            handle.write(f"{rec['name']}\t{rec['contig']}\t{rec['strand']}\t{ranges_text}\n")


def write_lightweight_nt_table(records: list[dict[str, Any]], fasta: Any, chrom_lengths: dict[str, int], path: pathlib.Path) -> None:
    intervals_by_chrom: dict[str, list[tuple[int, int]]] = {}
    for rec in records:
        chrom = str(rec["contig"])
        chrom_len = chrom_lengths.get(chrom)
        if chrom_len is None:
            continue
        start = max(0, min(seg_start for seg_start, _seg_end in rec["segments"]) - 30)
        end = min(chrom_len, max(seg_end for _seg_start, seg_end in rec["segments"]) + 30)
        if start < end:
            intervals_by_chrom.setdefault(chrom, []).append((start, end))
    with path.open("w") as handle:
        handle.write("chrom\tpos\tnt\n")
        for chrom in sorted(intervals_by_chrom):
            for start, end in merge_intervals(intervals_by_chrom[chrom]):
                seq = fasta.fetch(chrom, start, end).upper()
                for offset, nt in enumerate(seq):
                    handle.write(f"{chrom}\t{start + offset}\t{nt}\n")


def write_lightweight_codons(records: list[dict[str, Any]], fasta: Any, path: pathlib.Path) -> None:
    columns = ["chrom", "start", "end", "gene", "codon_idx", "gene_strand", "codon"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for rec in records:
            chrom = str(rec["contig"])
            if rec["strand"] == "+":
                ordered_segments = sorted(rec["segments"])
                positions = [pos for start, end in ordered_segments for pos in range(start, end)]
            else:
                ordered_segments = sorted(rec["segments"], reverse=True)
                positions = [pos for start, end in ordered_segments for pos in range(end - 1, start - 1, -1)]
            for idx in range(0, len(positions) - 2, 3):
                first = positions[idx]
                if rec["strand"] == "+":
                    start = first
                    end = first + 3
                    codon = fasta.fetch(chrom, start, end).upper()
                else:
                    start = first - 2
                    end = first + 1
                    codon = reverse_complement(fasta.fetch(chrom, start, end))
                writer.writerow(
                    {
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "gene": rec["name"],
                        "codon_idx": idx // 3,
                        "gene_strand": rec["strand"],
                        "codon": codon,
                    }
                )


def build_lightweight_indexes(args: argparse.Namespace) -> int:
    import pysam

    args.outdir.mkdir(parents=True, exist_ok=True)
    records, starts = read_cds_records(args.gtf)
    fasta = pysam.FastaFile(str(args.fasta))
    chrom_lengths = {name: length for name, length in zip(fasta.references, fasta.lengths)}
    records = [rec for rec in records if str(rec["contig"]) in chrom_lengths]
    if not records:
        raise ValueError(f"no CDS records from {args.gtf} matched contigs in {args.fasta}")

    write_lightweight_cds_bed(records, args.outdir / f"{args.prefix}.cds.bed")
    write_lightweight_start_bed(records, starts, args.outdir / f"{args.prefix}.start.bed")
    write_lightweight_pos_ranges(records, chrom_lengths, args.outdir / f"{args.prefix}.pos_ranges.txt")
    write_lightweight_nt_table(records, fasta, chrom_lengths, args.outdir / f"{args.prefix}.nt_table.txt")
    codons_bed = args.outdir / f"{args.prefix}.codons.bed"
    write_lightweight_codons(records, fasta, codons_bed)
    pd.read_csv(codons_bed, sep="\t").to_csv(args.outdir / f"{args.prefix}.codons.df", sep="\t", index=False)
    print(f"prepared lightweight scikit-ribo indexes at {args.outdir} from {len(records)} CDS transcripts")
    return 0


def safe_convert_gtf(worker: Any) -> None:
    """Convert the GTF as scikit-ribo does, skipping genes with invalid BED12 data."""

    import gffutils
    import pybedtools as pbt

    gtf_df = pd.read_table(worker.gtf, comment="#", header=None)
    gtf_df = gtf_df[(gtf_df.iloc[:, 2] != "gene") & (gtf_df.iloc[:, 2] != "transcript")]
    gtf_df = gtf_df[(gtf_df.iloc[:, 0] != "chrM") & (gtf_df.iloc[:, 0] != "chrMT") & (gtf_df.iloc[:, 0] != "MT")]
    gtf_dedup = pathlib.Path(worker.output) / f"{worker.prefix}.dedup.gtf"
    gtf_df.to_csv(path_or_buf=gtf_dedup, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)

    worker.bedtool = pbt.BedTool(str(gtf_dedup))
    worker.db = gffutils.create_db(str(gtf_dedup), ":memory:", force=True, verbose=False)
    candidate_gene_names = sorted(
        {
            gene["gene_id"][0]
            for gene in worker.db.features_of_type("CDS", order_by=("seqid", "start"))
            if "gene_id" in gene.attributes
        }
    )

    skipped: list[tuple[str, str]] = []
    worker.geneBed12s = []
    worker.geneNames = []
    for gene_name in candidate_gene_names:
        try:
            gene_bed12 = worker.db.bed12(gene_name, name_field="gene_id")
        except ValueError as exc:
            skipped.append((gene_name, str(exc)))
            continue
        row = gene_bed12.split("\t")
        worker.geneBed12s.append([row[0], int(row[1]), int(row[2])] + row[3:])
        worker.geneNames.append(gene_name)

    if not worker.geneBed12s:
        raise ValueError("no valid scikit-ribo BED12 CDS records could be built from the GTF")

    if skipped:
        print(
            f"[warning]\tskipped {len(skipped)} gene-level BED12 records that gffutils could not convert",
            file=sys.stderr,
        )
        for gene_name, message in skipped[:10]:
            print(f"[warning]\t{gene_name}: {message}", file=sys.stderr)
        if len(skipped) > 10:
            print(f"[warning]\t... {len(skipped) - 10} more skipped genes", file=sys.stderr)

    worker.geneBed12s = sorted(worker.geneBed12s, key=lambda x: (x[0], x[1], x[2]))
    with (pathlib.Path(worker.output) / f"{worker.prefix}.cds.bed").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerows(worker.geneBed12s)


def ensure_codons_df(index_dir: pathlib.Path, prefix: str) -> pathlib.Path:
    codons_bed = index_dir / f"{prefix}.codons.bed"
    codons_df = index_dir / f"{prefix}.codons.df"
    if not codons_bed.exists():
        raise FileNotFoundError(f"scikit-ribo build did not create {codons_bed}")
    frame = pd.read_csv(codons_bed, sep="\t")
    frame.to_csv(codons_df, sep="\t", index=False)
    return codons_df


def build_indexes(args: argparse.Namespace) -> int:
    if not args.force and all(path.exists() and path.stat().st_size > 0 for path in required_index_paths(args.outdir, args.prefix)):
        print(f"scikit-ribo indexes already prepared at {args.outdir}")
        return 0
    if args.builder == "lightweight":
        return build_lightweight_indexes(args)

    GtfPreProcess, _BamProcess, _PredictAsite, _VisualizeAsite = import_scikit_ribo(args.src_dir)
    args.outdir.mkdir(parents=True, exist_ok=True)

    worker = GtfPreProcess(str(args.gtf), str(args.fasta), args.prefix, str(args.outdir))
    try:
        worker.convertGtf()
    except ValueError as exc:
        print(f"[warning]\tscikit-ribo GTF conversion failed; retrying with invalid genes skipped: {exc}", file=sys.stderr)
        safe_convert_gtf(worker)
    worker.getChrLen()
    worker.getStartCodon()
    worker.getSeq()
    worker.getCodons()
    worker.getNts()
    worker.createCodonTable()
    ensure_codons_df(args.outdir, args.prefix)
    print(f"prepared scikit-ribo indexes at {args.outdir}")
    return 0



def load_target_positions(test_path: pathlib.Path | None, sample_id: str) -> dict[str, set[int]] | None:
    if test_path is None:
        return None
    sample = pd.read_parquet(
        test_path,
        columns=["contig", "candidate_ref_pos"],
        filters=[("sample_id", "==", str(sample_id))],
    )
    if sample.empty:
        return {}
    targets: dict[str, set[int]] = {}
    for row in sample[["contig", "candidate_ref_pos"]].drop_duplicates().itertuples(index=False):
        targets.setdefault(str(getattr(row, "contig")), set()).add(int(getattr(row, "candidate_ref_pos")))
    return targets


def filter_bed_to_target_positions(
    source_path: pathlib.Path,
    target_path: pathlib.Path,
    targets: dict[str, set[int]],
) -> int:
    kept = 0
    with source_path.open() as source, target_path.open("w") as target:
        for line in source:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            chrom_targets = targets.get(fields[0])
            if not chrom_targets:
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue
            if fields[5] == "+":
                possible_positions = range(start + 6, start + 16)
            elif fields[5] == "-":
                possible_positions = range(end - 16, end - 6)
            else:
                continue
            if any(position in chrom_targets for position in possible_positions):
                target.write(line)
                kept += 1
    if kept == 0:
        raise ValueError(f"no scikit-ribo reads can contribute to locked targets in {source_path}")
    return kept


def sort_bed_file(path: pathlib.Path) -> None:
    if not path.exists() or path.stat().st_size == 0:
        return
    sorted_path = path.with_name(f".{path.name}.sorted.tmp")
    sorted_path.unlink(missing_ok=True)
    env = os.environ.copy()
    env["LC_ALL"] = "C"
    try:
        subprocess.run(
            [
                "sort",
                "--parallel=1",
                "-S",
                "512M",
                "-T",
                str(path.parent),
                "-k1,1",
                "-k2,2n",
                "-k3,3n",
                "-o",
                str(sorted_path),
                str(path),
            ],
            check=True,
            env=env,
        )
        sorted_path.replace(path)
    except Exception:
        sorted_path.unlink(missing_ok=True)
        raise



def compatible_bam_process(base_cls: type) -> type:
    class CompatibleBamProcess(base_cls):  # type: ignore[misc, valid-type]
        def makeTrainingData(self) -> pd.DataFrame:  # noqa: N802 - keep scikit-ribo method name
            self.bedtool = self.bedtool.__class__(self.prefix + ".bed")
            training = self.bedtool.intersect(self.startCodons, wa=True, wb=True, sorted=True, s=True)
            print("[status]\tFinished intersecting the bedtool with start codons", file=sys.stderr)
            col_names = [
                "chrom",
                "start",
                "end",
                "name",
                "read_length",
                "strand",
                "sc_chrom",
                "sc_start",
                "sc_end",
                "gene",
                "sc_score",
                "gene_strand",
            ]
            training = training.to_dataframe(names=col_names)
            if not self.RelE:
                training["asite"] = np.where(
                    training["gene_strand"] == "+",
                    training["sc_start"] - training["start"] + 3,
                    training["end"] - training["sc_end"] + 3,
                )
            else:
                training["asite"] = np.where(
                    training["gene_strand"] == "+",
                    training["end"] - training["sc_start"] - 3,
                    training["sc_end"] - training["start"] - 3,
                )

            tmp_a = pd.merge(training, self.posOffsets, left_on=["chrom", "start"], right_on=["chrom", "pos"])
            tmp_b = pd.merge(training, self.negOffsets, left_on=["chrom", "end"], right_on=["chrom", "pos"])
            training = pd.concat([tmp_a, tmp_b])
            training.rename(columns={"offset": "5_offset"}, inplace=True)

            tmp_a = pd.merge(training, self.posOffsets, left_on=["chrom", "end"], right_on=["chrom", "pos"])
            tmp_b = pd.merge(training, self.negOffsets, left_on=["chrom", "start"], right_on=["chrom", "pos"])
            training = pd.concat([tmp_a, tmp_b])
            training.rename(columns={"offset": "3_offset"}, inplace=True)
            training = training.drop(columns=["pos", "pos_x", "pos_y"], errors="ignore")

            if not self.RelE:
                training = training[((training["asite"] >= 9) & (training["asite"] <= 18))]
                training = training[(training["asite"] >= training["read_length"] / 2 - 1)]
            else:
                training = training[((training["asite"] >= 1) & (training["asite"] <= 5))]

            self.nts = pd.read_table(self.ntsFn, header=0)
            training["pos_-1"] = np.where(training["gene_strand"] == "+", training["start"] - 1, training["end"])
            training["pos_0"] = np.where(training["gene_strand"] == "+", training["start"], training["end"] - 1)
            training["pos_n-1"] = np.where(training["gene_strand"] == "+", training["end"] - 1, training["start"])
            training["pos_n"] = np.where(training["gene_strand"] == "+", training["end"], training["start"] - 1)

            for source_col, nt_col in (
                ("pos_-1", "nt_-1"),
                ("pos_0", "nt_0"),
                ("pos_n-1", "nt_n-1"),
                ("pos_n", "nt_n"),
            ):
                training = pd.merge(training, self.nts, left_on=["chrom", source_col], right_on=["chrom", "pos"])
                training = training.drop(columns=[source_col, "pos"], errors="ignore")
                training.rename(columns={"nt": nt_col}, inplace=True)

            return training[
                ["asite", "read_length", "5_offset", "3_offset", "gene_strand", "nt_-1", "nt_0", "nt_n-1", "nt_n"]
            ]

    return CompatibleBamProcess

def compatible_predictor(base_cls: type, threads: int, random_state: int | None) -> type:
    class CompatiblePredictAsite(base_cls):  # type: ignore[misc, valid-type]
        def rfFit(self) -> None:  # noqa: N802 - keep scikit-ribo method name
            self.colNames = list(self.training.columns.values)
            self.colNames.remove("asite")
            dummy_training = pd.get_dummies(self.training[self.colNames]).astype(np.float32)
            self.dummyTrainingColNames = list(dummy_training.columns.values)
            self.X = sparse.csr_matrix(dummy_training.to_numpy(dtype=np.float32, copy=False))
            self.y = np.asarray(self.training["asite"])

            class_counts = pd.Series(self.y).value_counts()
            if len(class_counts) < 2:
                raise ValueError("scikit-ribo A-site training data has fewer than two classes")

            cpus = max(1, int(threads))
            self.clf = RandomForestClassifier(
                n_jobs=cpus,
                max_depth=10,
                min_samples_split=100,
                random_state=random_state,
            )
            self.clf = self.clf.fit(self.X, self.y)
            self.importance = self.clf.feature_importances_

            min_class_count = int(class_counts.min())
            cv_features = min(5, min_class_count)
            if cv_features >= 2:
                try:
                    self.selector = RFECV(self.clf, step=1, cv=cv_features)
                    self.selector = self.selector.fit(self.X, self.y)
                    selected_x = self.selector.transform(self.X)
                except Exception as exc:  # pragma: no cover - depends on sklearn/data edge cases
                    print(f"[warning]\tRFECV failed; using all features: {exc}", file=sys.stderr)
                    self.selector = IdentitySelector(self.X.shape[1])
                    selected_x = self.X
            else:
                self.selector = IdentitySelector(self.X.shape[1])
                selected_x = self.X

            print(
                "[result]\tOptimal number of features by recursive selection: %d" % self.selector.n_features_,
                flush=True,
            )
            self.reducedClf = RandomForestClassifier(
                n_jobs=cpus,
                max_depth=10,
                min_samples_split=100,
                random_state=random_state,
            )
            self.reducedClf = self.reducedClf.fit(selected_x, self.y)

            cv_scores = min(10, min_class_count)
            if cv_scores >= 2:
                scores = cross_val_score(self.reducedClf, selected_x, self.y, cv=cv_scores, scoring="accuracy")
                print("[result]\tAccuracy: %0.3f (+/- %0.3f)" % (scores.mean(), scores.std() * 2), flush=True)
            else:
                print("[result]\tAccuracy: skipped; not enough examples per A-site class", flush=True)


        def recoverAsite(self) -> pd.DataFrame:  # noqa: N802 - keep scikit-ribo method name
            if not self.RelE:
                self.cds.loc[self.cds["gene_strand"] == "+", "a_start"] = self.cds["start"] + self.cds["asite"]
                self.cds.loc[self.cds["gene_strand"] == "+", "a_end"] = self.cds["a_start"] + 3
                self.cds.loc[self.cds["gene_strand"] == "-", "a_end"] = self.cds["end"] - self.cds["asite"]
                self.cds.loc[self.cds["gene_strand"] == "-", "a_start"] = self.cds["a_end"] - 3
            else:
                self.cds.loc[self.cds["gene_strand"] == "+", "a_start"] = self.cds["end"] - self.cds["asite"]
                self.cds.loc[self.cds["gene_strand"] == "+", "a_end"] = self.cds["a_start"] + 3
                self.cds.loc[self.cds["gene_strand"] == "-", "a_end"] = self.cds["start"] + self.cds["asite"]
                self.cds.loc[self.cds["gene_strand"] == "-", "a_start"] = self.cds["a_end"] - 3

            self.cds.drop(["start", "end"], axis=1, inplace=True)
            reads = self.cds.groupby(["chrom", "a_start", "a_end", "strand"])
            reads = reads.size().reset_index(name="ribosome_count")
            cds_idx = pd.read_table(self.cdsIdxFn, header=0)
            orf_hits = pd.merge(
                cds_idx,
                reads,
                left_on=["chrom", "start", "end", "gene_strand"],
                right_on=["chrom", "a_start", "a_end", "strand"],
            )
            orf_hits = orf_hits[["chrom", "start", "end", "gene_strand", "ribosome_count"]].drop_duplicates()
            total_count = np.sum(reads.ribosome_count)
            num_orf_hits = np.sum(orf_hits.ribosome_count)
            orf_hits_percentage = num_orf_hits / total_count * 100 if total_count else 0.0
            print("[result]\tTotal number of Riboseq reads:", str(int(total_count)), flush=True)
            print("[result]\tNumber of reads in canonical ORF:", str(int(num_orf_hits)), flush=True)
            print("[result]\tPercentage of reads in canonical ORF:", str(round(orf_hits_percentage, 2)) + "%", flush=True)

            ribo_count = pd.merge(
                cds_idx,
                reads,
                how="left",
                left_on=["chrom", "start", "end", "gene_strand"],
                right_on=["chrom", "a_start", "a_end", "strand"],
            )
            ribo_count.drop(["a_start", "a_end", "strand"], axis=1, inplace=True)
            ribo_count["ribosome_count"] = ribo_count["ribosome_count"].fillna(0).astype(int)
            target_positions = getattr(self, "target_positions", None)
            if target_positions is not None:
                psite_coords = np.where(
                    ribo_count["gene_strand"] == "+",
                    ribo_count["start"].astype(int) - 3,
                    ribo_count["end"].astype(int) - 1 + 3,
                )
                keep = [
                    int(coord) in target_positions.get(str(chrom), set())
                    for chrom, coord in zip(ribo_count["chrom"], psite_coords)
                ]
                ribo_count = ribo_count.loc[keep].copy()
            ribo_count = ribo_count.sort_values(by=["chrom", "start", "end"])
            ribo_count.to_csv(
                path_or_buf=self.output + "/" + "riboseq_input.txt",
                sep="\t",
                header=True,
                index=False,
            )
            return ribo_count

        def rfPredict(self) -> None:  # noqa: N802 - keep scikit-ribo method name
            uniq_features = self.cds[self.colNames].groupby(self.colNames).size().reset_index(name="features")
            uniq_features = uniq_features[self.colNames]
            dummy_features = pd.get_dummies(uniq_features)
            dummy_features = dummy_features.reindex(columns=self.dummyTrainingColNames, fill_value=0).astype(np.float32)
            cds_x = sparse.csr_matrix(dummy_features.to_numpy(dtype=np.float32, copy=False))
            selected_cds_x = self.selector.transform(cds_x)
            uniq_features["asite"] = self.reducedClf.predict(selected_cds_x)
            self.cds = pd.merge(self.cds, uniq_features, how="left", on=self.colNames)

    return CompatiblePredictAsite


def run_asite(args: argparse.Namespace) -> int:
    _GtfPreProcess, BamProcess, PredictAsite, VisualizeAsite = import_scikit_ribo(args.src_dir)
    BamProcess = compatible_bam_process(BamProcess)
    riboseq_input = args.outdir / "riboseq_input.txt"
    if not args.force and riboseq_input.exists() and riboseq_input.stat().st_size > 0:
        print(f"scikit-ribo A-site profile already exists at {riboseq_input}")
        return 0

    for path in required_index_paths(args.index_dir, args.prefix):
        if not path.exists() or path.stat().st_size == 0:
            raise FileNotFoundError(f"missing scikit-ribo index file: {path}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    aln = BamProcess(
        bam=str(args.bam),
        mapq=args.mapq,
        minRL=args.min_read_len,
        maxRL=args.max_read_len,
        RelE=args.rele,
        output=str(args.outdir),
        startCodons=str(args.index_dir / f"{args.prefix}.start.bed"),
        cds=str(args.index_dir / f"{args.prefix}.cds.bed"),
        posRanges=str(args.index_dir / f"{args.prefix}.pos_ranges.txt"),
        nts=str(args.index_dir / f"{args.prefix}.nt_table.txt"),
    )
    aln.filterRegion()
    aln.posIndex()
    aln.filterBam()
    full_bed = args.outdir / "riboseq.bed"
    sort_bed_file(full_bed)
    target_positions = load_target_positions(args.target_test, args.outdir.name)
    training_data = aln.makeTrainingData()
    if target_positions is None:
        cds_data = aln.makeCdsData()
    else:
        target_bed = args.outdir / "riboseq.target.bed"
        kept_reads = filter_bed_to_target_positions(full_bed, target_bed, target_positions)
        print(f"[status]\tRetained {kept_reads} reads that can contribute to locked targets", file=sys.stderr)
        original_prefix = aln.prefix
        aln.prefix = str(target_bed.with_suffix(""))
        try:
            cds_data = aln.makeCdsData()
        finally:
            aln.prefix = original_prefix
    if training_data.empty:
        raise ValueError("scikit-ribo produced no start-codon training rows")
    if cds_data.empty:
        raise ValueError("scikit-ribo produced no CDS prediction rows")

    if args.plots:
        VisualizeAsite(training_data, args.rele, str(args.outdir)).plot()

    predictor_cls = compatible_predictor(PredictAsite, threads=args.threads, random_state=args.random_state)
    model = predictor_cls(training_data, cds_data, "rf", args.rele, args.prefix, str(args.outdir), str(args.index_dir))
    model.target_positions = target_positions
    model.rfFit()
    if args.plots:
        model.rfImportance()
        model.rocCurve()
    model.rfPredict()
    model.recoverAsite()
    if not riboseq_input.exists() or riboseq_input.stat().st_size == 0:
        raise FileNotFoundError(f"scikit-ribo did not create {riboseq_input}")
    print(f"wrote scikit-ribo A-site profile to {riboseq_input}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run the scikit-ribo A-site predictor for p-site offset prediction baselines.")
    parser.add_argument("--src-dir", type=pathlib.Path, default=pathlib.Path("tools/scikit_ribo/src"))
    subparsers = parser.add_subparsers(dest="command", required=True)

    build = subparsers.add_parser("build", help="Build scikit-ribo genomic indexes.")
    build.add_argument("--gtf", type=pathlib.Path, required=True)
    build.add_argument("--fasta", type=pathlib.Path, required=True)
    build.add_argument("--prefix", required=True)
    build.add_argument("--outdir", type=pathlib.Path, required=True)
    build.add_argument("--builder", choices=("lightweight", "scikit"), default="lightweight")
    build.add_argument("--force", action="store_true")
    build.set_defaults(func=build_indexes)

    run = subparsers.add_parser("run", help="Run scikit-ribo A-site prediction without TE modeling.")
    run.add_argument("--bam", type=pathlib.Path, required=True)
    run.add_argument("--index-dir", type=pathlib.Path, required=True)
    run.add_argument("--prefix", required=True)
    run.add_argument("--outdir", type=pathlib.Path, required=True)
    run.add_argument("--mapq", type=int, default=20)
    run.add_argument("--min-read-len", type=int, default=10)
    run.add_argument("--max-read-len", type=int, default=35)
    run.add_argument("--threads", type=int, default=4)
    run.add_argument("--random-state", type=int, default=1)
    run.add_argument("--target-test", type=pathlib.Path)
    run.add_argument("--rele", action="store_true")
    run.add_argument("--plots", action="store_true")
    run.add_argument("--force", action="store_true")
    run.set_defaults(func=run_asite)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
