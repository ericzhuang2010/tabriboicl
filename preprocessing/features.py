from __future__ import annotations

import argparse
import bisect
import csv
import dataclasses
import pathlib
import re
from collections import Counter, defaultdict
from typing import Iterable

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pysam
from Bio.Seq import Seq


# P-site offset extraction turns each aligned ribosome-profiling read into a small set of
# possible P-site offsets. The default offset window is intentionally broad:
# downstream models decide which offset looks most plausible for a read.
DEFAULT_OFFSETS = tuple(range(8, 21))

# Missing integer annotations are stored as -1 so the Parquet columns can keep
# fixed integer types. Empty strings play the same role for string annotations.
MISSING_INT = -1


# Candidate rows are the model-facing table: one row per read-offset candidate.
# Keep this order in sync with CANDIDATE_SCHEMA so DataFrame construction writes
# stable Parquet columns.
CANDIDATE_COLUMNS = [
    "read_id",  # Sample-prefixed read name; unique within the emitted dataset.
    "sample_id",  # Sample/run identifier copied from the metadata table.
    "study_id",  # Source study identifier copied from the metadata table.
    "organism",  # Organism name copied from the metadata table.
    "protocol",  # Ribosome-profiling protocol copied from the metadata table.
    "nuclease",  # Nuclease/enzyme metadata for the library.
    "library_method",  # Library preparation method metadata.
    "library_name",  # Library name metadata.
    "read_length",  # Query/read length in nucleotides.
    "strand",  # Strand-aware read orientation, "+" or "-".
    "contig",  # Reference sequence/chromosome containing the alignment, e.g 2L, 3R
    "mapq",  # Alignment mapping quality.
    "five_prime_pos",  # Strand-aware genomic coordinate of the read 5-prime end.
    "three_prime_pos",  # Strand-aware genomic coordinate of the read 3-prime end.
    "candidate_offset",  # Tested 5-prime-to-P-site offset in read-oriented bases.
    "candidate_ref_pos",  # Genomic coordinate of the candidate P-site.
    "candidate_frame",  # CDS frame at the candidate P-site, or -1 outside CDS.
    "five_prime_frame",  # CDS frame at the read 5-prime end, or -1 outside CDS.
    "three_prime_frame",  # CDS frame at the read 3-prime end, or -1 outside CDS.
    "distance_to_cds_start",  # Candidate CDS-coordinate distance from CDS start.
    "distance_to_cds_stop",  # Candidate distance to the final CDS base.
    "cds_percentile",  # Candidate CDS position scaled from start 0.0 to end 1.0.
    "transcript_id",  # Selected transcript overlapping the candidate P-site.
    "gene_id",  # Gene ID for the selected candidate transcript.
    "gene_name",  # Gene name for the selected candidate transcript.
    "landmark",  # Candidate landmark label such as cds_start or stop-last-codon.
    "five_prime_kmer",  # Strand-oriented k-mer centered on the read 5-prime end.
    "three_prime_kmer",  # Strand-oriented k-mer centered on the read 3-prime end.
    "candidate_site_kmer",  # Strand-oriented k-mer centered on the candidate P-site.
    "candidate_p_codon",  # Strand-oriented codon starting at the candidate P-site.
    "candidate_a_codon",  # Strand-oriented A-site codon one codon downstream.
    "candidate_e_codon",  # Strand-oriented E-site codon one codon upstream.
    "local_coverage_pm15",  # Mean read coverage in the local candidate window.
    "local_5p_end_density",  # Count of 5-prime read ends in the local window.
    "local_3p_end_density",  # Count of 3-prime read ends in the local window.
    "label",  # Silver label: 1 winning landmark offset, 0 other offsets, else null.
    "label_tier",  # Label provenance or ambiguity reason for the candidate row.
]


# Read rows are the compact per-read table. They are written for every usable
# read, even when candidate rows are not written in "labeled" output mode.
READ_COLUMNS = [
    "read_id",  # Unique read key: "{sample_id}|{BAM query_name}".
    "sample_id",  # SRA run ID, for example SRR1039768.
    "study_id",  # Study/project accession from sample table or RunInfo metadata.
    "organism",  # Organism name attached to this run.
    "protocol",  # Protocol label used for grouping, usually matching nuclease.
    "nuclease",  # Inferred or explicit nuclease/library type for this run.
    "library_method",  # Library selection/strategy from RunInfo or sample table.
    "library_name",  # SRA library/sample name when available.
    "read_length",  # Query/read length inferred from the BAM alignment.
    "strand",  # Read alignment strand: '+' for forward, '-' for reverse.
    "contig",  # Reference sequence/chromosome name where the read aligns.
    "mapq",  # BAM mapping quality for this primary alignment.
    "five_prime_pos",  # Strand-aware genomic coordinate of the read 5-prime end.
    "three_prime_pos",  # Strand-aware genomic coordinate of the read 3-prime end.
    "five_prime_frame",  # CDS frame at the read 5-prime end, or -1 outside CDS.
    "three_prime_frame",  # CDS frame at the read 3-prime end, or -1 outside CDS.
    "silver_offset",  # Winning candidate offset if exactly one landmark hit exists.
    "silver_landmark",  # Landmark hit by silver_offset, or ambiguity/missing marker.
    "label_tier",  # Per-read label status: labeled landmark, ambiguous, or empty.
]


# Explicit schemas keep Parquet types small and predictable across samples.
# This matters because the per-sample files are later concatenated directly.
CANDIDATE_SCHEMA = pa.schema(
    [
        ("read_id", pa.string()),
        ("sample_id", pa.string()),
        ("study_id", pa.string()),
        ("organism", pa.string()),
        ("protocol", pa.string()),
        ("nuclease", pa.string()),
        ("library_method", pa.string()),
        ("library_name", pa.string()),
        ("read_length", pa.int16()),
        ("strand", pa.string()),
        ("contig", pa.string()),
        ("mapq", pa.int16()),
        ("five_prime_pos", pa.int32()),
        ("three_prime_pos", pa.int32()),
        ("candidate_offset", pa.int16()),
        ("candidate_ref_pos", pa.int32()),
        ("candidate_frame", pa.int8()),
        ("five_prime_frame", pa.int8()),
        ("three_prime_frame", pa.int8()),
        ("distance_to_cds_start", pa.int32()),
        ("distance_to_cds_stop", pa.int32()),
        ("cds_percentile", pa.float32()),
        ("transcript_id", pa.string()),
        ("gene_id", pa.string()),
        ("gene_name", pa.string()),
        ("landmark", pa.string()),
        ("five_prime_kmer", pa.string()),
        ("three_prime_kmer", pa.string()),
        ("candidate_site_kmer", pa.string()),
        ("candidate_p_codon", pa.string()),
        ("candidate_a_codon", pa.string()),
        ("candidate_e_codon", pa.string()),
        ("local_coverage_pm15", pa.float32()),
        ("local_5p_end_density", pa.int32()),
        ("local_3p_end_density", pa.int32()),
        ("label", pa.int8()),
        ("label_tier", pa.string()),
    ]
)


# The read table stores only read-level metadata plus the silver landmark label.
READ_SCHEMA = pa.schema(
    [
        ("read_id", pa.string()),
        ("sample_id", pa.string()),
        ("study_id", pa.string()),
        ("organism", pa.string()),
        ("protocol", pa.string()),
        ("nuclease", pa.string()),
        ("library_method", pa.string()),
        ("library_name", pa.string()),
        ("read_length", pa.int16()),
        ("strand", pa.string()),
        ("contig", pa.string()),
        ("mapq", pa.int16()),
        ("five_prime_pos", pa.int32()),
        ("three_prime_pos", pa.int32()),
        ("five_prime_frame", pa.int8()),
        ("three_prime_frame", pa.int8()),
        ("silver_offset", pa.int16()),
        ("silver_landmark", pa.string()),
        ("label_tier", pa.string()),
    ]
)


@dataclasses.dataclass(frozen=True)
class SampleContext:
    """Metadata attached to every row emitted for one SRA run."""

    sample_id: str
    organism: str
    ref_key: str
    study_id: str
    protocol: str
    nuclease: str
    library_method: str
    library_name: str


@dataclasses.dataclass
class CdsSegment:
    """One half-open genomic CDS interval: [start, end)."""

    start: int
    end: int


@dataclasses.dataclass
class TranscriptModel:
    """CDS-only transcript model used for coordinate and frame annotation.

    GTF stores CDS features as genomic segments. This class keeps those segments
    in transcript order so a genomic coordinate can be converted into a CDS
    coordinate, regardless of transcript strand.
    """

    transcript_id: str
    gene_id: str
    gene_name: str
    contig: str
    strand: str
    segments: list[CdsSegment]

    def __post_init__(self) -> None:
        # Transcript order is genomic ascending on the plus strand and genomic
        # descending on the minus strand.
        if self.strand == "+":
            self.segments.sort(key=lambda segment: segment.start)
        else:
            self.segments.sort(key=lambda segment: segment.start, reverse=True)

        # prefix_lengths[i] is the CDS length before segment i in transcript
        # orientation. It lets pos_to_cds_coord return a transcript-relative
        # coordinate without flattening the whole CDS into per-base arrays.
        self.prefix_lengths: list[int] = []
        total = 0
        for segment in self.segments:
            self.prefix_lengths.append(total)
            total += segment.end - segment.start
        self.cds_length = total

    def pos_to_cds_coord(self, pos: int) -> int | None:
        """Map a genomic position to 0-based CDS coordinate, if inside CDS."""

        for prefix, segment in zip(self.prefix_lengths, self.segments):
            if segment.start <= pos < segment.end:
                if self.strand == "+":
                    return prefix + (pos - segment.start)
                return prefix + (segment.end - 1 - pos)
        return None

    def landmark_for_coord(self, coord: int) -> str:
        """Return a training landmark name for CDS start or stop-codon start."""

        if coord == 0:
            return "cds_start"
        if self.cds_length >= 3 and coord == self.cds_length - 3:
            return "cds_stop_last_codon"
        return ""


@dataclasses.dataclass(frozen=True)
class CdsAnnotation:
    """Annotation returned for a single genomic position.

    A position may not overlap any CDS; in that case the default values are used
    and frame/distances are set to MISSING_INT.
    """

    transcript_id: str = ""
    gene_id: str = ""
    gene_name: str = ""
    frame: int = MISSING_INT
    distance_to_cds_start: int = MISSING_INT
    distance_to_cds_stop: int = MISSING_INT
    cds_percentile: float = -1.0
    landmark: str = ""


class CdsIndex:
    """Binned lookup from genomic coordinate to overlapping CDS annotations."""

    def __init__(self, models: list[TranscriptModel], bin_size: int = 1000):
        self.models = models
        self.bin_size = bin_size
        self.bins: dict[str, dict[int, list[int]]] = defaultdict(lambda: defaultdict(list))
        for model_idx, model in enumerate(models):
            for segment in model.segments:
                # Store transcript IDs in every genomic bin touched by any CDS
                # segment. Querying a position then checks only nearby models.
                first_bin = segment.start // bin_size
                last_bin = (segment.end - 1) // bin_size
                for bin_id in range(first_bin, last_bin + 1):
                    self.bins[model.contig][bin_id].append(model_idx)

    def annotations_at(self, contig: str, pos: int) -> list[CdsAnnotation]:
        """Return all CDS annotations that overlap contig:pos."""

        model_ids = self.bins.get(contig, {}).get(pos // self.bin_size, [])
        if not model_ids:
            return []

        annotations: list[CdsAnnotation] = []
        seen: set[int] = set()
        for model_idx in model_ids:
            # A transcript can be inserted more than once if multiple CDS
            # segments share this bin. Deduplicate before testing the position.
            if model_idx in seen:
                continue
            seen.add(model_idx)
            model = self.models[model_idx]
            coord = model.pos_to_cds_coord(pos)
            if coord is None:
                continue
            annotations.append(
                CdsAnnotation(
                    transcript_id=model.transcript_id,
                    gene_id=model.gene_id,
                    gene_name=model.gene_name,
                    frame=coord % 3,
                    distance_to_cds_start=coord,
                    distance_to_cds_stop=model.cds_length - 1 - coord,
                    cds_percentile=float(coord / max(model.cds_length - 1, 1)),
                    landmark=model.landmark_for_coord(coord),
                )
            )
        return annotations

    def best_annotation(self, contig: str, pos: int) -> CdsAnnotation:
        """Pick one deterministic annotation when transcripts overlap."""

        annotations = self.annotations_at(contig, pos)
        if not annotations:
            return CdsAnnotation()

        def rank(annotation: CdsAnnotation) -> tuple[int, int, int, str]:
            # Prefer landmark hits because they create labels; otherwise prefer
            # frame 0, then positions near CDS boundaries, then transcript ID
            # for deterministic tie-breaking.
            landmark_rank = 0 if annotation.landmark else 1
            frame_rank = 0 if annotation.frame == 0 else 1
            boundary_distance = min(
                annotation.distance_to_cds_start,
                annotation.distance_to_cds_stop,
            )
            return (landmark_rank, frame_rank, boundary_distance, annotation.transcript_id)

        return min(annotations, key=rank)


class PositionCountIndex:
    """Fast range counter for sparse genomic positions.

    Used for 5-prime and 3-prime read-end densities. Instead of scanning every
    observed read end for every candidate P-site, this stores sorted end
    positions and prefix sums per contig.
    """

    def __init__(self, counter_by_contig: dict[str, Counter[int]]):
        self.positions: dict[str, np.ndarray] = {}
        self.prefix: dict[str, np.ndarray] = {}
        for contig, counter in counter_by_contig.items():
            if not counter:
                continue
            # Example:
            #   counter: {100: 3, 105: 1, 120: 7}
            #   positions: [100, 105, 120]
            #   prefix: [0, 3, 4, 11]
            # Then query(100, 121) is prefix[3] - prefix[0] == 11.
            positions = np.fromiter(sorted(counter), dtype=np.int64)
            counts = np.fromiter((counter[int(pos)] for pos in positions), dtype=np.int64)
            self.positions[contig] = positions
            self.prefix[contig] = np.concatenate(([0], np.cumsum(counts)))

    def query(self, contig: str, start: int, end: int) -> int:
        """Count observed positions in the half-open interval [start, end)."""

        positions = self.positions.get(contig)
        if positions is None:
            return 0
        prefix = self.prefix[contig]
        lo = bisect.bisect_left(positions, start)
        hi = bisect.bisect_left(positions, end)
        return int(prefix[hi] - prefix[lo])


class CoverageIndex:
    """Fast mean-coverage lookup for read-alignment intervals.

    The input is a difference array per contig: +1 at each read start and -1 at
    each read end. This class converts those sparse change points into coverage
    integrals, so mean coverage over any window is an O(log n) query.
    """

    def __init__(self, diff_by_contig: dict[str, Counter[int]]):
        self.breaks: dict[str, np.ndarray] = {}
        self.integrals: dict[str, np.ndarray] = {}
        self.coverages: dict[str, np.ndarray] = {}
        for contig, diff in diff_by_contig.items():
            if not diff:
                continue
            breaks: list[int] = []
            integrals: list[int] = []
            coverages: list[int] = []
            total = 0
            current = 0
            previous = 0
            for pos in sorted(diff):
                # Between two breakpoints, coverage is constant. Accumulate the
                # area under the coverage curve before applying the change at
                # the new breakpoint.
                if pos > previous:
                    total += current * (pos - previous)
                current += diff[pos]
                breaks.append(pos)
                integrals.append(total)
                coverages.append(current)
                previous = pos
            self.breaks[contig] = np.asarray(breaks, dtype=np.int64)
            self.integrals[contig] = np.asarray(integrals, dtype=np.int64)
            self.coverages[contig] = np.asarray(coverages, dtype=np.int64)

    def _integral_at(self, contig: str, pos: int) -> int:
        """Return integral coverage from coordinate 0 up to coordinate pos."""

        breaks = self.breaks.get(contig)
        if breaks is None or len(breaks) == 0:
            return 0
        idx = bisect.bisect_right(breaks, pos) - 1
        if idx < 0:
            return 0
        return int(
            self.integrals[contig][idx]
            + self.coverages[contig][idx] * (pos - self.breaks[contig][idx])
        )

    def mean_coverage(self, contig: str, start: int, end: int) -> float:
        """Average read coverage over [start, end), clamped away from negatives."""

        start = max(0, start)
        end = max(start + 1, end)
        covered = self._integral_at(contig, end) - self._integral_at(contig, start)
        return float(covered / (end - start))


def parse_attributes(raw: str) -> dict[str, str]:
    """Parse the quoted key-value attribute field used by these GTF files."""

    return {match.group(1): match.group(2) for match in re.finditer(r'(\S+) "([^"]*)"', raw)}


def parse_cds_index(gtf_path: pathlib.Path) -> CdsIndex:
    """Load CDS features from a GTF and build the coordinate lookup index."""

    transcript_rows: dict[str, dict[str, object]] = {}
    with gtf_path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            attrs = parse_attributes(fields[8])
            transcript_id = attrs.get("transcript_id")
            if not transcript_id:
                continue
            row = transcript_rows.setdefault(
                transcript_id,
                {
                    "transcript_id": transcript_id,
                    "gene_id": attrs.get("gene_id", ""),
                    "gene_name": attrs.get("gene_name", ""),
                    "contig": fields[0],
                    "strand": fields[6],
                    "segments": [],
                },
            )
            # GTF coordinates are 1-based inclusive. Convert to the 0-based
            # half-open convention used by pysam and the rest of this module.
            row["segments"].append(CdsSegment(int(fields[3]) - 1, int(fields[4])))

    # Build one transcript model per transcript_id. The model constructor sorts
    # each transcript's segments in transcript orientation.
    models = [
        TranscriptModel(
            transcript_id=str(row["transcript_id"]),
            gene_id=str(row["gene_id"]),
            gene_name=str(row["gene_name"]),
            contig=str(row["contig"]),
            strand=str(row["strand"]),
            segments=list(row["segments"]),
        )
        for row in transcript_rows.values()
        if row["segments"]
    ]
    return CdsIndex(models)


def read_sample_table(path: pathlib.Path) -> list[dict[str, str]]:
    """Read the tab-delimited run table that drives batch extraction."""

    with path.open(newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def load_runinfo(metadata_dir: pathlib.Path) -> dict[str, dict[str, str]]:
    """Load optional SRA RunInfo CSV metadata keyed by Run ID."""

    out: dict[str, dict[str, str]] = {}
    for path in metadata_dir.glob("*.runinfo.csv"):
        if path.name.startswith("."):
            continue
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                run = row.get("Run", "")
                if run.startswith("SRR"):
                    out[run] = row
    return out


def infer_nuclease(study_id: str, library_name: str, library_selection: str) -> str:
    """Best-effort protocol/nuclease label for known p-site offset prediction studies."""

    lower_name = library_name.lower()
    lower_selection = library_selection.lower()
    if "mrna" in lower_name or "polya" in lower_selection:
        return "rna_seq_control"
    if study_id in {"SRP033366", "PRJNA842067"}:
        return "rnase_i"
    if study_id == "SRP067542":
        return "mnase"
    if study_id == "SRP419250":
        if "ottr" in lower_name:
            return "ottr"
        if lower_name.startswith("yp") or "p1" in lower_name:
            return "p1"
        if lower_name.startswith("yr") or "rnase" in lower_name:
            return "rnase_i"
    return "unknown"


def sample_contexts(
    sample_rows: Iterable[dict[str, str]],
    runinfo: dict[str, dict[str, str]],
) -> list[SampleContext]:
    """Merge sample-table rows with optional RunInfo metadata."""

    contexts: list[SampleContext] = []
    for row in sample_rows:
        sample_id = row["run"]
        info = runinfo.get(sample_id, {})
        study_id = row.get("study_key") or info.get("SRAStudy") or ""
        library_name = info.get("LibraryName") or info.get("SampleName") or row.get("library_name") or ""
        library_method = (
            info.get("LibrarySelection")
            or info.get("LibraryStrategy")
            or row.get("library_method")
            or "unknown"
        )
        # Prefer explicit sample-table fields, then fall back to study-specific
        # inference from SRA metadata.
        nuclease = row.get("nuclease") or infer_nuclease(
            study_id,
            library_name,
            library_method,
        )
        protocol = row.get("protocol") or nuclease
        contexts.append(
            SampleContext(
                sample_id=sample_id,
                organism=row.get("organism") or info.get("ScientificName") or "",
                ref_key=row["ref_key"],
                study_id=study_id,
                protocol=protocol,
                nuclease=nuclease,
                library_method=library_method,
                library_name=library_name,
            )
        )
    return contexts


def parse_offsets(raw: str) -> tuple[int, ...]:
    """Parse either a range such as '8-20' or a comma list such as '12,15,18'."""

    if "-" in raw:
        start, end = raw.split("-", 1)
        return tuple(range(int(start), int(end) + 1))
    return tuple(int(item) for item in raw.split(",") if item)


def ensure_fasta_index(fasta_path: pathlib.Path) -> None:
    """Create genome.fa.fai when pysam needs random access and no index exists."""

    if not pathlib.Path(str(fasta_path) + ".fai").exists():
        pysam.faidx(str(fasta_path))


def oriented_sequence(
    fasta: pysam.FastaFile,
    contig: str,
    start: int,
    end: int,
    strand: str,
) -> str:
    """Fetch sequence in read/transcript orientation, padding off-contig bases."""

    length = max(0, end - start)
    contig_length = fasta.get_reference_length(contig)
    fetch_start = max(0, start)
    fetch_end = min(contig_length, end)
    left_pad = "N" * max(0, -start)
    right_pad = "N" * max(0, end - contig_length)
    seq = left_pad + fasta.fetch(contig, fetch_start, fetch_end).upper() + right_pad
    if len(seq) < length:
        seq += "N" * (length - len(seq))
    # Feature strings are oriented relative to the read/candidate strand. A
    # minus-strand candidate therefore sees reverse-complemented sequence.
    if strand == "-":
        seq = str(Seq(seq).reverse_complement())
    return seq


def centered_kmer(
    fasta: pysam.FastaFile,
    contig: str,
    pos: int,
    strand: str,
    radius: int,
) -> str:
    """Return a strand-oriented k-mer centered on pos."""

    return oriented_sequence(fasta, contig, pos - radius, pos + radius + 1, strand)


def codon_from_first_base(
    fasta: pysam.FastaFile,
    contig: str,
    first_base_pos: int,
    strand: str,
) -> str:
    """Return the strand-oriented codon whose first base is first_base_pos."""

    if strand == "+":
        return oriented_sequence(fasta, contig, first_base_pos, first_base_pos + 3, strand)
    return oriented_sequence(fasta, contig, first_base_pos - 2, first_base_pos + 1, strand)


def five_three_prime_positions(read: pysam.AlignedSegment) -> tuple[int, int, str]:
    """Return strand-aware 5-prime and 3-prime genomic positions for a read."""

    strand = "-" if read.is_reverse else "+"
    if strand == "+":
        return int(read.reference_start), int(read.reference_end - 1), strand
    return int(read.reference_end - 1), int(read.reference_start), strand


def reference_pos_for_candidate(read: pysam.AlignedSegment, offset: int) -> int:
    """Map a read-oriented candidate offset to a genomic reference position."""

    query_len = read.query_length or read.infer_query_length() or read.query_alignment_length
    if not query_len:
        five_prime_pos, _, strand = five_three_prime_positions(read)
        return five_prime_pos + offset if strand == "+" else five_prime_pos - offset

    # Use the read CIGAR to map query coordinates to reference coordinates. This
    # preserves candidate positions across indels instead of assuming a perfect
    # ungapped alignment.
    query_to_ref = {
        query_pos: ref_pos
        for query_pos, ref_pos in read.get_aligned_pairs(matches_only=False)
        if query_pos is not None
    }
    if read.is_reverse:
        qpos = query_len - 1 - offset
        direction = -1
    else:
        qpos = offset
        direction = 1
    ref_pos = query_to_ref.get(qpos)
    if ref_pos is not None:
        return int(ref_pos)

    # If the exact query base is an insertion/soft-clipped base with no reference
    # coordinate, walk outward to the nearest aligned query base.
    for step in range(1, query_len + 1):
        next_qpos = qpos + direction * step
        if 0 <= next_qpos < query_len and query_to_ref.get(next_qpos) is not None:
            return int(query_to_ref[next_qpos])
        previous_qpos = qpos - direction * step
        if 0 <= previous_qpos < query_len and query_to_ref.get(previous_qpos) is not None:
            return int(query_to_ref[previous_qpos])

    # Extremely defensive fallback for unusual reads with no usable aligned
    # pairs after pysam inference.
    five_prime_pos, _, strand = five_three_prime_positions(read)
    return five_prime_pos + offset if strand == "+" else five_prime_pos - offset


def is_primary_mapped(read: pysam.AlignedSegment, min_mapq: int) -> bool:
    """Filter to mapped primary alignments at or above the MAPQ threshold."""

    return (
        not read.is_unmapped
        and not read.is_secondary
        and not read.is_supplementary
        and read.mapping_quality >= min_mapq
    )


def build_density_indexes(
    bam_path: pathlib.Path,
    min_mapq: int,
) -> tuple[CoverageIndex, PositionCountIndex, PositionCountIndex, int]:
    """Scan a BAM once to build local-density feature indexes."""

    diff_by_contig: dict[str, Counter[int]] = defaultdict(Counter)
    five_prime_by_contig: dict[str, Counter[int]] = defaultdict(Counter)
    three_prime_by_contig: dict[str, Counter[int]] = defaultdict(Counter)
    usable_reads = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if not is_primary_mapped(read, min_mapq):
                continue
            usable_reads += 1
            contig = read.reference_name
            five_prime_pos, three_prime_pos, _strand = five_three_prime_positions(read)
            # End-density features count strand-aware read ends near a candidate
            # site. Coverage features count every base overlapped by a read.
            five_prime_by_contig[contig][five_prime_pos] += 1
            three_prime_by_contig[contig][three_prime_pos] += 1
            diff_by_contig[contig][int(read.reference_start)] += 1
            diff_by_contig[contig][int(read.reference_end)] -= 1

    return (
        CoverageIndex(diff_by_contig),
        PositionCountIndex(five_prime_by_contig),
        PositionCountIndex(three_prime_by_contig),
        usable_reads,
    )


class ParquetSink:
    """Buffered Parquet writer with fixed schema and chunked flushing."""

    def __init__(self, path: pathlib.Path, schema: pa.Schema, columns: list[str]):
        self.path = path
        self.schema = schema
        self.columns = columns
        self.writer: pq.ParquetWriter | None = None
        self.rows: list[dict[str, object]] = []
        self.row_count = 0
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def append(self, rows: Iterable[dict[str, object]]) -> None:
        """Append rows to the in-memory buffer and flush large batches."""

        for row in rows:
            self.rows.append(row)
        if len(self.rows) >= 250_000:
            self.flush()

    def flush(self) -> None:
        """Write buffered rows to the Parquet stream."""

        if not self.rows:
            return
        frame = pd.DataFrame.from_records(self.rows, columns=self.columns)
        table = pa.Table.from_pandas(frame, schema=self.schema, preserve_index=False)
        if self.writer is None:
            self.writer = pq.ParquetWriter(
                self.path,
                self.schema,
                compression="zstd",
                use_dictionary=True,
            )
        self.writer.write_table(table)
        self.row_count += len(self.rows)
        self.rows = []

    def close(self) -> None:
        """Flush and close the writer; remove empty files left by no-row sinks."""

        self.flush()
        if self.writer is not None:
            self.writer.close()
        elif self.path.exists():
            self.path.unlink()


def candidate_rows_for_read(
    read: pysam.AlignedSegment,
    context: SampleContext,
    fasta: pysam.FastaFile,
    cds_index: CdsIndex,
    coverage_index: CoverageIndex,
    five_prime_index: PositionCountIndex,
    three_prime_index: PositionCountIndex,
    offsets: tuple[int, ...],
    kmer_radius: int,
    local_window: int,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    """Create all candidate-offset rows plus the compact read row for one read."""

    contig = read.reference_name
    five_prime_pos, three_prime_pos, strand = five_three_prime_positions(read)
    read_length = read.query_length or read.infer_query_length() or read.query_alignment_length
    read_id = f"{context.sample_id}|{read.query_name}"

    # The read-end frames are features for downstream models. They are not the
    # labels; labels are assigned only from candidate offsets that hit landmarks.
    five_prime_annotation = cds_index.best_annotation(contig, five_prime_pos)
    three_prime_annotation = cds_index.best_annotation(contig, three_prime_pos)

    # Shared read-level fields copied into every candidate row for this read.
    row_base = {
        "read_id": read_id,
        "sample_id": context.sample_id,
        "study_id": context.study_id,
        "organism": context.organism,
        "protocol": context.protocol,
        "nuclease": context.nuclease,
        "library_method": context.library_method,
        "library_name": context.library_name,
        "read_length": int(read_length or 0),
        "strand": strand,
        "contig": contig,
        "mapq": int(read.mapping_quality),
        "five_prime_pos": int(five_prime_pos),
        "three_prime_pos": int(three_prime_pos),
        "five_prime_frame": int(five_prime_annotation.frame),
        "three_prime_frame": int(three_prime_annotation.frame),
    }

    rows: list[dict[str, object]] = []
    positive_offsets: dict[int, str] = {}
    for offset in offsets:
        # p_pos is the candidate P-site coordinate for this offset. A-site and
        # E-site codons are one codon downstream/upstream in transcript
        # orientation, which flips sign on the minus strand.
        p_pos = reference_pos_for_candidate(read, offset)
        a_pos = p_pos + 3 if strand == "+" else p_pos - 3
        e_pos = p_pos - 3 if strand == "+" else p_pos + 3
        annotation = cds_index.best_annotation(contig, p_pos)

        # A candidate is "positive" only when it lands exactly on a known CDS
        # landmark and the landmark is in coding frame 0.
        if annotation.landmark and annotation.frame == 0:
            positive_offsets[offset] = annotation.landmark
        rows.append(
            {
                **row_base,
                "candidate_offset": int(offset),
                "candidate_ref_pos": int(p_pos),
                "candidate_frame": int(annotation.frame),
                "distance_to_cds_start": int(annotation.distance_to_cds_start),
                "distance_to_cds_stop": int(annotation.distance_to_cds_stop),
                "cds_percentile": float(annotation.cds_percentile),
                "transcript_id": annotation.transcript_id,
                "gene_id": annotation.gene_id,
                "gene_name": annotation.gene_name,
                "landmark": annotation.landmark,
                # Sequence features use strand-oriented sequence, so the same
                # column is comparable across plus- and minus-strand reads.
                "five_prime_kmer": centered_kmer(fasta, contig, five_prime_pos, strand, kmer_radius),
                "three_prime_kmer": centered_kmer(
                    fasta,
                    contig,
                    three_prime_pos,
                    strand,
                    kmer_radius,
                ),
                "candidate_site_kmer": centered_kmer(fasta, contig, p_pos, strand, kmer_radius),
                "candidate_p_codon": codon_from_first_base(fasta, contig, p_pos, strand),
                "candidate_a_codon": codon_from_first_base(fasta, contig, a_pos, strand),
                "candidate_e_codon": codon_from_first_base(fasta, contig, e_pos, strand),
                # Local density features summarize the pileup around the
                # candidate site in the same window for every offset.
                "local_coverage_pm15": float(
                    coverage_index.mean_coverage(
                        contig,
                        p_pos - local_window,
                        p_pos + local_window + 1,
                    )
                ),
                "local_5p_end_density": int(
                    five_prime_index.query(
                        contig,
                        p_pos - local_window,
                        p_pos + local_window + 1,
                    )
                ),
                "local_3p_end_density": int(
                    three_prime_index.query(
                        contig,
                        p_pos - local_window,
                        p_pos + local_window + 1,
                    )
                ),
                "label": None,
                "label_tier": "",
            }
        )

    if len(positive_offsets) == 1:
        # Tier 1 silver label: exactly one candidate offset explains a CDS
        # start/stop landmark. The winning offset is positive and all other
        # offsets for this read are negatives.
        silver_offset, silver_landmark = next(iter(positive_offsets.items()))
        for row in rows:
            row["label"] = 1 if row["candidate_offset"] == silver_offset else 0
            row["label_tier"] = "start_stop_landmark"
    elif len(positive_offsets) > 1:
        # Multiple landmark-compatible offsets make the read unusable as a
        # single-positive training example, but we still count it in the summary.
        silver_offset = MISSING_INT
        silver_landmark = "ambiguous_landmark"
        for row in rows:
            row["label_tier"] = "ambiguous_landmark"
    else:
        # Most reads do not cross a CDS start/stop landmark at any candidate
        # offset. They stay unlabeled in the read table.
        silver_offset = MISSING_INT
        silver_landmark = ""

    # The read table is useful for summaries and grouping even when candidate
    # rows are filtered by output mode.
    read_row = {
        **row_base,
        "silver_offset": int(silver_offset),
        "silver_landmark": silver_landmark,
        "label_tier": "start_stop_landmark" if len(positive_offsets) == 1 else silver_landmark,
    }
    return rows, read_row


def extract_sample(
    context: SampleContext,
    ref_dir: pathlib.Path,
    bam_path: pathlib.Path,
    outdir: pathlib.Path,
    offsets: tuple[int, ...],
    output_mode: str,
    min_mapq: int,
    max_reads: int | None,
    kmer_radius: int,
    local_window: int,
) -> dict[str, object]:
    """Extract all p-site offset prediction features for one sample/BAM."""

    gtf_path = ref_dir / "annotation.gtf"
    fasta_path = ref_dir / "genome.fa"
    ensure_fasta_index(fasta_path)

    # Reference annotation and sequence are sample-specific through ref_key
    # because the current p-site offset prediction set mixes Drosophila and yeast runs.
    cds_index = parse_cds_index(gtf_path)
    fasta = pysam.FastaFile(str(fasta_path))

    # First BAM pass: build coverage and read-end density indexes. The second
    # pass can then query these indexes cheaply for every candidate offset.
    coverage_index, five_prime_index, three_prime_index, density_reads = build_density_indexes(
        bam_path,
        min_mapq,
    )

    table_dir = outdir / "tables"
    candidate_path = (
        table_dir / "all_candidates_by_sample" / f"{context.sample_id}.candidate_offsets.parquet"
    )
    landmark_path = (
        table_dir / "candidates_by_sample" / f"{context.sample_id}.candidate_offsets.parquet"
    )
    read_path = table_dir / "reads_by_sample" / f"{context.sample_id}.reads.parquet"

    # output_mode controls only candidate-row sinks. The read table is always
    # written because it is compact and useful for diagnostics.
    candidate_sink = (
        ParquetSink(candidate_path, CANDIDATE_SCHEMA, CANDIDATE_COLUMNS)
        if output_mode in {"all", "both"}
        else None
    )
    landmark_sink = (
        ParquetSink(landmark_path, CANDIDATE_SCHEMA, CANDIDATE_COLUMNS)
        if output_mode in {"labeled", "both"}
        else None
    )
    read_sink = ParquetSink(read_path, READ_SCHEMA, READ_COLUMNS)

    processed_reads = 0
    labeled_reads = 0
    ambiguous_reads = 0
    candidate_rows = 0
    landmark_rows = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if not is_primary_mapped(read, min_mapq):
                continue
            if max_reads is not None and processed_reads >= max_reads:
                break

            # Second BAM pass: emit one read row and one candidate row per
            # requested offset. With default offsets 8-20, that is 13 candidate
            # rows per processed read before output-mode filtering.
            rows, read_row = candidate_rows_for_read(
                read,
                context,
                fasta,
                cds_index,
                coverage_index,
                five_prime_index,
                three_prime_index,
                offsets,
                kmer_radius,
                local_window,
            )
            processed_reads += 1
            candidate_rows += len(rows)
            if read_row["label_tier"] == "start_stop_landmark":
                labeled_reads += 1
                landmark_rows += len(rows)
                if landmark_sink is not None:
                    # In "labeled" mode, only these labeled landmark candidate
                    # rows are written to the candidate-feature table.
                    landmark_sink.append(rows)
            elif read_row["label_tier"] == "ambiguous_landmark":
                ambiguous_reads += 1

            if candidate_sink is not None:
                # In "all" or "both" mode, every candidate row for every
                # processed read is written.
                candidate_sink.append(rows)
            read_sink.append([read_row])

    for sink in (candidate_sink, landmark_sink, read_sink):
        if sink is not None:
            sink.close()
    fasta.close()

    # This summary row becomes tables/feature_summary.tsv.
    return {
        "sample_id": context.sample_id,
        "study_id": context.study_id,
        "organism": context.organism,
        "ref_key": context.ref_key,
        "bam_path": str(bam_path),
        "density_reads": density_reads,
        "processed_reads": processed_reads,
        "candidate_rows": candidate_rows,
        "labeled_reads": labeled_reads,
        "landmark_candidate_rows": landmark_rows,
        "ambiguous_landmark_reads": ambiguous_reads,
        "output_mode": output_mode,
        "candidate_path": str(candidate_path) if candidate_sink is not None else "",
        "landmark_path": str(landmark_path) if landmark_sink is not None else "",
        "read_path": str(read_path),
    }


def combine_parquet_files(paths: list[pathlib.Path], output_path: pathlib.Path) -> str:
    """Concatenate existing per-sample Parquet files into one combined file."""

    existing_paths = [path for path in paths if path.exists()]
    if not existing_paths:
        return ""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    writer: pq.ParquetWriter | None = None
    try:
        for path in existing_paths:
            table = pq.read_table(path)
            if writer is None:
                # Use the first table's schema as the combined output schema.
                writer = pq.ParquetWriter(
                    output_path,
                    table.schema,
                    compression="zstd",
                    use_dictionary=True,
                )
            writer.write_table(table)
    finally:
        if writer is not None:
            writer.close()
    return str(output_path)


def main(argv: list[str] | None = None) -> int:
    """CLI entry point for batch p-site offset prediction feature extraction."""

    parser = argparse.ArgumentParser(
        description="Extract p-site offset candidate features from aligned BAMs."
    )
    parser.add_argument("--sample-table", type=pathlib.Path, required=True)
    parser.add_argument("--metadata-dir", type=pathlib.Path, default=pathlib.Path("data/metadata"))
    parser.add_argument("--refs-dir", type=pathlib.Path, default=pathlib.Path("data/references"))
    parser.add_argument("--align-dir", type=pathlib.Path, default=pathlib.Path("data/intermediates/alignments/genome_bams"))
    parser.add_argument(
        "--outdir",
        type=pathlib.Path,
        default=pathlib.Path("data/datasets/p_site_offset_prediction/landmark_riboseq"),
    )
    parser.add_argument("--offsets", default="8-20")
    parser.add_argument(
        "--output-mode",
        choices=("labeled", "all", "both"),
        default="labeled",
        help="labeled writes only Tier 1 landmark candidate rows; all writes every candidate row.",
    )
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument("--max-reads-per-sample", type=int)
    parser.add_argument("--kmer-radius", type=int, default=3)
    parser.add_argument("--local-window", type=int, default=15)
    parser.add_argument(
        "--skip-combined",
        action="store_true",
        help="Do not write combined cross-sample Parquet files.",
    )
    args = parser.parse_args(argv)

    # Build one SampleContext per table row before touching BAMs, so metadata
    # fallback/inference is handled consistently across samples.
    sample_rows = read_sample_table(args.sample_table)
    contexts = sample_contexts(sample_rows, load_runinfo(args.metadata_dir))
    offsets = parse_offsets(args.offsets)
    args.outdir.mkdir(parents=True, exist_ok=True)

    summaries: list[dict[str, object]] = []
    for context in contexts:
        # The aligner stage writes one coordinate-sorted BAM per SRA run using
        # this fixed naming convention.
        bam_path = args.align_dir / f"{context.sample_id}.Aligned.sortedByCoord.out.bam"
        if not bam_path.exists():
            print(f"missing BAM for {context.sample_id}: {bam_path}")
            continue
        ref_dir = args.refs_dir / context.ref_key
        print(f"extracting {context.sample_id} with {context.ref_key} annotations")
        summaries.append(
            extract_sample(
                context=context,
                ref_dir=ref_dir,
                bam_path=bam_path,
                outdir=args.outdir,
                offsets=offsets,
                output_mode=args.output_mode,
                min_mapq=args.min_mapq,
                max_reads=args.max_reads_per_sample,
                kmer_radius=args.kmer_radius,
                local_window=args.local_window,
            )
        )

    if summaries:
        # The summary is small and human-readable; it records row counts and
        # paths for every sample that was actually processed.
        summary_path = args.outdir / "tables" / "feature_summary.tsv"
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(summaries).to_csv(summary_path, sep="\t", index=False)
        print(f"wrote {summary_path}")
        if not args.skip_combined:
            # Combined outputs are convenient for downstream splitting/training.
            # Empty path fields are skipped, which is why labeled mode combines
            # landmarks and reads but not all candidates.
            landmark_paths = [
                pathlib.Path(str(row["landmark_path"]))
                for row in summaries
                if row.get("landmark_path")
            ]
            read_paths = [
                pathlib.Path(str(row["read_path"]))
                for row in summaries
                if row.get("read_path")
            ]
            candidate_paths = [
                pathlib.Path(str(row["candidate_path"]))
                for row in summaries
                if row.get("candidate_path")
            ]
            combined_landmarks = combine_parquet_files(
                landmark_paths,
                args.outdir / "tables" / "candidate_offsets.parquet",
            )
            combined_reads = combine_parquet_files(
                read_paths,
                args.outdir / "tables" / "reads.parquet",
            )
            combined_candidates = combine_parquet_files(
                candidate_paths,
                args.outdir / "tables" / "all_candidate_offsets.parquet",
            )
            for path in (combined_landmarks, combined_reads, combined_candidates):
                if path:
                    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
