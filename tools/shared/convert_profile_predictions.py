#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import pathlib
from collections import defaultdict
from dataclasses import dataclass

import pandas as pd
import pysam


@dataclass(frozen=True)
class TranscriptInfo:
    utr5_len: int = 0
    cds_len: int = 0
    mix_start: int | None = None


def read_query_length(read: pysam.AlignedSegment) -> int:
    return int(read.query_length or read.infer_query_length() or read.query_alignment_length or 0)


def reference_pos_for_candidate(read: pysam.AlignedSegment, offset: int) -> int | None:
    query_len = read_query_length(read)
    if not query_len:
        return None

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

    for step in range(1, query_len + 1):
        next_qpos = qpos + direction * step
        if 0 <= next_qpos < query_len and query_to_ref.get(next_qpos) is not None:
            return int(query_to_ref[next_qpos])
        previous_qpos = qpos - direction * step
        if 0 <= previous_qpos < query_len and query_to_ref.get(previous_qpos) is not None:
            return int(query_to_ref[previous_qpos])
    return None


def read_sample_table(path: pathlib.Path) -> dict[str, str]:
    samples: dict[str, str] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample = str(row.get("run", "")).strip()
            ref_key = str(row.get("ref_key", "")).strip()
            if sample:
                samples[sample] = ref_key
    return samples


def read_ribodeblur_annotation(annotation_dir: pathlib.Path, ref_key: str) -> dict[str, TranscriptInfo]:
    path = annotation_dir / f"{ref_key}.ribowaltz_annotation.tsv"
    frame = pd.read_csv(path, sep="\t")
    return {
        str(row.transcript): TranscriptInfo(utr5_len=int(row.l_utr5), cds_len=int(row.l_cds))
        for row in frame.itertuples(index=False)
    }


def strip_ribodeblur_suffix(transcript_id: str) -> str:
    suffix = " (no deblur)"
    if transcript_id.endswith(suffix):
        return transcript_id[: -len(suffix)]
    return transcript_id


def iter_profile_values(line: str, start: int):
    index = 0
    cursor = start
    length = len(line)
    while cursor < length:
        while cursor < length and line[cursor].isspace():
            cursor += 1
        if cursor >= length:
            break
        end = cursor + 1
        while end < length and not line[end].isspace():
            end += 1
        yield index, float(line[cursor:end])
        index += 1
        cursor = end


def ribodeblur_profile_index(
    *,
    read: pysam.AlignedSegment,
    offset: int,
    transcript_info: dict[str, TranscriptInfo],
) -> tuple[str, int] | None:
    transcript_id = str(read.reference_name)
    info = transcript_info.get(transcript_id)
    if info is None:
        return None

    local_pos = reference_pos_for_candidate(read, offset)
    if local_pos is None:
        return None
    cds_index = int(local_pos) - info.utr5_len
    if cds_index < 0:
        return None
    if info.cds_len and cds_index >= info.cds_len:
        return None
    return transcript_id, cds_index


def read_ribodeblur_profile_scores(
    path: pathlib.Path,
    requested: dict[str, set[int]],
    window: int,
) -> dict[tuple[str, int], float]:
    if not requested or not path.exists() or path.stat().st_size == 0:
        return {}

    scores: dict[tuple[str, int], float] = {}
    with path.open() as handle:
        for line in handle:
            sep = line.find(": ")
            if sep < 0:
                continue
            transcript_id = strip_ribodeblur_suffix(line[:sep])
            needed_positions = requested.get(transcript_id)
            if not needed_positions:
                continue
            max_needed = max(needed_positions) + window
            transcript_scores = {position: 0.0 for position in needed_positions}
            for value_index, value in iter_profile_values(line, sep + 2):
                if value_index > max_needed:
                    break
                for position in range(value_index - window, value_index + window + 1):
                    if position in transcript_scores and value > transcript_scores[position]:
                        transcript_scores[position] = value
            scores.update(
                ((transcript_id, position), score)
                for position, score in transcript_scores.items()
            )
    return scores


def read_riboba_tx_map(path: pathlib.Path) -> dict[str, TranscriptInfo]:
    if not path.exists() or path.stat().st_size == 0:
        return {}
    frame = pd.read_csv(path, sep="\t")
    return {
        str(row.tx_name): TranscriptInfo(
            utr5_len=int(getattr(row, "utr5_len", 0)),
            cds_len=int(getattr(row, "cds_len", 0)),
            mix_start=int(row.mix_start),
        )
        for row in frame.itertuples(index=False)
    }


def read_riboba_psite_profile(path: pathlib.Path) -> dict[int, float]:
    by_pos: dict[int, float] = defaultdict(float)
    if not path.exists() or path.stat().st_size == 0:
        return {}
    frame = pd.read_csv(path, sep="\t")
    # Older RiboBA exports escaped tabs instead of writing tab separators.
    # Accept those too so cached profiles don't need to be regenerated.
    if len(frame.columns) == 1 and "\\t" in str(frame.columns[0]):
        frame = pd.read_csv(path, sep=r"\\t", engine="python")
    if frame.empty:
        return {}
    required = {"pos", "weight"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    for row in frame.itertuples(index=False):
        pos = int(row.pos)
        weight = float(row.weight)
        by_pos[pos] += weight
    return dict(by_pos)


def primary_mapped(read: pysam.AlignedSegment) -> bool:
    return not read.is_unmapped and not read.is_secondary and not read.is_supplementary


def choose_alignments(
    *,
    bam_path: pathlib.Path,
    sample_id: str,
    needed_read_ids: set[str],
    usable_transcripts: set[str],
) -> dict[str, pysam.AlignedSegment]:
    alignments: dict[str, pysam.AlignedSegment] = {}
    if not bam_path.exists() or bam_path.stat().st_size == 0:
        return alignments

    with pysam.AlignmentFile(bam_path) as bam:
        for read in bam.fetch(until_eof=True):
            if not primary_mapped(read):
                continue
            read_id = f"{sample_id}|{read.query_name}"
            if read_id not in needed_read_ids:
                continue
            current = alignments.get(read_id)
            if current is None:
                alignments[read_id] = read
                continue
            current_usable = current.reference_name in usable_transcripts
            new_usable = read.reference_name in usable_transcripts
            if new_usable and not current_usable:
                alignments[read_id] = read
    return alignments


def window_max(values: list[float], index: int, radius: int) -> float:
    if index < 0 or index >= len(values):
        return 0.0
    start = max(0, index - radius)
    end = min(len(values), index + radius + 1)
    if start >= end:
        return 0.0
    return max(values[start:end])


def ribodeblur_score(
    *,
    read: pysam.AlignedSegment,
    offset: int,
    profiles: dict[str, list[float]],
    transcript_info: dict[str, TranscriptInfo],
    window: int,
) -> float:
    transcript_id = str(read.reference_name)
    profile = profiles.get(transcript_id)
    info = transcript_info.get(transcript_id)
    if profile is None or info is None:
        return 0.0

    local_pos = reference_pos_for_candidate(read, offset)
    if local_pos is None:
        return 0.0
    cds_index = int(local_pos) - info.utr5_len
    if info.cds_len and cds_index >= info.cds_len:
        return 0.0
    return float(window_max(profile, cds_index, window))


def riboba_score(
    *,
    read: pysam.AlignedSegment,
    offset: int,
    by_pos: dict[int, float],
    transcript_info: dict[str, TranscriptInfo],
    window: int,
) -> float:
    info = transcript_info.get(str(read.reference_name))
    if info is None or info.mix_start is None:
        return 0.0
    local_pos = reference_pos_for_candidate(read, offset)
    if local_pos is None:
        return 0.0

    global_pos = info.mix_start + int(local_pos)
    best = 0.0
    for pos in range(global_pos - window, global_pos + window + 1):
        best = max(best, float(by_pos.get(pos, 0.0)))
    return best


def normalize_scores(scores: list[float], pseudocount: float) -> list[float]:
    if not scores:
        return []
    adjusted = [max(float(score), 0.0) + pseudocount for score in scores]
    total = sum(adjusted)
    if total <= 0:
        return [1.0 / len(scores)] * len(scores)
    return [score / total for score in adjusted]


def predictions_for_method(args: argparse.Namespace, test: pd.DataFrame) -> pd.DataFrame:
    samples = read_sample_table(args.sample_table)
    needed = {
        sample_id: set(group["read_id"].astype(str).unique())
        for sample_id, group in test.groupby("sample_id", sort=False)
    }
    candidates_by_read = {
        str(read_id): [int(value) for value in group["candidate_offset"].astype(int).tolist()]
        for read_id, group in test.groupby("read_id", sort=False)
    }

    rows: list[dict[str, object]] = []
    annotation_cache: dict[str, dict[str, TranscriptInfo]] = {}
    riboba_tx_cache: dict[str, dict[str, TranscriptInfo]] = {}
    for sample_id, read_ids in sorted(needed.items()):
        ref_key = samples.get(sample_id, "")
        if args.method == "ribodeblur":
            if not ref_key:
                continue
            transcript_info = annotation_cache.setdefault(
                ref_key,
                read_ribodeblur_annotation(args.annotation_dir, ref_key),
            )
            profile_path = args.baseline_root / sample_id / f"{sample_id}.profile"
            if not profile_path.exists() or profile_path.stat().st_size == 0:
                continue
            bam_path = args.baseline_root / sample_id / f"{sample_id}.cds_transcripts.bam"
            alignments = choose_alignments(
                bam_path=bam_path,
                sample_id=sample_id,
                needed_read_ids=read_ids,
                usable_transcripts=set(transcript_info),
            )
            profile_requests: dict[str, set[int]] = defaultdict(set)
            profile_indexes: dict[tuple[str, int], tuple[str, int]] = {}
            for read_id, read in alignments.items():
                for offset in candidates_by_read.get(read_id, []):
                    profile_index = ribodeblur_profile_index(
                        read=read,
                        offset=offset,
                        transcript_info=transcript_info,
                    )
                    if profile_index is None:
                        continue
                    transcript_id, cds_index = profile_index
                    profile_requests[transcript_id].add(cds_index)
                    profile_indexes[(read_id, offset)] = profile_index

            profile_scores = read_ribodeblur_profile_scores(
                profile_path,
                dict(profile_requests),
                args.window,
            )
            for read_id, read in alignments.items():
                offsets = candidates_by_read.get(read_id, [])
                scores = [
                    profile_scores.get(profile_indexes.get((read_id, offset)), 0.0)
                    for offset in offsets
                ]
                probs = normalize_scores(scores, args.pseudocount)
                rows.extend(
                    {"read_id": read_id, "candidate_offset": offset, "p_offset": prob}
                    for offset, prob in zip(offsets, probs)
                )
        else:
            if not ref_key:
                continue
            ref_root = args.baseline_root / "refs" / ref_key
            tx_map = riboba_tx_cache.setdefault(
                ref_key,
                read_riboba_tx_map(ref_root / "riboba_tx_map.tsv"),
            )
            profile_path = (
                ref_root
                / "downstream"
                / sample_id
                / f"{sample_id}.riboba_psite_profile.tsv"
            )
            by_pos = read_riboba_psite_profile(profile_path)
            if not tx_map or not by_pos:
                continue
            bam_path = (
                ref_root
                / "prepare_bam"
                / "bam"
                / f"{sample_id}_cds_txsorted.bam"
            )
            alignments = choose_alignments(
                bam_path=bam_path,
                sample_id=sample_id,
                needed_read_ids=read_ids,
                usable_transcripts=set(tx_map),
            )
            for read_id, read in alignments.items():
                offsets = candidates_by_read.get(read_id, [])
                scores = [
                    riboba_score(
                        read=read,
                        offset=offset,
                        by_pos=by_pos,
                        transcript_info=tx_map,
                        window=args.window,
                    )
                    for offset in offsets
                ]
                if max(scores, default=0.0) <= 0.0:
                    # No RiboBA call: leave this read missing so strict evaluation
                    # counts the coverage loss rather than inventing an offset.
                    continue
                probs = normalize_scores(scores, args.pseudocount)
                rows.extend(
                    {"read_id": read_id, "candidate_offset": offset, "p_offset": prob}
                    for offset, prob in zip(offsets, probs)
                )

    if not rows:
        return pd.DataFrame(columns=["read_id", "candidate_offset", "p_offset"])
    out = pd.DataFrame.from_records(rows, columns=["read_id", "candidate_offset", "p_offset"])
    return out.drop_duplicates(["read_id", "candidate_offset"])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Convert profile baselines into p-site offset candidate predictions.")
    parser.add_argument("--method", choices=("ribodeblur", "riboba"), required=True)
    parser.add_argument("--test", type=pathlib.Path, required=True)
    parser.add_argument("--sample-table", type=pathlib.Path, required=True)
    parser.add_argument("--baseline-root", type=pathlib.Path, required=True)
    parser.add_argument("--annotation-dir", type=pathlib.Path, default=pathlib.Path("data/datasets/p_site_offset_prediction/landmark_riboseq/views/annotations"))
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--window", type=int, default=1)
    parser.add_argument("--pseudocount", type=float, default=1e-9)
    args = parser.parse_args(argv)

    test = pd.read_parquet(args.test, columns=["read_id", "sample_id", "candidate_offset"])
    predictions = predictions_for_method(args, test)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    predictions.to_parquet(args.output, index=False)
    print(f"wrote {len(predictions)} {args.method} predictions to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
