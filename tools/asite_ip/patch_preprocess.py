#!/usr/bin/env python3
"Apply TabRiboICL's bounded-memory patch to converted Asite-IP preprocessing."

from __future__ import annotations

import argparse
from pathlib import Path


FUNCTION_START = "def samparser_transcriptome("
FUNCTION_END = "\ndef create_cds_counts_transcriptome("
CALL_ORIGINAL = "samparser_transcriptome(SAM, min_frag, max_frag)"
CALL_PATCHED = "samparser_transcriptome(SAM, min_frag, max_frag, idx_file=annotations)"


BOUNDED_MEMORY_FUNCTION = r"""def samparser_transcriptome(
        sfile, frag_min, frag_max, idx_file=None, three_prime=False):
    # Upstream allocates a frag-range list for every observed position before
    # filtering sparse genes. Use the same filter first, then retain compact
    # counts only for the CDS window used by create_cds_counts_transcriptome.
    if idx_file is None:
        raise ValueError("idx_file is required for bounded-memory transcriptome parsing")

    annotation = {}
    with open(idx_file) as annotation_file:
        for line in annotation_file:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 3:
                continue
            annotation[fields[0]] = (int(fields[1]), int(fields[2]))

    total_count = {}
    frag_range = frag_max - frag_min + 1
    read_count = 0
    mul_read_count = 0
    discarded_reads = 0
    progress_every = 1000000

    if three_prime:
        print("Mapping reads by 3' end...")
    else:
        print("Mapping reads by 5' end...")
    print("Pass 1/2: identifying genes that pass the upstream coverage filter")

    counter = 0
    with open(sfile) as sam_file:
        for line in sam_file:
            counter += 1
            if counter % progress_every == 0:
                print("Pass 1/2 SAM lines parsed: {0}".format(counter))
            if not line or line[0] == "@":
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 11:
                continue
            try:
                sam_flag = int(fields[1])
            except ValueError:
                continue
            read_length = len(fields[9])
            if not frag_min <= read_length <= frag_max:
                continue

            gene = fields[2]
            multi_align = len(fields) > 19 and fields[19].strip() != "NH:i:1"
            if sam_flag == 0:
                if multi_align:
                    mul_read_count += 1
                else:
                    total_count[gene] = total_count.get(gene, 0) + 1
                    read_count += 1
            elif sam_flag == 16:
                discarded_reads += 1

    active_genes = set()
    bounds = {}
    for gene, count in total_count.items():
        annotation_row = annotation.get(gene)
        if annotation_row is None:
            continue
        start_pos, cds_len = annotation_row
        if start_pos in (0, -1) or count < cds_len:
            continue
        active_genes.add(gene)
        # Conservative endpoints include all positions retained downstream.
        bounds[gene] = (max(1, start_pos - 50), start_pos + cds_len + 51)

    print(
        "Pass 1/2 complete: {0} reads across {1} genes; "
        "{2} genes pass the coverage filter".format(
            read_count, len(total_count), len(active_genes)
        )
    )
    print("Pass 2/2: accumulating compact CDS-window counts")

    dict_count = {gene: {} for gene in active_genes}
    dict_mul_count = {}
    counter = 0
    with open(sfile) as sam_file:
        for line in sam_file:
            counter += 1
            if counter % progress_every == 0:
                print("Pass 2/2 SAM lines parsed: {0}".format(counter))
            if not line or line[0] == "@":
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 11:
                continue

            gene = fields[2]
            if gene not in active_genes:
                continue
            try:
                sam_flag = int(fields[1])
            except ValueError:
                continue
            if sam_flag not in (0, 256):
                continue

            read_length = len(fields[9])
            if not frag_min <= read_length <= frag_max:
                continue
            position = int(fields[3])
            if three_prime:
                position += read_length - 1
            lower_bound, upper_bound = bounds[gene]
            if position < lower_bound or position > upper_bound:
                continue

            position_counts = dict_count[gene].get(position)
            if position_counts is None:
                position_counts = array("I", [0]) * (frag_range + 1)
                dict_count[gene][position] = position_counts

            length_index = read_length - frag_min
            multi_align = len(fields) > 19 and fields[19].strip() != "NH:i:1"
            if sam_flag == 0 and not multi_align:
                position_counts[length_index] += 1
            elif sam_flag == 0 and multi_align:
                position_counts[-1] += 1
                gene_multi = dict_mul_count.setdefault(gene, {})
                multi_counts = gene_multi.get(position)
                if multi_counts is None:
                    multi_counts = array("I", [0]) * (frag_range + 1)
                    gene_multi[position] = multi_counts
                multi_counts[length_index] += 1
            elif sam_flag == 256:
                # Preserve the upstream secondary-alignment marker. The
                # per-length count exists only when NH marks it non-unique.
                position_counts[-1] += 1
                if multi_align:
                    gene_multi = dict_mul_count.setdefault(gene, {})
                    multi_counts = gene_multi.get(position)
                    if multi_counts is None:
                        multi_counts = array("I", [0]) * (frag_range + 1)
                        gene_multi[position] = multi_counts
                    multi_counts[length_index] += 1

    retained_positions = sum(len(values) for values in dict_count.values())
    print(
        "SAM file parsed for total {0} reads mapping onto {1} genes.".format(
            read_count, len(total_count)
        )
    )
    print(
        "{0} reads are multiple aligned; {1} reverse-strand reads were discarded.".format(
            mul_read_count, discarded_reads
        )
    )
    print(
        "Retained {0} genes and {1} CDS-window positions in memory.".format(
            len(active_genes), retained_positions
        )
    )
    return dict_count, dict_mul_count, total_count

"""


def patch_source(text: str) -> str:
    if "from array import array\n" not in text:
        anchor = "import os\n"
        if anchor not in text:
            raise ValueError("expected import anchor not found")
        text = text.replace(anchor, anchor + "from array import array\n", 1)

    start = text.find(FUNCTION_START)
    end = text.find(FUNCTION_END, start)
    if start < 0 or end < 0:
        raise ValueError("expected transcriptome parser function boundaries not found")
    text = text[:start] + BOUNDED_MEMORY_FUNCTION + text[end + 1 :]

    if CALL_ORIGINAL not in text:
        raise ValueError("expected transcriptome parser call not found")
    return text.replace(CALL_ORIGINAL, CALL_PATCHED, 1)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("converted_preprocess", type=Path)
    args = parser.parse_args()
    path = args.converted_preprocess
    path.write_text(patch_source(path.read_text()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
