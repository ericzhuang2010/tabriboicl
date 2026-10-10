#!/usr/bin/env python3
"Apply TabRiboICL's compact-count patch to the converted Asite-IP runtime."

from __future__ import annotations

import argparse
from pathlib import Path


CLASS_ANCHOR = "CODON_TYPES = ["
LOADER_START = "                count_list = list(map(int, line_list[3].split(',')))"
LOADER_END = '\n\n        sys.stdout.write("Parsed read counts for fragment size {0}.\\t\\r".format(fsize))'


COMPACT_CLASS = r'''class CompactReadCounts:
    # Array-backed integer-position mapping used by the unchanged IP algorithm.
    __slots__ = ("start", "values")

    def __init__(self, start, values):
        self.start = start
        self.values = values

    @classmethod
    def from_text(cls, start, gene_len, fsize, raw_counts, three_prime=False):
        values = array("I")
        first_position = None
        position = start
        for raw_value in raw_counts.split(","):
            if three_prime:
                include = 1 <= position <= gene_len + fsize
            else:
                include = -fsize <= position <= gene_len
            if include:
                if first_position is None:
                    first_position = position
                values.append(int(raw_value))
            position += 1
            if position == 0:
                position += 1
        if first_position is None:
            first_position = start
        return cls(first_position, values)

    def __getitem__(self, position):
        if position == 0:
            raise KeyError(position)
        offset = position - self.start
        if self.start < 0 and position > 0:
            offset -= 1
        if offset < 0 or offset >= len(self.values):
            raise KeyError(position)
        return self.values[offset]

    def __iter__(self):
        position = self.start
        for _value in self.values:
            yield position
            position += 1
            if position == 0:
                position += 1

    def __len__(self):
        return len(self.values)

    def __reduce__(self):
        return self.__class__, (self.start, self.values)


'''


COMPACT_LOADER = '''                compact_counts = CompactReadCounts.from_text(
                    start_index,
                    gene_len,
                    fsize,
                    line_list[3],
                    three_prime=three_prime,
                )
                reads_dict[fsize][gene_name] = compact_counts
                if gene_name not in total_reads:
                    total_reads[gene_name] = 0
                total_reads[gene_name] += sum(compact_counts.values)'''


def patch_source(text: str) -> str:
    if "from array import array\n" not in text:
        import_anchor = "import os\n"
        if import_anchor not in text:
            raise ValueError("expected import anchor not found")
        text = text.replace(import_anchor, import_anchor + "from array import array\n", 1)

    if "class CompactReadCounts:" not in text:
        class_index = text.find(CLASS_ANCHOR)
        if class_index < 0:
            raise ValueError("expected class insertion anchor not found")
        text = text[:class_index] + COMPACT_CLASS + text[class_index:]

    loader_start = text.find(LOADER_START)
    loader_end = text.find(LOADER_END, loader_start)
    if loader_start < 0 or loader_end < 0:
        raise ValueError("expected read-count loader block not found")
    text = text[:loader_start] + COMPACT_LOADER + text[loader_end:]

    minimum = "extra_s = min(dict_reads)"
    replacement = (
        "extra_s = (dict_reads.start if hasattr(dict_reads, 'start') "
        "else min(dict_reads))"
    )
    if minimum not in text:
        raise ValueError("expected minimum-position lookup not found")
    return text.replace(minimum, replacement)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("converted_run", type=Path)
    args = parser.parse_args()
    path = args.converted_run
    path.write_text(patch_source(path.read_text()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
