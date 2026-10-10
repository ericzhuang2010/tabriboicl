from __future__ import annotations

FEATURE_SETS: dict[str, list[str]] = {
    "full": [
        "candidate_offset",
        "read_length",
        "five_prime_frame",
        "three_prime_frame",
        "distance_to_cds_start",
        "distance_to_cds_stop",
        "cds_percentile",
        "five_prime_kmer",
        "three_prime_kmer",
        "candidate_site_kmer",
        "candidate_p_codon",
        "candidate_a_codon",
        "candidate_e_codon",
        "local_coverage_pm15",
        "local_5p_end_density",
        "local_3p_end_density",
        "organism",
        "nuclease",
        "library_method",
        "study_id",
    ],
    "no_boundary": [
        "candidate_offset",
        "read_length",
        "strand",
        "mapq",
        "five_prime_kmer",
        "three_prime_kmer",
        "local_coverage_pm15",
        "local_5p_end_density",
        "local_3p_end_density",
        "organism",
        "nuclease",
        "library_method",
        "study_id",
    ],
    "minimal": [
        "candidate_offset",
        "read_length",
        "organism",
        "nuclease",
        "library_method",
        "study_id",
    ],
}


def select_feature_columns(
    train_columns: set[str],
    test_columns: set[str],
    feature_set: str,
) -> list[str]:
    if feature_set not in FEATURE_SETS:
        valid = ", ".join(sorted(FEATURE_SETS))
        raise ValueError(f"unknown feature set '{feature_set}'. Valid feature sets: {valid}")
    return [
        column
        for column in FEATURE_SETS[feature_set]
        if column in train_columns and column in test_columns
    ]
