"""Refresh SRA metadata and build run lists.

Fetching RunInfo contacts NCBI. The regular benchmark scripts use the saved
metadata and don't need to run this step.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import urllib.parse
import urllib.request


RUNINFO_URL = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo"
DEFAULT_METADATA_DIR = pathlib.Path("data/metadata")


def fetch_runinfo(accession: str) -> bytes:
    query = urllib.parse.urlencode({"acc": accession})
    with urllib.request.urlopen(f"{RUNINFO_URL}?{query}", timeout=120) as response:
        return response.read()


def fetch_runinfo_manifest(manifest: pathlib.Path, outdir: pathlib.Path) -> int:
    outdir.mkdir(parents=True, exist_ok=True)
    with manifest.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = row["runinfo_accession"].strip()
            if not accession or accession == "manual":
                continue
            study_key = row["study_key"].strip()
            target = outdir / f"{study_key}.runinfo.csv"
            data = fetch_runinfo(accession)
            target.write_bytes(data)
            print(f"wrote {target}")
    return 0


def extract_runs_from_runinfo(path: pathlib.Path) -> list[str]:
    runs: set[str] = set()
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            run = row.get("Run", "").strip()
            if run.startswith("SRR"):
                runs.add(run)
    return sorted(runs)


def extract_run_lists(runinfo_dir: pathlib.Path, outdir: pathlib.Path) -> int:
    outdir.mkdir(parents=True, exist_ok=True)
    all_runs: set[str] = set()
    for path in sorted(runinfo_dir.glob("*.runinfo.csv")):
        runs = extract_runs_from_runinfo(path)
        target = outdir / path.name.replace(".runinfo.csv", ".runs.txt")
        target.write_text("\n".join(runs) + ("\n" if runs else ""))
        all_runs.update(runs)
        print(f"{target}\t{len(runs)}")

    all_target = outdir / "all_runs.txt"
    all_target.write_text("\n".join(sorted(all_runs)) + ("\n" if all_runs else ""))
    print(f"{all_target}\t{len(all_runs)}")
    return 0


def read_runs(path: pathlib.Path) -> set[str]:
    with path.open(newline="") as handle:
        return {
            row["Run"].strip()
            for row in csv.DictReader(handle)
            if row.get("Run", "").startswith("SRR")
        }


def build_accession_run_list(accessions: pathlib.Path, runinfo_dir: pathlib.Path, output: pathlib.Path) -> int:
    selected = [
        line.strip()
        for line in accessions.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    runs: set[str] = set()
    missing: list[pathlib.Path] = []
    for accession in selected:
        path = runinfo_dir / f"{accession}.runinfo.csv"
        if not path.exists():
            missing.append(path)
            continue
        runs.update(read_runs(path))

    if missing:
        for path in missing:
            print(f"missing {path}")
        return 1

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(sorted(runs)) + ("\n" if runs else ""))
    print(f"wrote {output} with {len(runs)} runs")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Maintenance commands for p-site offset prediction metadata."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    fetch = subparsers.add_parser(
        "fetch-runinfo",
        help="Fetch NCBI SRA RunInfo CSVs; requires network access.",
    )
    fetch.add_argument("manifest", type=pathlib.Path)
    fetch.add_argument("--outdir", type=pathlib.Path, default=DEFAULT_METADATA_DIR)
    fetch.set_defaults(func=lambda args: fetch_runinfo_manifest(args.manifest, args.outdir))

    extract = subparsers.add_parser("extract-runs", help="Extract SRR run IDs from RunInfo CSV files.")
    extract.add_argument("runinfo_dir", type=pathlib.Path)
    extract.add_argument("--outdir", type=pathlib.Path, default=DEFAULT_METADATA_DIR)
    extract.set_defaults(func=lambda args: extract_run_lists(args.runinfo_dir, args.outdir))

    build = subparsers.add_parser("build-run-list", help="Build a combined run list from selected accessions.")
    build.add_argument("accessions", type=pathlib.Path)
    build.add_argument("--runinfo-dir", type=pathlib.Path, default=DEFAULT_METADATA_DIR)
    build.add_argument("--output", type=pathlib.Path, required=True)
    build.set_defaults(
        func=lambda args: build_accession_run_list(args.accessions, args.runinfo_dir, args.output)
    )

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
