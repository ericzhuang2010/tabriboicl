#!/usr/bin/env bash
set -euo pipefail

TARGET=""
YES=0

usage() {
  printf 'Usage: %s {dmel|scer} --yes\n' "$0"
}

if [[ $# -lt 1 ]]; then
  usage
  exit 2
fi

TARGET="$1"
shift

while [[ $# -gt 0 ]]; do
  case "$1" in
    --yes)
      YES=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      usage
      exit 2
      ;;
  esac
done

if [[ "$YES" -ne 1 ]]; then
  printf 'Refusing reference download without explicit approval. Re-run with --yes only after approval.\n' >&2
  exit 3
fi

case "$TARGET" in
  dmel)
    mkdir -p data/references/dmel
    curl --fail --location --retry 3 -o data/references/dmel/genome.fa.gz \
      https://ftp.ensemblgenomes.org/pub/metazoa/release-52/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa.gz
    curl --fail --location --retry 3 -o data/references/dmel/cdna.fa.gz \
      https://ftp.ensemblgenomes.org/pub/metazoa/release-52/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.32.cdna.all.fa.gz
    curl --fail --location --retry 3 -o data/references/dmel/ncrna.fa.gz \
      https://ftp.ensemblgenomes.org/pub/metazoa/release-52/fasta/drosophila_melanogaster/ncrna/Drosophila_melanogaster.BDGP6.32.ncrna.fa.gz
    curl --fail --location --retry 3 -o data/references/dmel/annotation.gtf.gz \
      https://ftp.ensemblgenomes.org/pub/metazoa/release-52/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.52.gtf.gz
    gunzip -f data/references/dmel/*.gz
    ;;
  scer)
    mkdir -p data/references/scer
    curl --fail --location --retry 3 -o data/references/scer/genome.fa.gz \
      https://ftp.ensembl.org/pub/release-115/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
    curl --fail --location --retry 3 -o data/references/scer/cdna.fa.gz \
      https://ftp.ensembl.org/pub/release-115/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz
    curl --fail --location --retry 3 -o data/references/scer/ncrna.fa.gz \
      https://ftp.ensembl.org/pub/release-115/fasta/saccharomyces_cerevisiae/ncrna/Saccharomyces_cerevisiae.R64-1-1.ncrna.fa.gz
    curl --fail --location --retry 3 -o data/references/scer/annotation.gtf.gz \
      https://ftp.ensembl.org/pub/release-115/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.115.gtf.gz
    gunzip -f data/references/scer/*.gz
    ;;
  *)
    printf 'Unsupported reference target: %s\n' "$TARGET" >&2
    exit 2
    ;;
esac
