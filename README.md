## Introduction
This is the code repository for a paper titled 
"TabRiboICL for read-specific P-site prediction and a unified framework for benchmarking P-site assignment in ribosome"
by Eric Zhuang and Manoj Kandpal.

A preprint is available at https://zenodo.org/records/24237397

## Install

Use Python 3.11 or 3.12 and run commands from the repository root:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[predictor,comparators]'
tabriboicl --help
```

Python dependencies are in [pyproject.toml](pyproject.toml). Raw-read processing
also needs the bioinformatics programs listed in [environment.yml](environment.yml);
the Linux setup below installs both.

## Predict and evaluate

Replace the input paths below with your prepared Parquet candidate tables.
Training requires candidate labels; evaluation also requires labelled test
candidates. Generate these tables with `tabriboicl extract` and `tabriboicl split`
after preparing alignments and references, as shown below.

```bash
tabriboicl predict --model tabriboicl \
  --train data/train.parquet \
  --test data/test.parquet \
  --output results/tabriboicl/candidate_predictions.parquet \
  --top-output results/tabriboicl/read_predictions.parquet \
  --max-train-reads 13000 \
  --kv-cache false

tabriboicl evaluate \
  --test data/test.parquet \
  --predictions results/tabriboicl/candidate_predictions.parquet \
  --metrics-output results/tabriboicl/metrics.tsv
```

TabICL may download pretrained weights on first use. Each command accepts
`--help`. Use `tabriboicl baselines` for fixed-offset and read-length baselines;
`predict --model` also supports `xgboost` and `lightgbm`.

## Full paper benchmark

Use Linux with Rscript, Git, curl, tar, and Python already installed. This setup
replaces the virtual-environment setup above and installs the external-tool
dependencies. The data preparation step downloads the configured sequencing runs.

```bash
bash scripts/setup_linux_p_site_offset_prediction.sh --download-refs
source scripts/local_env.sh

mkdir -p external_tools/src
git clone https://github.com/nabeel-bioinfo/Asite_IP_method \
  external_tools/src/Asite_IP_method
export ASITE_IP_SRC="$PWD/external_tools/src/Asite_IP_method"

python scripts/init_dataset.py
bash scripts/run_p_site_offset_prediction_data_prep.sh --threads 8

tabriboicl extract \
  --sample-table data/datasets/p_site_offset_prediction/landmark_riboseq/samples.tsv
tabriboicl split

TABICL_KV_CACHE=false bash scripts/run_p_site_offset_prediction_evals.sh
```

