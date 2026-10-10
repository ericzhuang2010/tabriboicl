# Configuration

- `predictors.json` is the shared registry for all 11 predictors, including
  display names, input/output schemas, implementation paths, and native
  prediction commands. External-tool compatibility code remains under `tools/`.
- `p_site_offset_prediction_landmark_riboseq.tsv` is the full retained paper
  matrix: 11 predictors across one within-sample split, eight study holdouts,
  two organism holdouts, and five protocol-category holdouts.
- `predictors/tabriboicl.json` defines the complete predictor and distinguishes
  TabRiboICL from its underlying TabICL classifier.
- `samples.tsv` contains only public accessions and preprocessing metadata for
  the 31 Ribo-seq samples. Explicit protocol/library fields preserve the
  metadata used by the feature extractor without shipping raw RunInfo files.
- `dataset.yaml` is copied into the data bundle by `scripts/init_dataset.py`.

The `nuclease` column includes OTTR library preparation and unknown values;
it is the legacy code name for the manuscript's mixed protocol categories.
