# Immunogenicity Prediction Subtool

This subtool predicts immunogenicity scores for peptide-HLA pairs.
It can run as a standalone tool, and it can also be called by other pipelines.

## Run

```bash
# Method 1
python -m mimicneoai.immunogenicity_prediction.immunogenicity_prediction -c /path/to/immunogenicity_prediction_configure.yaml

# Method 2 (unified CLI)
mimicneoai immunogenicity-prediction -c /path/to/immunogenicity_prediction_configure.yaml
```

## Example Test Data

Example files are organized under:

```text
mimicneoai/example/immunogenicity_prediction/
├── config/
│   ├── microbial.yaml
│   ├── mutation_derived.yaml
│   └── cryptic.yaml
├── input/
│   └── input_peptide_hla.csv
├── models/
│   └── hla_prot.fasta  # retained for the legacy example
├── run_cryptic_ensemble.py
└── output/
    ├── microbial_predictions.csv
    ├── mutation_derived_predictions.csv
    └── cryptic_predictions.csv
```

### Input CSV format (`input/input_peptide_hla.csv`)

Required columns:
- `hla`
- `peptide`

Optional column:
- `hla_type` (not required for prediction; can be kept for annotation only)

Current example content:

```csv
hla,hla_type,peptide
HLA-A*02:01,HLA-I,ILDAIELAV
HLA-A*02:01,HLA-I,MLAAKTTVPV
HLA-A*02:07,HLA-I,ILDAIELAV
HLA-A*02:01,HLA-I,ALGNPEDFPV
HLA-A*02:01,HLA-I,ALGDWRAEV
HLA-A*02:01,HLA-I,GMLAAKTTV
```

## Source-specific examples

The runtime model files are recovered separately and are not included in the
repository. Place them under:

```text
mimicneoai/immunogenicity_prediction/models/default/
├── microbial/model.pth
├── mutation_derived/model.pth
└── cryptic/member_01.pth ... member_10.pth
```

Run the microbial and mutation-derived examples with the standard CLI:

```bash
python -m mimicneoai.immunogenicity_prediction.immunogenicity_prediction \
  -c mimicneoai/example/immunogenicity_prediction/config/microbial.yaml

python -m mimicneoai.immunogenicity_prediction.immunogenicity_prediction \
  -c mimicneoai/example/immunogenicity_prediction/config/mutation_derived.yaml
```

Run the cryptic example with its bundled ensemble runner:

```bash
python mimicneoai/example/immunogenicity_prediction/run_cryptic_ensemble.py \
  -c mimicneoai/example/immunogenicity_prediction/config/cryptic.yaml
```

For pipeline code, use `run_default_inference` with `microbial`,
`mutation_derived`, or `cryptic`. It resolves the recovered runtime payloads
and averages the ten cryptic members. `run_inference` remains available for an
explicit model path and for legacy configurations.

To verify a recovered runtime payload against a locked Fig. 2 prediction file
without copying manuscript data into this repository:

```bash
python -m mimicneoai.immunogenicity_prediction.tests.verify_locked_prediction_reproduction \
  --reference-tsv /path/to/three_antigen_specific_models.test_predictions.tsv \
  --antigen-class microbial --device cuda
```

The verifier defaults to an absolute score tolerance of `1e-3` to allow the
small R/float differences observed when recalculating peptide features across
runtime environments.

The legacy `microbial_pred_test.yaml` and the files in the example `models/`
directory are retained for backwards compatibility.

## Config Template

General template:
- `mimicneoai/configures/immunogenicity_prediction_configure.yaml`

Required keys:
- `path.input_csv`
- `path.model_path`
- `path.hla_fasta`
- `path.output_csv`

## Reuse in Other Pipelines

This module can also be reused by other pipelines through:

- `mimicneoai.functions.immunogenicity_runner.predict_immunogenicity_df`
- `mimicneoai.functions.immunogenicity_runner.predict_immunogenicity_csv`
