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
│   └── microbial_pred_test.yaml
├── input/
│   └── input_peptide_hla.csv
├── models/
│   ├── hla_prot.fasta
│   ├── MimicNeoAI_Microbial_Pred.pth
│   ├── MimicNeoAI_Cryptic_Pred.pth
│   └── MimicNeoAI_Mutation-derived_Pred.pth
└── output/
    └── predictions.csv
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

## Test Config (Microbial_Pred)

For immunogenicity testing, use the `Microbial_Pred` model:

- `mimicneoai/example/immunogenicity_prediction/config/microbial_pred_test.yaml`

Config content:

```yaml
path:
  input_csv: "mimicneoai/example/immunogenicity_prediction/input/input_peptide_hla.csv"
  output_csv: "mimicneoai/example/immunogenicity_prediction/output/predictions.csv"
  model_path: "mimicneoai/example/immunogenicity_prediction/models/MimicNeoAI_Microbial_Pred.pth"
  hla_fasta: "mimicneoai/example/immunogenicity_prediction/models/hla_prot.fasta"
  export_onnx: ""

args:
  batch_size: 512
  num_processes: 8
  verbose: true
  device: "auto"
  export_onnx_only: false
  onnx_opset: 17

io:
  peptide_col: "peptide"
  hla_col: "hla"
  score_col: "immunogenicity_score"
```

Run with:

```bash
python -m mimicneoai.immunogenicity_prediction.immunogenicity_prediction \
  -c mimicneoai/example/immunogenicity_prediction/config/microbial_pred_test.yaml

# or
mimicneoai immunogenicity-prediction \
  -c mimicneoai/example/immunogenicity_prediction/config/microbial_pred_test.yaml
```

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
