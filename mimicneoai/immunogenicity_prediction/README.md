# Immunogenicity Prediction Subtool

This subtool predicts immunogenicity scores for peptide-HLA pairs.

## Run

```bash
# Method 1
python -m mimicneoai.immunogenicity_prediction.immunogenicity_prediction -c /path/to/immunogenicity_prediction_configure.yaml

# Method 2 (unified CLI)
mimicneoai immunogenicity-prediction -c /path/to/immunogenicity_prediction_configure.yaml
```

## Config Template

Use:

- `mimicneoai/configures/immunogenicity_prediction_configure.yaml`

Required inputs:

- `path.input_csv`
- `path.model_path`
- `path.hla_fasta`
- `path.output_csv`

This module can also be reused by other pipelines through:

- `mimicneoai.functions.immunogenicity_runner.predict_immunogenicity_df`
- `mimicneoai.functions.immunogenicity_runner.predict_immunogenicity_csv`
