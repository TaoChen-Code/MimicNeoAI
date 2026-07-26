# Immunogenicity Prediction Subtool

This subtool predicts immunogenicity scores for peptide-HLA pairs.
It can run as a standalone tool, and it can also be called by other pipelines.

## Module boundary

- `api.py` and `runtime/` are the stable runtime surface for pipelines.
- `cli.py` is the YAML-configured command-line entry point;
  `immunogenicity_prediction.py` remains a compatibility module.
- `core.py` remains the compatibility implementation during the incremental
  migration; do not add new pipeline dependencies to it directly.
- `train_immunogenicity.py` and `benchmark/` are training/research tools, not
  required by formal inference.

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

Set `args.include_input_qc: true` to retain input QC annotations in output:
normalized HLA, peptide length, non-canonical amino-acid flags, duplicate
peptide-HLA flags, and HLA resource-match status. This augments output only;
it does not rewrite or silently exclude input records.

When called from a source pipeline, resolve `num_processes` with
`resolve_immunogenicity_num_processes(configure)`. It inherits the pipeline
YAML's `args.thread` or `args.threads`.

Pipeline entry points run immunogenicity scoring through a separate Python
interpreter when `others.immunogenicity_python_bin` is set, or when
`MIMICNEOAI_IMMUNOGENICITY_PYTHON_BIN` is exported. Use this for dedicated CPU
and GPU PyTorch environments. If neither is set, the current pipeline Python is
used.

The helper script below creates a dedicated runtime environment without
modifying the main pipeline Python:

```bash
# GPU runtime, CUDA 11.8 PyTorch wheel.
scripts/install_immunogenicity_runtime.sh gpu-cu118 /workspace/pkgs/MimicNeoAI_immunogenicity_gpu

# CPU fallback runtime.
scripts/install_immunogenicity_runtime.sh cpu /workspace/pkgs/MimicNeoAI_immunogenicity_cpu
```

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

For formal source-class inference, prefer:

- `mimicneoai.functions.immunogenicity_runner.predict_default_immunogenicity_df`
- `mimicneoai.functions.immunogenicity_runner.predict_default_immunogenicity_csv`

These select the recovered default weights, use checkpoint-compatible local
pseudo-sequences by default, emit input QC annotations, and apply the cryptic
ensemble automatically.
