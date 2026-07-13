# Native Binding Prediction Backend

MimicNeoAI provides a local binding-prediction backend shared by mutation-derived,
cryptic, and microbial antigen pipelines. The packaged pipeline configuration
keeps `pvactools` as the default backend. Set
`others.binding_prediction_backend: mimicneoai` explicitly to use the native
backend.

## Pipeline boundary

- Mutation-derived antigens use external pVACtools only for VCF conversion and
  WT/MT protein FASTA generation. MimicNeoAI builds mutation-covering peptide
  windows, predicts binding, and writes a pVACseq-compatible merged table.
- Cryptic and microbial antigens read their peptide FASTA directly, build
  de-duplicated peptide-HLA-algorithm tasks, and write a pVACbind-compatible
  merged table.
- The common runner does not perform antigen-specific biological filtering.

## Configuration

Common native-backend options:

```yaml
others:
  binding_prediction_backend: "mimicneoai"
  binding_prediction_algorithms: "MHCflurry MHCflurryEL MHCnuggetsI MHCnuggetsII NNalign NetMHCpan NetMHCpanEL NetMHCIIpan NetMHCIIpanEL"
  binding_prediction_workers: 8
```

Cryptic and microbial pipelines also support:

```yaml
  binding_prediction_max_task_rows: 5000000
  binding_prediction_force_large_samples: false
```

The task limit is evaluated after peptide-window construction but before
`binding_tasks.tsv` is materialized. Oversized samples retain epitope windows,
the task-count estimate, and a manifest. Prediction runs only when
`binding_prediction_force_large_samples` is explicitly enabled.

## Output layout

Mutation-derived output:

```text
07.binding_prediction_mimicneoai/
├── 00_input_vcf
├── 01_pvactools_sources
├── 02_epitope_tasks
├── 03_binding_predictions
├── 04_merged_epitopes
└── archive
```

Cryptic and microbial native output:

```text
<binding_step>_mimicneoai/
├── mimicneoai_epitope_tasks
├── mimicneoai_binding_predictions
├── combined
└── <sample>.mimicneoai_binding.summary.json
```

## Resume rules

- Mutation pVACtools source files record the VCF, flank length, pass-only mode,
  and pVACtools image identity. A later mismatch stops with an explicit request
  to use a new output directory instead of silently reusing converter/FASTA
  files.
- Non-mutation epitope windows are reused only when the peptide FASTA identity
  and requested peptide lengths match their manifest.
- Binding tasks are reused only when the epitope-window identity, HLA file, and
  algorithm sets match their manifest.
- Predictor chunks are reused only when normalized rows exactly match the
  current peptide, HLA, algorithm, MHC class, and peptide length requests.
- Changed inputs rebuild the affected stage instead of silently accepting stale
  output.

File identity uses resolved path, size, and nanosecond modification time. Use a
new output directory when replacing an input while preserving all three values.

## Output semantics

- `IC50_SUMMARY_ALGORITHMS` contains binding-affinity algorithms only.
- EL/presentation algorithms contribute their score and percentile fields but
  never contribute to Best/Median IC50 or fold-change summaries.
- Mutation missense and in-frame events retain corresponding WT predictions.
- Frameshift WT peptide, WT IC50, WT percentile, and fold-change fields remain
  blank.
- Unsupported predictor-HLA pairs remain in the normalized long table with
  `status=skipped` and `error=unsupported_allele_by_predictor`.
- Predictor failures remain in the long table with `status=error`.

## HLA support discovery

The runner builds support catalogs from the installed predictor resources:

- MHCflurry model `allele_sequences.csv` or `--list-supported-alleles`;
- NetMHCpan `-listMHC`;
- NetMHCIIpan `-list`;
- MHCnuggets model and training-resource intersection;
- IEDB MHC-I/MHC-II allele-info pickle files.

`binding_predictions.summary.json` records `allele_support_matrix` and
`predictor_runtime`, including catalog sources, allele counts, configured
executables, scripts, and environments. Catalog discovery failure is fail-open:
the task is attempted and the adapter records the real execution result.

## Regression tests

```bash
PYTHONPATH=. python -m unittest discover -v \
  mimicneoai/functions/binding_prediction/tests
```

The dependency-light suite does not invoke external predictors. It covers all
three pipeline backend branches, mutation event classes, cryptic and microbial
FASTA fixtures, HLA-II pairing, unsupported alleles, scale gating, resume input
signatures, normalized error states, and merged-table summary semantics.
