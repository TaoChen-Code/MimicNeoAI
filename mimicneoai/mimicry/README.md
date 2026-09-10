# Molecular Mimicry

The MimicNeoAI molecular-mimicry module compares quality-controlled 8–11-aa
candidate peptides from different antigen sources within each patient. It first
identifies sequence-similar pairs on the basis of positional identity across the
central region and across the full peptide. When HLA-binding predictions are
available, pairs in which both peptides are predicted to bind the same patient
HLA-I allele are classified as predicted mimicry candidates. The module operates
on candidate sets retained by the antigen-discovery workflows and does not rerun
antigen discovery.

## Sequence Matching and Shared HLA-I Binding

The analysis distinguishes sequence-similar pairs from predicted mimicry
candidates. A pair is evaluated for sequence similarity only when both peptides:

- belong to the same patient but different antigen sources;
- have the same length of 8, 9, 10, or 11 amino acids;
- contain only the 20 canonical amino acids; and
- have passed source-specific quality control.

A peptide pair is classified as a sequence-similar pair when all three
requirements are met:

1. central position-wise identity is at least 60% over P4 through
   Pmin(8, L-1);
2. the same central window contains at least three consecutive exact matches;
3. the full peptide contains at least four position-wise exact matches.

Passing pairs are reported as sequence-similar pairs. The machine-readable
`primary_call` field records this result as
`sequence_based_mimicry_v12_candidate`. HLA type, predicted binding,
immunogenicity, immunopeptidomics evidence, and the legacy longest common
substring do not contribute to the sequence classification. Exact P2 and
PΩ matches and BLOSUM62 scores are retained as descriptive annotations.

Microbial--Mutation and Microbial--Cryptic comparisons comprise the
microbial--host analyses. Cryptic--Mutation comparisons are reported separately
as cross-source tumour-antigen sequence similarity. Sequence similarity and
shared predicted HLA-I binding nominate peptide pairs for further evaluation
but do not establish natural co-presentation or T-cell cross-reactivity.

## Command Line

Create a compact run configuration:

```yaml
method: "sequence_mimicry_v1.2"
input_manifest: "/analysis/mimicry_inputs.tsv"
output_dir: "/analysis/mimicry"

source_pairs:
  - "microbial--mutation-derived"
  - "microbial--cryptic"
  - "cryptic--mutation-derived"

workers: 8
shard_count: 2
```

Then run:

```bash
mimicneoai mimicry -c mimicry.yaml
```

The method identifier records the implemented thresholds for reproducibility,
so these values are not repeated in the run configuration. Changes to the
eligibility thresholds require a new method version.

`workers` limits concurrent search processes. `shard_count` divides each
patient, source-pair, and peptide-length comparison into deterministic,
resumable units. More shards may improve utilization for a single large
comparison, but each shard rebuilds its sequence index. A value of 1 or 2 is
therefore preferable for multi-patient analyses; increase it only when one
comparison does not provide enough parallel work.

## Input Manifest

The tab-delimited input manifest has one row per patient and antigen source:

| patient_id | antigen_source | input_format | peptide_path | provenance_path | run_manifest_path | binding_path |
|---|---|---|---|---|---|---|
| Patient_1 | microbial | microbial_core_v1 | `<microbial-core.tsv>` | `<parent-map.tsv>` | `<run-manifest.json>` | `<binding.tsv>` |
| Patient_1 | cryptic | cryptic_final_core_v1_1 | `<cryptic-core.tsv>` | `<final-sidecar.tsv>` | `<run-manifest.json>` | `<binding.tsv>` |
| Patient_1 | mutation-derived | mutation_epitope_windows_v1 | `<epitope-windows.tsv>` | | | `<binding.tsv>` |

Relative paths are resolved from the manifest directory. Supported formats
are:

- `microbial_core_v1`: paired microbial Core records with
  `peptide_qc_status=core`;
- `cryptic_final_core_v1_1`: final external-normal-QC records with
  `peptide_core_status=core`;
- `mutation_epitope_windows_v1`: mutation-covering MT windows and event-level
  WT controls;
- `generic_peptide_table_v1`: a portable table containing `peptide`,
  `peptide_id`, `mhc_class`, and `qc_status=core`.

For mutation input, the generic format may additionally contain `event_id`,
`wt_peptide`, and `wt_control_status`. Native mutation output uses
`mt_epitope_seq`, `wt_epitope_seq`, and `covers_mutation` directly.

When a provenance table is supplied, every eligible peptide must map to at
least one provenance row. Explicitly incomplete upstream manifests and Core
manifests with `binding_eligible=false` fail closed.

## Outputs

```text
mimicry/
├── mimicry_pairs.tsv.gz
├── mimicry_mutation_wt_evidence.tsv.gz
├── mimicry_member_provenance.tsv.gz
├── mimicry_input_qc.tsv
├── mimicry_stagewise_qc.tsv
├── run_manifest.json
└── .parts/
```

`mimicry_pairs.tsv.gz` contains the sequence-similar pairs identified by the
analysis. When binding support is enabled, their shared HLA-I binding status
and predicted mimicry classification are recorded in
`mimicry_hla_support.tsv.gz`. The stagewise table records the number of all
evaluable, same-length cross-source pairs without materializing failed
Cartesian pairs. Member provenance is written once for each participating
source occurrence and retains the original source row as JSON.

For source pairs containing mutation-derived peptides, event-level WT
classification is reported separately as `MT_specific_sequence_mimic`,
`WT_compatible_sequence_mimic`, `WT_not_evaluable`, or
`WT_not_evaluable_for_neo_specificity`.

The manifest records resolved policy values, source file identities, code
identities, output hashes, worker settings, and elapsed time. Existing results
are reused only when the complete input and output contracts match.

## Predicted Mimicry Candidates with Shared HLA-I Binding

Set `binding_support.enabled: true` to evaluate sequence-similar pairs using
existing MimicNeoAI merged binding tables listed in `binding_path`. The module
does not run or reimplement a predictor. A pair is classified as a predicted
mimicry candidate only when both peptides meet the binding criteria for the
same patient HLA-A, HLA-B, or HLA-C allele:

```yaml
binding_support:
  enabled: true
```

- Best and Median IC50 are both below 500 nM;
- Best and Median percentile rank are both below 2.

The result is written to `mimicry_hla_support.tsv.gz`. Missing metrics or
unsupported alleles are marked `binding_not_evaluable`, not as negative
binding evidence. The HLA-A/B/C allele universe must also agree between both
source results for a patient. A mismatch is reported as
`binding_not_evaluable_hla_contract_mismatch`. This assessment does not change
the set of sequence-similar pairs.

## Python API

```python
from mimicneoai.mimicry import compute_pair_metrics, is_sequence_mimic

is_sequence_mimic("AAAAAAAA", "CAAAAACC")
metrics = compute_pair_metrics("AAAAAAAA", "CAAAAACC")
```

The stable public API is exported from `mimicneoai.mimicry`. Modules beginning
with an underscore are implementation details.
