# Molecular Mimicry

The MimicNeoAI molecular-mimicry module identifies equal-length, HLA class
I-length peptide pairs with sequence similarity concentrated in the central,
putative TCR-facing region. It operates on final peptide Core outputs from two
or more antigen workflows and does not rerun antigen discovery. In this
documentation, a peptide Core is the unique candidate set retained after a
source workflow's formal quality-control steps.

## Evidence Boundary

The current public policy is `sequence_mimicry_v1.2`. A pair is evaluated only
when both peptides:

- belong to the same patient but different antigen sources;
- have the same length of 8, 9, 10, or 11 amino acids;
- contain only the 20 canonical amino acids; and
- passed their source-specific upstream Core QC.

A pair is called `sequence_based_mimicry_v12_candidate` when all three
sequence requirements are met:

1. central position-wise identity is at least 60% over P4 through
   Pmin(8, L-1);
2. the same central window contains at least three consecutive exact matches;
3. the full peptide contains at least four position-wise exact matches.

HLA type, predicted binding, immunogenicity, mass-spectrometry evidence, and
the legacy longest common substring do not determine this primary sequence
call. Sequence similarity alone does not demonstrate shared HLA presentation
or T-cell cross-reactivity.

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

Thresholds are part of the named frozen policy and are intentionally absent
from the run configuration. A change to an eligibility threshold requires a
new policy version.

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

`mimicry_pairs.tsv.gz` contains only pairs passing v1.2. The stagewise table
records the number of all evaluable, same-length cross-source pairs without
materializing failed Cartesian pairs. Member provenance is written once for
each participating source occurrence and retains the original source row as
JSON.

For source pairs containing mutation-derived peptides, event-level WT
classification is reported separately as `MT_specific_sequence_mimic`,
`WT_compatible_sequence_mimic`, `WT_not_evaluable`, or
`WT_not_evaluable_for_neo_specificity`.

The manifest records resolved policy values, source file identities, code
identities, output hashes, worker settings, and elapsed time. Existing results
are reused only when the complete input and output contracts match.

## Shared HLA-I Evidence

Set `binding_support.enabled: true` to annotate sequence pairs using existing
MimicNeoAI merged binding tables listed in `binding_path`. The module does not
run or reimplement a predictor. It uses the Best/Median values produced by the
shared binding backend and reports support only when both peptide members pass
for the same classical HLA-A, HLA-B, or HLA-C allele:

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
`binding_not_evaluable_hla_contract_mismatch`. This annotation never removes
or changes a v1.2 sequence candidate.

## Python API

```python
from mimicneoai.mimicry import compute_pair_metrics, is_sequence_mimic

is_sequence_mimic("AAAAAAAA", "CAAAAACC")
metrics = compute_pair_metrics("AAAAAAAA", "CAAAAACC")
```

The stable public API is exported from `mimicneoai.mimicry`. Modules beginning
with an underscore are implementation details.
