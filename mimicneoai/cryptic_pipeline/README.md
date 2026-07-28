# Cryptic (sORF-Encoded) Antigen Pipeline

Discover and prioritize sORF-encoded peptides from known and novel transcripts.

## Overview
1. QC
2. STAR alignment
3. Known/novel lncRNA-sORF discovery
4. Tumor/control quantification
5. HLA typing
6. Aberrantly expressed sORF peptide extraction
7. ORF genome annotation and ORF-level filtering
8. Cryptic Core QC
9. Optional junction QC inside Cryptic Core v1.1
10. External-normal exact sequence QC
11. Binding prediction

## Installation

```bash
pip install mimicneoai
# or
conda install -c conda-forge -c bioconda mimicneoai
```

## External Dependencies (in `PATH`)

| Tool | Recommended/validated version | Purpose |
|---|---|---|
| `fastp` | `0.22.0` | FASTQ QC |
| `STAR` | `2.5.3a` | RNA alignment |
| `samtools` | `1.5` | BAM/SAM/FASTQ processing |
| `stringtie` | `3.0.1` | Transcript assembly (novel branch) |
| `gffcompare` | `0.12.10` | Transcript annotation comparison |
| `gffread` | `0.12.7` | GTF/FASTA conversion |
| `minimap2` | `2.30-r1287` | Contig-to-reference alignment |
| `bcftools` | `1.11` | Variant calling/filtering/consensus (known branch) |
| `TransDecoder.LongOrfs` | `5.5.0` | ORF calling |
| `salmon` | `1.10.0` | Expression quantification |
| `bowtie2` | `2.4.1` | HLA typing pre-alignment |
| `hlahd.sh` | `1.7.0` | HLA typing |
| `apptainer` | `1.4.2` | Trinity/pVACtools container execution |
| `pVACtools` (container) | `4.2.1` | Binding prediction (`pvacbind`) |
| `awk` | system tool | FASTQ header normalization in HLA typing |

Optional (only when `others.trinity_mode: native`):
- `Trinity`

## Database and Paths

Shared documentation:
- [`mimicneoai/configures/Database_and_Paths.md`](../configures/Database_and_Paths.md)

## Configuration

Use the canonical template directly:
- [`mimicneoai/configures/cryptic_configure.yaml`](../configures/cryptic_configure.yaml)

Key switches are defined in this template (for example `others.known`, `others.novel`, `others.salmon_quant_control`, `others.hla_binding_pred`). Please follow template key names exactly.

## Run

```bash
# Method 1
python -m mimicneoai.cryptic_pipeline.cryptic -c /path/to/cryptic_configure.yaml

# Method 2 (unified CLI)
mimicneoai cryptic -c /path/to/cryptic_configure.yaml
```

## Output Structure

Pipeline outputs are written under:

```text
<output_dir>/Cryptic/<tumor_sample>/
├── 00-clean
├── 01-star
├── 02-known
├── 03-novel
├── 04-salmon_quant
├── 05-hla_typing
├── 06-aeSEPs
├── 07-orf_genome_annotation
├── 08-orf_filter
├── 08b-cryptic_core_qc
├── 08c-external_normal_qc
├── 09-hla_binding_pred_mimicneoai
├── 10-immunogenicity_prediction_mimicneoai
└── 023-shared
```

Notable subfolders:
- `04-salmon_quant/salmon_index`, `salmon_quant`, `salmon_quant_control`
- `01-star/star-provenance-freeze`: optional frozen STAR provenance for
  paired junction QC. The directory name is stable; the policy version is
  recorded inside the manifests.
- `07-orf_genome_annotation`: maps selected ORF/CDS records back to the reference genome
- `08-orf_filter`: writes the ORF-filtered aeSEP FASTA used by binding prediction
- `08b-cryptic_core_qc`: materializes the strict pre-binding Cryptic Core
  with parent sidecars, HLA-I/HLA-II peptide-core FASTA files, stagewise counts,
  and an input/config manifest.
- `08c-external_normal_qc`: applies frozen external-normal resources. Policy
  `cryptic_external_normal_qc_v1.0` uses exact peptide sequence plus HLA class.
  Policy `cryptic_external_normal_qc_v1.1` keeps those exact-match rules and
  adds normal smORF coordinate/frame evidence. It keeps the complete source-Core
  landscape and writes the tumor-restricted primary Core FASTA used by
  downstream binding.
- `09-hla_binding_pred/<tumor_sample>/pvacbind` for the pVACbind backend, or
  `09-hla_binding_pred_mimicneoai/<tumor_sample>/` for the MimicNeoAI backend

## Notes

- Tumor/control samples should be provided as `Tumor,Control` in `samples`.
- The pipeline is resumable when manifest, input, code, configuration, and
  output signatures match. Empty Core FASTA files are valid for zero-candidate
  samples and can be reused through the manifest.
- `others.orf_genome_annotation` and `others.orf_filter` default to enabled;
  binding prediction uses the ORF-filtered aeSEP FASTA when `orf_filter` is enabled.
- `others.cryptic_core_qc` defaults to enabled in the template. It accepts only
  `novel` and `noncoding` aeSEP sources, requires ORF-filtered parent records,
  preserves excluded rows, and writes pre-tiled peptide-core FASTA files for
  HLA-I 8-11 aa and HLA-II 13-17 aa.
- `candidate_selection.mode: all` keeps the complete Cryptic Core and is the
  default. `ranked_cap` first applies strict source/expression, coordinate,
  reference-translation, junction, and human-reference QC at the parent/peptide
  levels, then ranks parent ORFs and caps unique peptide sequences independently
  for HLA-I and HLA-II, for example `max_hla_i_peptides: 400000` and
  `max_hla_ii_peptides: 400000`. Cap-external parents or peptides are recorded
  as `not_selected_due_to_analysis_cap` / deferred; this is a compute-routing
  state, not binding-negative evidence or QC failure. Ranked mode keeps a
  complete, deterministic ranked parent stream and materializes only the
  selected/boundary peptide candidates in 08b. When 08c removes selected
  candidates, it refills by continuing that ranked parent stream, applying the
  same human-reference, coordinate, junction and external-normal evidence
  checks to each refill peptide before it can enter the final peptide Core.
- `junction_qc.enabled: true` enables production junction support QC for
  `cryptic_core_qc_v1.1`. It consumes the explicit STAR pair table produced by
  `freeze_star_provenance.py`; downstream code should read SJ paths from that
  table rather than infer unprefixed filenames. The primary policy is
  `junction_qc_v1.0`, requiring all required parent junctions to have tumor
  unique split reads >=2. It also records threshold sensitivity for 1, 2, 3,
  and 5 reads. Intronless parents are retained as not applicable. Matched-normal
  junctions are annotation-only and are reported as unique-read threshold
  fields, not as full-chain support.
  `cryptic_peptide_junction_evidence.tsv` covers retained Core peptide windows
  with stable peptide IDs; human-reference-excluded peptide windows remain in
  the QC/excluded tables and are not assigned peptide-level junction evidence.
- `rna_variant_editing_qc.enabled` remains `false` for the current strict main
  analysis. Reference-translation mismatches and RNA-dependent parents are kept
  in provisional/excluded sidecars and do not enter this round's strict primary
  Core. The optional `cryptic_rna_variant_editing_qc_v1.0` code can validate a
  normalized known-branch RNA VCF and AD-derived read evidence for exploratory
  RNA-supported sequence reconstruction, but REDIportal is not a required
  pipeline dependency in the main configuration and is not used to label events
  as somatic, non-editing, or tumor-specific mutation. If this optional branch
  is enabled without a formal editing resource/provenance contract, it writes
  `run_status=complete_exploratory`, `binding_eligible=false`, and cannot be
  routed into binding.
- `others.cryptic_external_normal_qc` defaults to disabled in the generic
  template because the frozen resource package is project-specific. Formal
  project configs should enable it and provide `external_normal_resources`, or
  rely on populated `paths.yaml` resource entries. Missing or hash-mismatched
  resources fail closed. `allow_missing_external_normal_resources: true` is
  exploratory only and cannot be routed into binding.
- External-normal QC v1.1 only uses coordinate/frame evidence when the candidate
  parent has exactly one primary ORF-genome alignment in the complete
  `orf2genome.bam`, no secondary/supplementary ambiguity, MAPQ above the
  configured threshold, canonical GRCh38 contig coordinates, consistent
  CDS/block length, and strand-aware genomic translation from the configured
  GRCh38 FASTA reproduces the parent peptide. Reference translation mismatch is
  recorded as RNA-variant-aware not evaluable, not as normal evidence.
- In v1.1, a unique peptide is removed from the strict primary Core only when at
  least one trusted candidate parent is
  `normal_smorf_coordinate_frame_concordant`. Partial overlap, frame-discordant
  overlap, incompatible junction chains, low MAPQ, secondary/supplementary
  alignments, noncanonical contigs, and reference-translation mismatch are
  annotated as not evaluable or non-excluding coordinate evidence.
- `external_normal_status` distinguishes exact-match exclusions,
  coordinate/frame exclusions, and peptides with both evidence types; use
  `external_normal_qc_reasons` for the specific resource-level reason.
- `others.binding_prediction_backend` defaults to `mimicneoai` in the template.
  The local backend estimates task scale before materializing the task table;
  see the [native binding backend documentation](../functions/binding_prediction/README.md).
- For legacy pVACtools reruns, disable `others.cryptic_core_qc` or explicitly
  provide a parent FASTA. The strict peptide-core output is routed only to the
  MimicNeoAI backend to avoid accidental second-pass peptide tiling.
- When `08c` is enabled, binding consumes
  `08c-external_normal_qc/cryptic_tumor_restricted_primary_core.fasta`.
  It does not fall back to `08b`, `06-aeSEPs`, or the ORF-filtered parent FASTA.
  Binding is allowed only when the 08c manifest is formal complete,
  `binding_eligible=true`, and the final FASTA size/SHA256 matches
  `final_binding_fasta_identity`.
- `08c-external_normal_qc/cryptic_final_peptide_parent_sidecar.tsv` is the
  pre-binding downstream sidecar. It maps each retained unique peptide back to
  all supporting parent windows and records expression, coordinate/reference,
  junction, human-reference, external-normal, ranked-selection, and discovery
  FDR-status fields. It intentionally contains MHC class only; HLA allele and
  binding evidence are joined after step 09.
- With `others.binding_prediction_backend: mimicneoai`, set
  `others.binding_prediction_preset: full` for one-stage multialgorithm
  prediction or `fast` for EL-rank Stage 1 routing before formal local binding
  prediction. Omitting the preset preserves explicit length and algorithm
  settings in the YAML.
- Set `others.run_immunogenicity_prediction: true` to score peptide-HLA rows
  after binding. If binding is skipped by the scale gate, immunogenicity scoring
  is skipped and recorded in the summary rather than being treated as negative.
- Immunogenicity inference requires a Python environment with PyTorch and
  scikit-learn. Use `others.immunogenicity_python_bin` or the
  `MIMICNEOAI_IMMUNOGENICITY_PYTHON_BIN` environment variable to select a CPU
  or GPU runtime explicitly; otherwise the current pipeline Python is used.
