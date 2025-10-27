# MimicNeoAI

MimicNeoAI is a unified toolkit for discovering **microbial epitopes**, **cryptic epitopes**, and **mutation-derived neoantigens**, and for assessing **molecular mimicry** between microbial and tumor antigens.

It ships with three production pipelines:

- **Microbial pipeline** (`mimicneoai/microbial_pipeline`) – host‐read depletion, microbial profiling, microbial peptide extraction, HLA binding, immunogenicity.
- **Cryptic pipeline** (`mimicneoai/cryptic_pipeline`) – lncRNA/novel transcript reconstruction, sORF discovery, expression quantification, HLA binding, immunogenicity.
- **Mutation-derived pipeline** (`mimicneoai/mutation_derived_pipeline`) – QC, alignment, somatic variant calling/annotation, HLA typing, pVACseq-based neoantigen discovery.

> If you use MimicNeoAI, please cite the preprint listed in **Citation** below.

---

## Repository Layout
```bash
mimicneoai/
├─ configures/ # Example YAMLs for configuration and paths
│ ├─ cryptic_configure.yaml # cryptic pipeline example
│ ├─ microbial_configure.yaml # microbial pipeline example
│ ├─ mutation_derived_configure.yaml # mutation-derived pipeline example
│ └─ paths.yaml # common paths example
├─ cryptic_pipeline/ # cryptic (sORF) pipeline
├─ microbial_pipeline/ # microbial pipeline
└─ mutation_derived_pipeline/ # mutation-derived pipeline
```

---

## Quick Start

Each pipeline has its own README with **installation**, **configuration**, and **run** instructions:

- Microbial: [`mimicneoai/microbial_pipeline/README.md`](mimicneoai/microbial_pipeline/README.md)
- Cryptic: [`mimicneoai/cryptic_pipeline/README.md`](mimicneoai/cryptic_pipeline/README.md)
- Mutation-derived: [`mimicneoai/mutation_derived_pipeline/README.md`](mimicneoai/mutation_derived_pipeline/README.md)

Reference bundles and minimal test data:
- **Zenodo**: https://doi.org/10.5281/zenodo.15582924

---

## High-Level Features

- End-to-end automation with resumable steps and structured logs
- Modular YAML configs (sample list, runtime args, tool/resource paths)
- HLA binding and **immunogenicity** scoring (multi-tool support)
- Consistent output layout for downstream integration and figure generation

---

## Requirements

- Linux (tested on Ubuntu 20.04 x86_64)
- Conda or Mamba (recommended)
- Toolchain per pipeline (see pipeline READMEs), including:
  - `fastp`, `samtools`, `bwa/bowtie2`
  - `blast+` (microbial)
  - `GATK`/`VEP` (mutation-derived)
  - `HLA-HD` (HLA typing)
  - `pVACtools`, IEDB predictors, and optional ML predictors (MHCflurry, MHCnuggets, BigMHC, DeepImmuno)

---

## Configuration

Copy and edit YAMLs under `mimicneoai/configures/`:
- `*_configure.yaml`: runtime toggles, input/output roots, sample list
- `paths.yaml`: absolute paths to references and executables

Prefer absolute paths to avoid ambiguity.

---

## Reproducibility & Logging

- Each step writes command and status logs under the pipeline’s working directory.
- On reruns, **existing non-empty outputs are skipped**. If in doubt, delete incomplete products from the affected step directory before rerun.

---

## License

See [`LICENSE`](LICENSE).

---

## Citation

**MimicNeoAI: An integrated pipeline for identifying microbial epitopes and mimicry of tumor neoepitopes**  
Tao Chen, Wei Wang, Xiao Zuo, Yuxin Zhang, Mingwei Li, Zhilei Li, Yin He, Yanfei Zhou, Fang Ye, Bin Zhang, Qionghui Jiang, Huimin Liu, Lu Zhang, Jinman Fang, Yuanwei Zhang  
*bioRxiv* 2025.06.13.658292; doi: https://doi.org/10.1101/2025.06.13.658292
