# Mutation Epitope Prediction

Experimental workflow for mutation-derived epitope generation and binding
prediction.

The current production pipeline runs pVACseq end to end under
`07.binding_prediction`. This package is a side-channel implementation used to
test a MimicNeoAI-native prediction workflow:

1. Use external pVACtools only for VCF-derived mutation annotation and WT/MT
   protein sequence generation.
2. Build stable variant and transcript identifiers independent of pVACseq
   chunk-local `Index` values.
3. Generate mutation-covering epitope windows, extended peptides, and
   de-duplicated peptide-HLA-algorithm prediction tasks.
4. Run binding predictors through MimicNeoAI adapters and merge prediction
   outputs back to variant, transcript, and expression
   evidence.

This package must not import or vendor pVACtools source code. pVACtools should
be invoked as an external tool, typically through the validated Apptainer image.

Prototype output directory layout:

```text
07.binding_prediction_v2/
  00_input_vcf/
  01_pvactools_sources/
  02_epitope_tasks/
  bp/
  logs/
```

Planned CLI steps:

```text
00_prepare_pvacseq_sources.py
01_build_epitope_tasks.py
02_predict_and_merge.py
```

For the prototype layout, run step 00 with the parent output directory:

```text
07.binding_prediction_v2/
```

Run step 01 with:

```text
-o 07.binding_prediction_v2/02_epitope_tasks
```
