# HLA pseudo-sequence resources

This directory documents the optional NetMHC-derived HLA pseudo-sequence
resources used for immunogenicity model training.

The local resource files are intentionally kept under `local/`, which is
ignored by git, because the original NetMHCpan and NetMHCIIpan resources are
distributed by DTU under their own download terms.

Expected local layout:

```text
local/
  raw/
    netmhcpan/
      MHC_pseudo.dat
    netmhciipan/
      pseudosequence.2023.dat
  netmhcpan_class1_allele_to_pseudoseq.csv
  netmhciipan_class2_allele_to_pseudoseq.csv
```

Class I pseudo-sequences are derived from the official NetMHCpan data bundle
file `data/MHC_pseudo.dat`. Class II pseudo-sequences are derived from the
official NetMHCIIpan training resource file `pseudosequence.2023.dat`.

To regenerate the local CSV files after placing the raw files above:

```bash
python mimicneoai/immunogenicity_prediction/resources/hla_pseudoseq/prepare_netmhc_pseudoseq_csv.py
```

The generated CSV fields are:

```text
allele,pseudo_sequence,mhc_class,raw_allele,source_file
```

