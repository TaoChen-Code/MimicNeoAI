# Runtime immunogenicity model payloads

The runtime model payloads are intentionally not tracked by Git. Recover them
from the separately archived model bundle into this layout:

```text
default/
  microbial/model.pth
  mutation_derived/model.pth
  cryptic/member_01.pth ... member_10.pth
```

The expected SHA256 values are:

```text
microbial/model.pth                 0690b4a3ebf6d79ca7e8186e694b68651e44f1e1ea30d01ae84d00c5807fa9e3
mutation_derived/model.pth          6a55541a327b680013e3afb7cccdb14bdd9a7b0b8d100595492e52267a5ea013
cryptic/member_01.pth               7b9cc914610ca85424ea27ae783a7c1b233506023d8f02d214c5b060b5b41d18
cryptic/member_02.pth               fecc2bfe64f2eb6906cd06e1d61956bd143a91af814014ee17d07248410d96fd
cryptic/member_03.pth               6662dfc5b0c54f9151e83bdd77afbdcb5a410cdb22cc9df9ef846d5c7235b58e
cryptic/member_04.pth               21db89a6d2e4c8e86959222f04a1d574e9024beb9108890c68e894ea2bdde899
cryptic/member_05.pth               7c059e64ba1888409622f9d50f03053be813f484bc93287fedc13725136312dc
cryptic/member_06.pth               bb7a9e313a663087fa49c42e05cc76f969dac4bce84309645e30362707a8301c
cryptic/member_07.pth               8ddbe80dd82c3fa073f4e8214998895d87fbb4f6059d57d58d76b9e2f59b6928
cryptic/member_08.pth               3b6686a6ee60843cbc8a8c48674160aa22f5ce66c20595d50a23fb68e51fd7f4
cryptic/member_09.pth               258cb40fd82da568a18e25a787c9d998123d74867759df452fd36e1f31c114c9
cryptic/member_10.pth               1fceda96ee7146c918ebb855e2134704468d4eb30b8497cda16b785b7588e8b1
```

The source-specific runtime checkpoints require the separately supplied
NetMHC-derived HLA pseudo-sequence CSVs under
`../resources/hla_pseudoseq/local/`; see that directory's README.
