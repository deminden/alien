# Analysis examples

R code and report walkthroughs behind the figures in the documentation: a female-versus-male contrast in GTEx
thyroid and TCGA lung, scored with GMTs built from the study configs in [examples/](../). The scripts build any
missing GMT with the `alien` CLI. They are examples, not part of the `bioalien` package, and are excluded from the
source distribution.

```bash
Rscript examples/analysis/gtex_thyroid_sex_contrast.R --help
Rscript examples/analysis/tcga_lung_sex_contrast.R --help
```

Reports: [GTEx thyroid](sex_contrast/gtex_thyroid/analysis.md), [TCGA lung](sex_contrast/tcga_lung/analysis.md).
