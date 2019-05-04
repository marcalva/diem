
# DIEM v0.1

Clean out debris-containing droplets from single-cell based data, particularly from 
single-nucleus RNA-seq of frozen tissue.

## Installation

Use devtools to install as an R package

```R
library(devtools)
devtools::install()
```

## Usage

Shown below is the suggested workflow.

```R
library(diem)

x <- read_10x(path)
sce <- create_SCE(x)
sce <- subset_dropls(sce, min_c=30)
sce <- normalize(sce)
sce <- get_var_genes(sce)
sce <- get_pcs(sce, n_pcs=5)
sce <- get_bg_cand(sce, bg_max=150)
sce <- run_em_pcs(sce, n_pcs=5)
sce <- call_targets(sce)
sce <- summary_results(sce)
```

