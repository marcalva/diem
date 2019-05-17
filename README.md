
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
sce <- set_limits(sce, min_bg_count=0, max_bg_count=100, min_tg_count=200, min_tg_gene=200)
sce <- get_de_genes(sce, n_genes=1000)
sce <- set_expression(sce)
sce <- normalize(sce)
sce <- get_bgscore(sce)
sce <- get_pcs(sce)
sce <- run_em_bg_score_pcs(sce, n_pcs=5, verbose=TRUE)
sce <- call_targets(sce, p=0.5)
sce <- get_mt_malat1(sce)
```
