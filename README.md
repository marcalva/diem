
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
sce <- get_de_genes(sce, low_count=c(0,150), high_count=c(150, Inf), verbose=TRUE)
sce <- set_expression(sce, gene_range=c(200, Inf), simf=2, verbose=TRUE)
sce <- normalize(sce, verbose=TRUE)
sce <- get_bgscore(sce)
sce <- get_pcs(sce, n_pcs=5, verbose=TRUE)
sce <- run_em_pcs(sce, n_pcs=5, verbose=TRUE)
sce <- call_targets(sce, p=p)
sce <- get_mt_malat1(sce)
targets <- call_targets(sce)
```
