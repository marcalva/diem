
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

# Variables
path <- "data/s444/raw_gene_bc_matrices/GRCh38.primary_assembly.pli.premrna.genome/"
expected_target <- 5000
min_bg <- 25
max_bg <- 150
n_pcs <- 5
seedn <- 10
min_count <- 200
p <- 0.95

# Run the code
x <- read_10x(path)
sce <- create_SCE(x)
sce <- subset_dropls(sce)
sce <- normalize(sce)
sce <- get_var_genes(sce)
sce <- get_pcs(sce, n_pcs)
sce <- subset_n_dropls(sce, 2*expected_target, slot_name=c("raw", "norm"))
sce <- get_bg_cand(sce, min_bg, max_bg, expected_target)
sce <- run_em_pcs(sce, n_pcs=n_pcs, seedn=seedn)
sce <- call_targets(sce, min_count=min_count, p=p)
```


