
# DIEM

[![Build Status](https://travis-ci.com/marcalva/diem.svg?branch=master)](https://travis-ci.com/marcalva/diem)

Remove debris-contaminated droplets from single-cell based data.

## Installation

Currently, we only support installation of the `diem` R package 
through devtools

```R
library(devtools)
devtools::install_github("marcalva/diem")
```

## Usage

Check out the vignette for a thorough tutorial on using `diem`. 

Shown below is a quick workflow for reading 10X data, filtering 
droplets using default parameters, and converting to a 
seurat object. Note you need Seurat installed to run the last step.

```R
library(diem)
library(Seurat)

counts <- read_10x("path/to/10x") # Read 10X data into sparse matrix
sce <- create_SCE(counts) # Create SCE object from counts

# Add MT% and MALAT1%
mt_genes <- grep(pattern="^mt-", x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)
sce <- get_gene_pct(x = sce, genes=mt_genes, name="pct.mt")
malat <- grep(pattern="^malat1$", x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)
sce <- get_gene_pct(x = sce, genes=malat, name="MALAT1")

barcode_rank_plot(sce)

# DIEM steps
sce <- set_debris_test_set(sce)
sce <- filter_genes(sce)
sce <- get_pcs(sce)
sce <- init(sce, k_init = 30)
sce <- get_dist(sce)

# Plot distances
plot_dist(sce)

fltr <- 0.1
sce <- rm_close(sce, fltr = fltr)
sce <- run_em(sce, fltr = fltr)

sce <- call_targets(sce)

seur <- convert_to_seurat(sce)
```

## Version History

February 19, 2020
* version 2.1.0
    * Additional function for extracting Alpha parameters for use with DE
    * Run multiple k_init values at the same time

February 18, 2020
* version 2.0.1
    * Patch fixes some installation issues with tests and docs

February 2, 2020
* version 2.0.0
    * Uses Dirichlet-multinomial
    * Initializes with k-means
    * Removes background centers using a likelihood strategy


