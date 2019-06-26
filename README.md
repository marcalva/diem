
# DIEM v0.4

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
sce <- diem(sce)
```
