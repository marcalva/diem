## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(diem)

# In the diem directory
counts <- read_10x("../tests/testdata/")
dim(counts)
class(counts)

## ------------------------------------------------------------------------
mb_small <- create_SCE(counts, name="MouseBrain")
dim(mb_small)
class(mb_small)

## ------------------------------------------------------------------------
drop_data <- droplet_data(mb_small)
head(drop_data)
summary(drop_data)

## ------------------------------------------------------------------------
mt_genes <- grep(pattern = "^mt-", x = rownames(mb_small@gene_data), 
                 ignore.case = TRUE, value = TRUE)
mb_small <- get_gene_pct(x = mb_small, genes = mt_genes, name = "pct.mt")
genes <- grep(pattern = "^malat1$", x = rownames(mb_small@gene_data), 
              ignore.case = TRUE, value = TRUE)
mb_small <- get_gene_pct(x = mb_small, genes = genes, name = "MALAT1")
drop_data <- droplet_data(mb_small)
summary(drop_data)

## ---- fig.width=14, fig.height=14----------------------------------------
datf <- drop_data
par(mfrow=c(2,2))
plot(datf$total_counts, datf$n_genes, 
     xlab="Total Counts", ylab="Number Genes")
plot(datf$n_genes, datf$pct.mt,
     xlab="Number Genes", ylab="MT%")
plot(datf$n_genes, datf$MALAT1, 
     xlab="Number Genes", ylab="MALAT1%")
plot(datf$pct.mt, datf$MALAT1,
     xlab="MT%", ylab="MALAT1%")

## ------------------------------------------------------------------------
barcode_rank_plot(mb_small, title = "MouseBrain")

## ------------------------------------------------------------------------
mb_small <- set_debris_test_set(mb_small, 
                                top_n = 2000, 
                                min_counts = 100, 
                                min_genes = 100)
length(mb_small@test_set)
length(mb_small@bg_set)

## ------------------------------------------------------------------------
mb_small <- filter_genes(mb_small, cpm_thresh = 10)
genes <- gene_data(mb_small)
summary(genes)

## ------------------------------------------------------------------------
mb_small <- set_cluster_set(mb_small, 
                            cluster_n = 200, 
                            order_by = "gene")
length(mb_small@cluster_set)

## ------------------------------------------------------------------------
mb_small <- initialize_clusters(mb_small, 
                                use_var = TRUE, 
                                n_var = 200, 
                                verbose = TRUE)

## ------------------------------------------------------------------------
mb_small <- run_em(mb_small, verbose = TRUE)
drop_data <- droplet_data(mb_small)
summary(drop_data)

## ------------------------------------------------------------------------
mb_small <- call_targets(mb_small, 
                         pp_thresh = 0.95, 
                         min_genes = 100)
drop_data <- droplet_data(mb_small)
summary(drop_data)

## ------------------------------------------------------------------------
clean <- get_clean_ids(mb_small)
removed <- setdiff(mb_small@test_set, clean)
length(clean)
length(removed)
length(removed) / (length(removed) + length(clean))

## ---- fig.show='hold'----------------------------------------------------
plot_umi_gene_call(mb_small, alpha = 1)

## ------------------------------------------------------------------------
library(ggplot2)
drop_data <- droplet_data(mb_small)
ggplot(drop_data, aes(x = n_genes, y = pct.mt, color = CleanProb)) + 
geom_point() + 
theme_minimal() + 
xlab("Number Genes") + ylab("MT%")

## ------------------------------------------------------------------------
counts <- raw_counts(mb_small)
counts_clean <- counts[,clean]
dim(counts_clean)

## ------------------------------------------------------------------------
seur <- convert_to_seurat(mb_small, 
                          targets = TRUE, 
                          meta = TRUE, 
                          min.features = 100)
dim(seur)
head(seur@meta.data)

