
# Test functions in filter.R

context("Filtering")

library(Matrix)

nr <- 100
nc <- 100

set.seed(1)
counts <- sample(0:30, replace = TRUE, size = 1e4)
counts <- matrix(counts, nrow = nr, ncol = nc)
rownames(counts) <- paste0("G", as.character(1:nr))
colnames(counts) <- paste0("C", as.character(1:nc))

sce <- create_SCE(counts)

test_that("Filtering works", {

          n_test <- 10
         sce_ts <- set_debris_test_set(sce, top_n = 1e4,
                                       min_counts = 0, min_genes = 0, 
                                       fix_debris = colnames(counts)[1:(ncol(sce) - n_test)])
         expect_equal(length(sce_ts@test_set), n_test)
         expect_equal(length(intersect(sce_ts@test_set, sce_ts@bg_set)), 0)

         sce_ts <- set_debris_test_set(sce, top_n=1e4, 
                                       min_counts=1500, min_genes=98)
         dd <- sce_ts@droplet_data
         testd <- rownames(dd)[dd$total_counts >= 1500 & dd$n_genes >= 98]
         expect_equal(length(sce_ts@test_set), length(testd))
         expect_equal(length(union(sce_ts@test_set, sce_ts@bg_set)), nc)
         expect_equal(length(intersect(sce_ts@test_set, sce_ts@bg_set)), 0)

         sce_ts <- set_cluster_set(sce_ts, cluster_n = 2)
         dd <- sce_ts@droplet_data[sce_ts@test_set,,drop=FALSE]
         mg <- dd$n_genes[order(dd$n_genes, decreasing=TRUE)][2]
         top2 <- rownames(dd)[dd$n_genes >= mg]
         expect_equal(sce_ts@cluster_set, top2)
         
         sce_ts <- set_cluster_set(sce_ts, cluster_n = 1000)
         expect_equal(sort(sce_ts@cluster_set), sort(sce_ts@test_set))

         sce_ts <- filter_genes(sce_ts, 0)
         expect_equal(sum(sce_ts@gene_data$exprsd), nr)
         expect_error(filter_genes(sce_ts, 1e7))
})

