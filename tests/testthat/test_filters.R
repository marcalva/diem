
# Test functions in filter.R

library(Matrix)

context("Filtering")

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
          sce_ts <- set_debris_test_set(sce, min_counts = 0, min_genes = 0, 
                                        debris_ids = colnames(counts)[1:(ncol(sce) - n_test)])
          expect_equal(length(sce_ts@test_set), n_test)
          expect_equal(length(intersect(sce_ts@test_set, sce_ts@bg_set)), 0)

          sce_ts <- set_debris_test_set(sce, 
                                        min_counts=1500, min_genes=98)
          dd <- sce_ts@droplet_data
          testd <- rownames(dd)[dd$total_counts >= 1500 & dd$n_genes >= 98]
          expect_equal(length(sce_ts@test_set), length(testd))
          expect_equal(length(union(sce_ts@test_set, sce_ts@bg_set)), nc)
          expect_equal(length(intersect(sce_ts@test_set, sce_ts@bg_set)), 0)

          sce_ts <- filter_genes(sce_ts, 0)
          expect_equal(sum(sce_ts@gene_data$exprsd), nr)
          expect_error(filter_genes(sce_ts, 1e7))
})

