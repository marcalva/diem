
# Tests functions in em.R

context("EM")

# rf <- system.file("data", "mb_small.rda", package="diem")
# load(rf)
expect_equal(ncol(mb_small), 2451)

eps <- 1e-8

test_that("EM flags errors when initialized", {
          expect_error(run_em(mb_small))
          mb_small <- set_debris_test_set(mb_small)
          expect_error(run_em(mb_small))
          mb_small <- filter_genes(mb_small)
          expect_error(run_em(mb_small))
})

mb_small <- set_debris_test_set(mb_small)
mb_small <- filter_genes(mb_small)
mb_small <- set_cluster_set(mb_small, cluster_n = 500)
mb_small <- initialize_clusters(mb_small, 
                                nn = 30, n_var = 2000, 
                                min_size = 10, verbose = FALSE)

test_that("EM works when initialized", {
          mb_small <- run_em(mb_small)
          expect_equal(length(mb_small@emo), 6)
          expect_equal(ncol(mb_small@emo$Mu), length(mb_small@ic$assignments))
          expect_equal(ncol(mb_small@emo$Z), length(mb_small@ic$assignments))
          expect_equal(length(mb_small@emo$Mc), length(mb_small@ic$assignments))
          expect_equal(ncol(mb_small@emo$PP), length(mb_small@ic$assignments))
          expect_true(mb_small@emo$converged)
          expect_true(all(rowSums(mb_small@emo$PP) >= (1 - eps) & rowSums(mb_small@emo$PP) <= (1 + eps)))
          expect_true(all(mb_small@droplet_data$CleanProb >= 0 & mb_small@droplet_data$CleanProb <= (1 + eps)))
          expect_true(all(mb_small@droplet_data$ClusterProb >= 0 & mb_small@droplet_data$ClusterProb <= (1 + eps)))
                                })

