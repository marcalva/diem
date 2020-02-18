
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

mb_small <- set_debris_test_set(mb_small, min_counts = 100)
mb_small <- filter_genes(mb_small)
mb_small <- get_pcs(mb_small)
mb_small <- init(mb_small, k_init = 15)
mb_small <- get_dist(mb_small)
mb_small <- rm_close(mb_small, fltr = 0.5)


test_that("EM works when initialized", {
          mb_small <- run_em(mb_small, fltr = 0.5)
          mb_small <- call_targets(mb_small, min_genes = 100)
          emo <- mb_small@kruns[["15"]]
          expect_equal(length(emo), 6)
          expect_equal(ncol(emo$params$Alpha), length(emo$params$Pi))
          expect_equal(ncol(emo$Z), length(emo$params$Pi))
          expect_equal(ncol(emo$llk), ncol(emo$Z))
          expect_equal(length(emo$Dist), length(emo$params$Pi))
          expect_true(emo$converged)
          expect_true(all(rowSums(emo$Z) >= (1 - eps) & rowSums(emo$Z) <= (1 + eps)))
          expect_equal(emo$removed, 0)
                                })

