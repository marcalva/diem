
# Tests functions in em.R

context("EM")

eps <- 1e-8

test_that("EM flags errors when initialized", {
          expect_error(run_em(mb_small))
          mb_small <- set_debris_test_set(mb_small, verbose=FALSE)
          expect_error(run_em(mb_small))
          mb_small <- filter_genes(mb_small, verbose=FALSE)
          expect_error(run_em(mb_small))
})

mb_small <- set_debris_test_set(mb_small, min_counts = 5, verbose = FALSE)
mb_small <- filter_genes(mb_small, verbose = FALSE)
mb_small <- get_pcs(mb_small)
mb_small <- init(mb_small, k_init = 15, verbose = FALSE)


test_that("EM works when initialized", {
          mb_small <- run_em(mb_small, k_init = "15", verbose=FALSE)
          mb_small <- call_targets(mb_small, min_genes = 5, k_init = "15", verbose=FALSE)
          emo <- mb_small@kruns[["15"]]
          expect_equal(length(emo), 5)
          expect_equal(ncol(emo$params$Alpha), length(emo$params$Pi))
          expect_true(emo$converged)
                                })

