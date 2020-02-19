
# Tests functions in init_clust.R

context("Initialization")

#rf <- system.file("data", "mb_small.rda", package="diem")
#load(rf)

expect_equal(ncol(mb_small), 2451)
mb_small <- set_debris_test_set(mb_small)
mb_small <- filter_genes(mb_small)
mb_small <- init(mb_small)

test_that("initialization works", {
         expect_true( all(dim(mb_small@pcs) == c(length(mb_small@test_set), 30)) )
         expect_equal(nrow(mb_small@norm), nrow(mb_small@counts))
         expect_false(any(is.na(Matrix::colSums(mb_small@norm))))

         mb_small <- init(mb_small, k_init = 10)
         ic <- mb_small@emo
         expect_equal( ncol(ic$params$Alpha), ncol(ic$Z) )
})

