
# Tests functions in init_clust.R

context("Initialization")

#rf <- system.file("data", "mb_small.rda", package="diem")
#load(rf)

expect_equal(ncol(mb_small), 2451)
mb_small <- set_debris_test_set(mb_small)
mb_small <- filter_genes(mb_small)

test_that("initialization works", {
         mb_small <- get_pcs(mb_small, n_var = 200)
         expect_true( all(dim(mb_small@pcs) == c(1752, 30)) )
         expect_equal(length(mb_small@vg), 200)
         expect_equal(nrow(mb_small@norm), 1001)
         expect_false(any(is.na(Matrix::colSums(mb_small@norm))))

         mb_small <- init(mb_small, k_init = 10)
         ic <- mb_small@kruns[["10"]]
         expect_equal( ncol(ic$params$Alpha), ncol(ic$Z) )
         expect_equal( ncol(ic$params$Alpha), length(ic$params$Pi) )

         mb_small <- get_dist(mb_small)
         d <- distances(mb_small)
         expect_equal( length(d), length(ic$params$Pi) )
         expect_equal( sum(d[!is.na(d)] > 1) , 0 )
         expect_equal( sum(d[!is.na(d)] < 0) , 0 )
         expect_equal( sum(d[!is.na(d)] >= 0) , length(ic$params$Pi) - 1 )

         mb_small <- rm_close(mb_small)
         ic <- mb_small@kruns[["10"]]
         expect_equal( ncol(ic$params$Alpha), length(ic$params$Pi) )
         d <- distances(mb_small)
         expect_equal( length(d), length(ic$params$Pi) )
         expect_equal( sum(d[!is.na(d)] > 1) , 0 )
         expect_equal( sum(d[!is.na(d)] < 0) , 0 )
         expect_equal( sum(d[!is.na(d)] >= 0) , length(ic$params$Pi) - 1 )
})

