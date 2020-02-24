
# Tests functions in init_clust.R

context("Initialization")

#rf <- system.file("data", "mb_small.rda", package="diem")
#load(rf)

mb_small <- set_debris_test_set(mb_small, verbose = FALSE)
mb_small <- filter_genes(mb_small, verbose = FALSE)

test_that("initialization works", {
         mb_small <- get_pcs(mb_small, n_var = 200, n_pcs = 30)
         expect_true( ncol(mb_small@pcs) == 30 ) 
         expect_equal(length(mb_small@vg), 200)

         mb_small <- init(mb_small, k_init = 15, verbose = FALSE)
         ic <- mb_small@kruns[["15"]]
         expect_equal( ncol(ic$params$Alpha), ncol(ic$llk) )
         expect_equal( ncol(ic$params$Alpha), length(ic$params$Pi) )

         mb_small <- get_dist(mb_small, verbose = FALSE)
         d <- distances(mb_small, k_init = "15")
         expect_equal( length(d), length(ic$params$Pi) )
         expect_equal( sum(d[!is.na(d)] > 1) , 0 )
         expect_equal( sum(d[!is.na(d)] < 0) , 0 )
         expect_equal( sum(d[!is.na(d)] >= 0) , length(ic$params$Pi) - 1 )

         mb_small <- rm_close(mb_small, verbose = FALSE)
         ic <- mb_small@kruns[["15"]]
         expect_equal( ncol(ic$params$Alpha), length(ic$params$Pi) )
         d <- distances(mb_small, k_init = "15")
         expect_equal( length(d), length(ic$params$Pi) )
         expect_equal( sum(d[!is.na(d)] > 1) , 0 )
         expect_equal( sum(d[!is.na(d)] < 0) , 0 )
         expect_equal( sum(d[!is.na(d)] >= 0) , length(ic$params$Pi) - 1 )
})

