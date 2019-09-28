
# Tests functions in init_clust.R

context("Initialization")

#rf <- system.file("data", "mb_small.rda", package="diem")
#load(rf)

expect_equal(ncol(mb_small), 2451)
mb_small <- set_debris_test_set(mb_small)
mb_small <- filter_genes(mb_small)

test_that("initialization works", {
         mb_small <- set_cluster_set(mb_small, cluster_n = 200, order_by="count")
         expect_equal(length(mb_small@cluster_set), 200)
         mb_small <- set_cluster_set(mb_small, cluster_n = 200)
         expect_equal(length(mb_small@cluster_set), 201)
         expect_error(set_cluster_set(mb_small, cluster_n = 0))
         mb_small <- set_cluster_set(mb_small, cluster_n = Inf)
         expect_equal(sort(mb_small@cluster_set), sort(mb_small@test_set))

         mb_small <- get_var_genes(mb_small, n_genes = 200)
         expect_equal(length(mb_small@vg), 200)

         mb_small <- normalize_data(mb_small, use_var=TRUE)
         expect_equal(nrow(mb_small@norm), 200)
         expect_false(any(is.na(Matrix::colSums(mb_small@norm))))

         mb_small <- set_cluster_set(mb_small, cluster_n = 500)
         mb_small <- initialize_clusters(mb_small, 
                                         nn = 30, 
                                         n_var = 2000, 
                                         min_size = 10, 
                                         verbose = FALSE)
         expect_lte(length(mb_small@ic$clusters), ncol(mb_small))
         expect_equal(nlevels(mb_small@ic$assignments), 2)
         expect_equal(length(mb_small@ic$assignments), length(table(mb_small@ic$clusters)))
         expect_error(initialize_clusters(mb_small, cluster_n = 0, verbose = FALSE))
         expect_error(initialize_clusters(mb_small, min_size = ncol(mb_small)+1, verbose = FALSE))
})

