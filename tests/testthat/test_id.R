
# Test functions in debris_id.R

context("Debris calling")

# library(Matrix)

nr <- 100
nc <- 800

test_thresh <- 25
pp_thresh <- 0.95
min_genes <- 25

db_size <- 50
cell_size <- 500

db_n <- 5000
d1_n <- 100
c1_n <- 100
c2_n <- 100

nc <- db_n + d1_n + c1_n + c2_n
n_clean <- sum(c1_n, c2_n)
n_db <- sum(db_n, d1_n)

#=============================================
# Set random matrix of 100 genes
#=============================================

set.seed(1)
db_prob1 <- rnorm(n = 20, mean = 90)
db_prob2 <- rnorm(n = 80, mean = 5)
db_prob <- c(db_prob1, db_prob2)
db_prob <- db_prob / sum(db_prob)

db_count <- rmultinom(n = db_n, size = db_size, prob = db_prob)

# Contaminated cluster
d1_count <- rmultinom(n = d1_n,  size = cell_size, prob = db_prob)


# Cluster 1
c1_prob1 <- rnorm(n = 10, mean = 10)
c1_prob2 <- rnorm(n = 50, mean = 50)
c1_prob3 <- rnorm(n = 40, mean = 10)
c1_prob <- c(c1_prob1, c1_prob2, c1_prob3)

c1_counts <- rmultinom(n = c1_n, size = cell_size, prob = c1_prob)

# Cluster 2
c2_prob1 <- rnorm(n = 10, mean = 10)
c2_prob2 <- rnorm(n = 70, mean = 10)
c2_prob3 <- rnorm(n = 20, mean = 50)
c2_prob <- c(c2_prob1, c2_prob2, c2_prob3)

ct2_counts <- rmultinom(n = c2_n, size = cell_size, prob = c2_prob)

counts <- do.call(cbind, list(db_count, d1_count, c1_counts, ct2_counts))
rownames(counts) <- paste0("G", as.character(1:nr))
colnames(counts) <- paste0("C", as.character(1:nc))

#=============================================
#=============================================

#=============================================
# Set labels
#=============================================

labels <- c(rep("BG", sum(db_n, d1_n)), 
            rep("CT1", c1_n), 
            rep("CT2", c2_n))
names(labels) <- colnames(counts)


sce <- create_SCE(counts)
sce <- set_debris_test_set(sce, min_counts = 100, min_genes = test_thresh)
sce <- filter_genes(sce, cpm_thresh = 0)
sce <- get_pcs(sce, min_genes = 10, n_var_genes = 30)
sce <- init(sce, k_init = 3)
sce <- run_em(sce, alpha_prior = 0.1, pi_prior = 0.1, verbose = FALSE)
sce <- assign_clusters(sce)
sce <- estimate_dbr_score(sce, thresh_genes = test_thresh)


# Add truth
si <- intersect(rownames(sce@test_data), names(labels))
sce@test_data[si,"CT"] <- labels[si]

ids <- rownames(sce@test_data)

test_that("Clean calling works",{
          sce_c <- call_targets(sce, min_genes = 0)
          dbk <- sce@test_data[,"CT"] == "BG"
          expect_equal(sum(sce_c@test_data[dbk,"Call"] == "Debris"), d1_n)
          expect_equal(sum(sce_c@test_data$Call == "Clean"), n_clean)

          nr1 <- get_removed_ids(sce_c, min_genes = 0)
          expect_equal(length(nr1), d1_n)

          test_cell <- "C5100"
          sce_c <- get_gene_pct(sce_c, "G1", name = "G")
          c1_g <- 100*sce_c@counts["G1",test_cell] / sum(sce_c@counts[,test_cell])
          expect_equal(sce_c@test_data[test_cell,"G"], c1_g)
})

