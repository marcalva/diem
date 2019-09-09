
# Test functions in debris_id.R

library(Matrix)

nr <- 100
nc <- 100

set.seed(10)
counts <- sample(0:30, replace=TRUE, size=1e4)
counts <- matrix(counts, nrow=nr, ncol=nc)
rownames(counts) <- paste0("G", as.character(1:nr))
colnames(counts) <- paste0("C", as.character(1:nc))

sce <- create_SCE(counts)
sce@droplet_data$CleanProb <- 0
sce@droplet_data[sce@droplet_data$n_genes >= 98, "CleanProb"] <- 1

mg <- 100
n_all <- table(sce@droplet_data$n_genes >= mg)[["TRUE"]]
names_pass <- rownames(sce@droplet_data)[ sce@droplet_data$n_genes >= mg ]

test_that("Clean calling works",{
          sce_c <- call_targets(sce, pp_thresh=0.95, min_genes=mg)
          expect_equal(sum(sce_c@droplet_data$Call == "Clean"), n_all)
          expect_equal(sort(get_clean_ids(sce_c)), sort(names_pass))
          nr1 <- get_removed_ids(sce_c, min_genes=99)
          nr2 <- rownames(sce@droplet_data)[sce@droplet_data$n_genes == 99]
          expect_equal(sort(nr1), sort(nr2))
          sce_c <- get_gene_pct(sce_c, "G1", name="G")
          expect_equal(sce_c@droplet_data["C1","G"], 100*15/1329)

          sce_c <- call_targets(sce, pp_thresh=1.1, min_genes=0)
          expect_equal(sum(sce_c@droplet_data$Call == "Clean"), 0)
          expect_equal(length(get_clean_ids(sce_c)), 0)
          rmd <- get_removed_ids(sce_c, min_genes = 0)
          expect_equal(length(rmd), nr)
          expect_error(get_gene_pct(sce_c, gene="doesntexist", "NAME"))
}

