test_that("scent_select() selects the correct columns", {
  # Creating test data
  gene1 <- c(0, 0, 0, 0, 1, 2, 3)
  gene2 <- c(5, 5, 3, 2, 0, 0, 0)
  gene3 <- c(2, 0, 2, 1, 3, 0, 1)
  gene4 <- c(3, 3, 3, 3, 3, 3, 3)
  gene5 <- c(0, 0, 0, 0, 5, 0, 0)
  gene_counts <- matrix(c(gene1, gene2, gene3, gene4, gene5), ncol = 5)
  rownames(gene_counts) <- paste0("cell", 1:7)
  colnames(gene_counts) <- paste0("gene", 1:5)

  expect_equal(
    scent_select(gene_counts, bit_threshold = 0.85),
    gene_counts[, c(1, 2, 5)]
  )

  expect_equal(
    scent_select(gene_counts, count_threshold = 2),
    gene_counts[, c(1, 5)]
  )

  expect_equal(
    scent_select(gene_counts, perc_threshold = 0.25),
    gene_counts[, c(1, 2, 3, 5)]
  )
})

test_that("scent_select_tidy() selects the correct columns", {
  # Creating test data
  library(tibble)
  gene1 <- c(0, 0, 0, 0, 1, 2, 3)
  gene2 <- c(5, 5, 3, 2, 0, 0, 0)
  gene3 <- c(2, 0, 2, 1, 3, 0, 1)
  gene4 <- c(3, 3, 3, 3, 3, 3, 3)
  gene5 <- c(0, 0, 0, 0, 5, 0, 0)
  gene_counts <- matrix(c(gene1, gene2, gene3, gene4, gene5), ncol = 5)
  rownames(gene_counts) <- paste0("cell", 1:7)
  colnames(gene_counts) <- paste0("gene", 1:5)
  gene_counts <- as_tibble(gene_counts)

  expect_equal(
    scent_select_tidy(gene_counts, bit_threshold = 0.85),
    gene_counts[, c(1, 2, 5)]
  )

  expect_equal(
    scent_select_tidy(gene_counts, count_threshold = 2),
    gene_counts[, c(1, 5)]
  )

  expect_equal(
    scent_select_tidy(gene_counts, perc_threshold = 0.25),
    gene_counts[, c(1, 2, 3, 5)]
  )
})

test_that("gene_counts sums correctly", {
  # Creating test data
  gene1 <- c(0, 0, 0, 0, 1, 2, 3)
  gene2 <- c(5, 5, 3, 2, 0, 0, 0)
  gene3 <- c(2, 0, 2, 1, 3, 0, 1)
  gene4 <- c(3, 3, 3, 3, 3, 3, 3)
  gene5 <- c(0, 0, 0, 0, 5, 0, 0)
  gene_counts <- matrix(c(gene1, gene2, gene3, gene4, gene5), ncol = 5)
  rownames(gene_counts) <- paste0("cell", 1:7)
  colnames(gene_counts) <- paste0("gene", 1:5)

  expect_equal(
    gene_counts(gene_counts),
    apply(gene_counts, 2, sum)
  )
})

test_that("rm_low_counts filters correctly", {
  # Creating test data
  gene1 <- c(0, 0, 0, 0, 1, 2, 3)
  gene2 <- c(5, 5, 3, 2, 0, 0, 0)
  gene3 <- c(2, 0, 2, 1, 3, 0, 1)
  gene4 <- c(3, 3, 3, 3, 3, 3, 3)
  gene5 <- c(0, 0, 0, 0, 5, 0, 0)
  gene_counts <- matrix(c(gene1, gene2, gene3, gene4, gene5), ncol = 5)
  rownames(gene_counts) <- paste0("cell", 1:7)
  colnames(gene_counts) <- paste0("gene", 1:5)

  expect_equal(
    rm_low_counts(gene_counts, count_threshold = 7),
    gene_counts[,2:4]
  )

  expect_equal(
    rm_low_counts(gene_counts, perc_threshold = 0.1),
    gene_counts[,1:4]
  )
})




test_that("rm_low_counts filters correctly", {
  # Creating test data
  library(tibble)
  gene1 <- c(0, 0, 0, 0, 1, 2, 3)
  gene2 <- c(5, 5, 3, 2, 0, 0, 0)
  gene3 <- c(2, 0, 2, 1, 3, 0, 1)
  gene4 <- c(3, 3, 3, 3, 3, 3, 3)
  gene5 <- c(0, 0, 0, 0, 5, 0, 0)
  gene_counts <- matrix(c(gene1, gene2, gene3, gene4, gene5), ncol = 5)
  rownames(gene_counts) <- paste0("cell", 1:7)
  colnames(gene_counts) <- paste0("gene", 1:5)
  gene_counts <- as_tibble(gene_counts)

  expect_equal(
    rm_low_counts_tidy(gene_counts, count_threshold = 7),
    gene_counts[,2:4]
  )

  expect_equal(
    rm_low_counts_tidy(gene_counts, perc_threshold = 0.1),
    gene_counts[,1:4]
  )
})
