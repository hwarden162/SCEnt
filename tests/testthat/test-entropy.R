test_that("gene_hom calculates correct entropies", {
  #Creating test data
  gene1 <- c(0,0,0,0,1,2,3)
  gene2 <- c(5,5,3,2,0,0,0)
  gene3 <- c(2,0,2,1,3,0,1)
  gene4 <- c(3,3,3,3,3,3,3)
  gene5 <- c(0,0,0,0,5,0,0)
  gene_counts <- matrix(c(gene1,gene2,gene3,gene4,gene5), ncol = 5)
  rownames(gene_counts) <- paste0("cell",1:7)
  colnames(gene_counts) <- paste0("gene",1:5)

  expect_equal(
    round(gene_hom(gene1), 6),
    1.459148
  )

  expect_equal(
    round(gene_hom(gene2), 6),
    1.908613
  )

  expect_equal(
    round(gene_hom(gene3), 6),
    2.19716
  )

  expect_equal(
    round(gene_hom(gene4), 6),
    2.807355
  )

  expect_equal(
    round(gene_hom(gene5), 6),
    0
  )

  expected_vector <- c(1.459148, 1.908613, 2.19716, 2.807355, 0)
  names(expected_vector) <- paste0("gene", 1:5)
  expect_equal(
    round(gene_hom(gene_counts, transpose = TRUE), 6),
    expected_vector
  )
})

test_that("gene_het calculates correct entropies", {
  #Creating test data
  gene1 <- c(0,0,0,0,1,2,3)
  gene2 <- c(5,5,3,2,0,0,0)
  gene3 <- c(2,0,2,1,3,0,1)
  gene4 <- c(3,3,3,3,3,3,3)
  gene5 <- c(0,0,0,0,5,0,0)
  gene_counts <- matrix(c(gene1,gene2,gene3,gene4,gene5), ncol = 5)
  rownames(gene_counts) <- paste0("cell",1:7)
  colnames(gene_counts) <- paste0("gene",1:5)

  expect_equal(
    round(gene_het(gene1), 6),
    1.348207
  )

  expect_equal(
    round(gene_het(gene2), 6),
    0.898742
  )

  expect_equal(
    round(gene_het(gene3), 6),
    0.610195
  )

  expect_equal(
    round(gene_het(gene4), 6),
    0
  )

  expect_equal(
    round(gene_het(gene5), 6),
    2.807355
  )

  expected_vector <- c(1.348207, 0.898742, 0.610195, 0, 2.807355)
  names(expected_vector) <- paste0("gene", 1:5)
  expect_equal(
    round(gene_het(gene_counts, transpose = TRUE), 6),
    expected_vector
  )
})
