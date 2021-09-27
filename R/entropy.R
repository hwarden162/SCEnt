#' Normalise Counts into a Distribution
#'
#' A function that takes frequency count data and normalises it into a
#' probability distribution.
#'
#' @param dist A vector of a frequency distribution.
#'
#' @return A vector of a probability distribution relative to the
#' frequencies.
normalise <- function(dist) {
  dist/sum(dist)
}

#' Find the Homogeneity of a Gene Within a Population
#'
#' @param expr A vector or matrix of gene expressions. For the matrix, genes
#' should be represented as rows and cells as columns.
#' @param unit The units to be parsed to the entropy function.
#' @param normalise A logical value representing whether the gene frequencies
#' should be normalised into a distribution.
#' @param transpose A legical value representing whether the matrix should be
#' transposed before any calculations are performed.
#'
#' @return A vector of the information contained in the distribution of each
#' gene. The higher this is, the more homogeneous the gene is within the cell
#' population.
#' @export
#'
#' @examples
gene_hom <- function(expr, unit = "log2", normalise = TRUE, transpose = FALSE) {
  if (is.vector(expr)) {
    if (normalise) {
      expr <- normalise(expr)
    }
    entropy::entropy(expr, unit = unit)
  } else if (is.matrix(expr)) {
    if (normalise) {
      expr <- t(apply(expr, 1, normalise))
    }
    apply(expr, 1, function(x) {entropy::entropy(x, unit = unit)})
  } else {
    stop()
  }
}

#' Find the Heterogeneity of a Gene Within a Population
#'
#' @param expr A vector or matrix of gene expressions. For the matrix, genes
#' should be represented as rows and cells as columns.
#' @param unit The units to be parsed to the entropy function.
#' @param normalise A logical value representing whether the gene frequencies
#' should be normalised into a distribution.
#' @param transpose A logical value representing whether the matrix should be
#' transposed before any calculations are performed.
#'
#' @return A vector of the information gained from the gene distribution
#' compared to the uniform distribution. The higher the value more
#' heterogeneous the cell is within the population.
#' @export
#'
#' @examples
gene_het <- function(expr, unit = "log2", normalise = TRUE, transpose = FALSE) {
  if (is.vector(expr)) {
    if (normalise) {
      expr <- normalise(expr)
    }
    num_cell <- length(expr)
    unif_dist <- rep(1/num_cell, num_cell)
    entropy::KL.plugin(expr, unif_dist, unit = unit)
  } else if (is.matrix(expr)) {
    if (transpose) {
      expr <- t(expr)
    }
    if (normalise) {
      expr <- t(apply(expr, 1, normalise))
    }
    num_cell <- dim(expr)[2]
    unif_dist <- rep(1/num_cell, num_cell)
    apply(expr, 1, function(x) {entropy::KL.plugin(x, unif_dist, unit = unit)})
  } else {
    stop()
  }
}
