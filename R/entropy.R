#' Normalise Counts into a Distribution
#'
#' A function that takes frequency count data and normalises it into a
#' probability distribution. Only available internally within SCEnt.
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
#' #Creating Data
#' gene1 <- c(0,0,0,0,1,2,3)
#' gene2 <- c(5,5,3,2,0,0,0)
#' gene3 <- c(2,0,2,1,3,0,1)
#' gene4 <- c(3,3,3,3,3,3,3)
#' gene5 <- c(0,0,0,0,5,0,0)
#' gene_counts <- matrix(c(gene1,gene2,gene3,gene4,gene5), ncol = 5)
#' rownames(gene_counts) <- paste0("cell",1:7)
#' colnames(gene_counts) <- paste0("gene",1:5)
#'
#' #Calculating Homogeneity For Each Gene
#' gene_hom(gene1)
#' gene_hom(gene2)
#' gene_hom(gene3)
#' gene_hom(gene4)
#' gene_hom(gene5)
#'
#' #Calculating Homogeneity For a Matrix
#' gene_hom(gene_counts)
gene_hom <- function(expr, unit = "log2", normalise = TRUE, transpose = FALSE) {
  #Checking if the expression is a vector
  if (is.vector(expr)) {
    #Normalising the data if needed
    if (normalise) {
      expr <- normalise(expr)
    }
    #Calculating the entroy of the expression distribution
    entropy::entropy(expr, unit = unit)
  #Checking if the expression is a matrix
  } else if (is.matrix(expr)) {
    #Transposing matrix if needed
    if (transpose) {
      expr <- t(expr)
    }
    #Normalising the data if needed
    if (normalise) {
      expr <- t(apply(expr, 1, normalise))
    }
    #Calculating the entropy of each expression distribution
    apply(expr, 1, function(x) {entropy::entropy(x, unit = unit)})
  } else {
    #throwing an error if a vector or a matrix was not passed in
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
#' #Creating Data
#' gene1 <- c(0,0,0,0,1,2,3)
#' gene2 <- c(5,5,3,2,0,0,0)
#' gene3 <- c(2,0,2,1,3,0,1)
#' gene4 <- c(3,3,3,3,3,3,3)
#' gene5 <- c(0,0,0,0,5,0,0)
#' gene_counts <- matrix(c(gene1,gene2,gene3,gene4,gene5), ncol = 5)
#' rownames(gene_counts) <- paste0("cell",1:7)
#' colnames(gene_counts) <- paste0("gene",1:5)
#'
#' #Calculating Heterogeneity For Each Gene
#' gene_het(gene1)
#' gene_het(gene2)
#' gene_het(gene3)
#' gene_het(gene4)
#' gene_het(gene5)
#'
#' #Calculating Heterogeneity For a Matrix
#' gene_het(gene_counts)
gene_het <- function(expr, unit = "log2", normalise = TRUE, transpose = FALSE) {
  #Checking if the expression is a matrix
  if (is.vector(expr)) {
    #Normalising the expression if needed
    if (normalise) {
      expr <- normalise(expr)
    }
    #Generating the uniform distribution vector
    num_cell <- length(expr)
    unif_dist <- rep(1/num_cell, num_cell)
    #Finding the KL Divergence of the expression distribution from the
    #uniform distribution
    entropy::KL.plugin(expr, unif_dist, unit = unit)
  #Checking if the data is a matrix
  } else if (is.matrix(expr)) {
    #Transposing the matrix if needed
    if (transpose) {
      expr <- t(expr)
    }
    #Normalising the data if needed
    if (normalise) {
      expr <- t(apply(expr, 1, normalise))
    }
    #Generating the uniform distribution vector
    num_cell <- dim(expr)[2]
    unif_dist <- rep(1/num_cell, num_cell)
    #Finding the KL Divergence of the expression distribution from the
    #uniform distribution for each exoression distribution
    apply(expr, 1, function(x) {entropy::KL.plugin(x, unif_dist, unit = unit)})
  } else {
    #Throwing an error if the input is not a vector or a matrix
    stop()
  }
}
