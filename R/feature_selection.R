#' Feature Selection by Gene Heterogeneity
#'
#' @param expr A matrix of gene expression data. Cells should be represented
#' as rows and genes should be represented as columns.
#' @param bit_threshold The threshold for the amount of bits of information a
#' gene must add to be selected as a feature. Only one threshold can be used
#' at a time.
#' @param count_threshold A number represented how many of the most
#' heterogeneous cells should be selected. Only one threshold can be used at
#' a time.
#' @param perc_threshold The percentile of the hetergeneity distribution above
#' which a gene should be to be selected as a feature.
#' @param unit The units to be used when calculating entropy.
#' @param normalise A logical value representing whether the gene counts
#' should be normalised into a probability distribution.
#' @param transpose  A logical value representing whether the matrix should
#' be transposed before having any operations computed on it.
#'
#' @return A matrix of gene expression values where genes with low
#' heterogeneity have been removed.
#' @export
#'
#' @examples
scent_select <- function(expr, bit_threshold = NULL, count_threshold = NULL, perc_threshold = NULL, unit = "log2", normalise = TRUE, transpose = FALSE) {
  if (transpose) {
    expr <- t(expr)
  }

  if (is.null(bit_threshold) + is.null(count_threshold) + is.null(perc_threshold) != 2) {
    stop()
  }

  het_vals <- gene_het(expr, unit, normalise, transpose = TRUE)

  if (!is.null(bit_threshold)) {
    indices <- het_vals > bit_threshold
    return(as.matrix(expr[,indices]))
  }

  if (!is.null(count_threshold)) {
    indices <- sort(het_vals, decreasing = TRUE, index.return = TRUE)$ix
    if (length(indices) < count_threshold) {
      return(expr)
    } else {
      indices <- sort(indices[1:count_threshold])
      return(as.matrix(expr[,indices]))
    }
  }

  if (!is.null(perc_threshold)) {
    indices <- het_vals >= stats::quantile(het_vals, perc_threshold)
    return(as.matrix(expr[,indices]))
  }
}

#' A Tidy Wrapper for Feature Selection by Heterogeneity
#'
#' @param expr A tibble of gene expression data. Cells should be represented
#' as rows and genes should be represented as columns.
#' @param bit_threshold The threshold for the amount of bits of information a
#' gene must add to be selected as a feature. Only one threshold can be used
#' at a time.
#' @param count_threshold A number represented how many of the most
#' heterogeneous cells should be selected. Only one threshold can be used at
#' a time.
#' @param perc_threshold The percentile of the hetergeneity distribution above
#' which a gene should be to be selected as a feature.
#' @param unit The units to be used when calculating entropy.
#' @param normalise A logical value representing whether the gene counts
#' should be normalised into a probability distribution.
#' @param transpose  A logical value representing whether the matrix should
#' be transposed before having any operations computed on it.
#'
#' @return A tibble of gene expression values where genes with low
#' heterogeneity have been removed.
#' @export
#'
#' @examples
scent_select_tidy <- function(expr, bit_threshold = NULL, count_threshold = NULL, perc_threshold = NULL, unit = "log2", normalise = TRUE, transpose = FALSE) {
  expr <- as.matrix(expr)
  reduced_expr <- scent_select(expr, bit_threshold, count_threshold, perc_threshold, unit, normalise, transpose)
  tibble::as_tibble(reduced_expr)
}
