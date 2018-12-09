#' JACA: Joint Association and Classification Analysis
#'
#' This package implements "Joint association and classification analysis of multi-view data".
#' Multi-view data refers to matched sets of measurements on the same subjects,
#' and JACA is a joint statistical framework that performs association and classification analysis simultaneously.
#'
#' @details  This package implements three main methods that can:
#' \itemize{
#'  \item Fit the JACA model.
#'  \item Use cross-validation function to automatically select parameters.
#'  \item Classify observations using the supplied model.
#' }
#'
#' @docType package
#' @name JACA-package
#'
#' @author Yunfeng Zhang and Irina Gaynanova
#'
#' @references Zhang, Y., & Gaynanova, I. (2018). Joint association and classification analysis of multi-view data. \emph{arXiv preprint arXiv:1811.08511}.
#'
#' @importFrom stats predict sd var
#' @importFrom utils combn
#' @importFrom Rcpp  evalCpp
#'
#' @useDynLib JACA
"_PACKAGE"
#'

