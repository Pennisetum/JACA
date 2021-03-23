#' Generate augmented data matrices Y' and X' for JACA based on supplied Z, X_list and alha
#'
#' @inheritParams jacaTrain
#'
#' @return A list with
#' \item{bigx}{Augmented matrix X'.}
#' \item{bigy}{Augmented matrix Y'.}
#' \item{coef}{A list of length D of scaling coefficients from standardization of X_list when forming X'.}
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' # Generate class indicator matrix Z
#' n = 10
#' Z=matrix(c(rep(1, n),rep(0, 2 * n)), byrow = FALSE, nrow = n)
#' for(i in 1:n){
#'   Z[i, ] = sample(Z[i, ])
#' }
#'
#' # Generate input data X_list
#' d = 2; p = 5
#' X_list = sapply(1:d, function(i) list(matrix(rnorm(n * p), n, p)))
#'
#' # Generate augmented X' and Y'
#' out = generateAugmentedXY(Z, X_list)
#' out$bigx
#' out$bigy
#'
generateAugmentedXY <- function(Z, X_list, alpha = 0.5, missing = FALSE){
  # Determine the number of views
  D = length(X_list)
  # Determine the total number of samples
  n = nrow(Z)

  # If there are missing views
  if (missing) {
    # Find out sample IDs with not NA class
    Z_nonNA_id = which(apply(Z, 1, function(x) sum(is.na(x)) == 0))
    # Find out sample IDs with not NA view d (for each view)
    X_nonNA_id = lapply(X_list, function(sublist) {
      which(apply(sublist, 1, function(x) sum(is.na(x)) == 0))
    })

    # compute the union ID (Y and X)
    # These are all samples that have class/1 view available
    inter_ID_YX = lapply(X_nonNA_id, function(x) sort(intersect(Z_nonNA_id, x)))
    # These are all samples that have at least 2 view pairs available
    inter_ID_XX = lapply(combn(1:D, 2, simplify = F), function(x) sort(intersect(X_nonNA_id[[x[1]]],
                                                                                 X_nonNA_id[[x[2]]])))

    # center and scale each view X_d
    # scale function can work with NA measurements
    coef = lapply(X_list, function(x) apply(x, 2, sd, na.rm = TRUE) * sqrt((nrow(na.omit(x)) - 1)/nrow(na.omit(x))))
    centeredX = lapply(1:D, function(i) scale(X_list[[i]], scale = coef[[i]], center = T))

    # generate augmented data matrices X' and Y'
    # generate augmented X with given alpha
    bigx = xprimeMissing(centeredX, inter_ID_YX, inter_ID_XX, alpha = alpha)
    # generate augmented Y
    transy = matrix(NA, ncol = ncol(Z) - 1, nrow = nrow(Z))
    transy[Z_nonNA_id, ] = transformYCpp(Z[Z_nonNA_id, ])
    bigy = matrix(0, ncol = ncol(transy), nrow = nrow(bigx))
    temp_y = lapply(inter_ID_YX, function(obs_idx) transy[obs_idx, , drop = F]/sqrt(n * D))
    bigy[1:(sum(sapply(inter_ID_YX, length))), ] = sqrt(alpha) * do.call(rbind, temp_y)

  } else {
    # There should be no missing views or missing classes
    if (any(anyNA(Z), sapply(X_list, anyNA))) {
      stop("Datasets contain missing values! Please set missing = TRUE")
    }

    # center and scale each view X_d
    coef = lapply(X_list, function(x) apply(x, 2, sd, na.rm = TRUE) * sqrt((n - 1)/n))
    centeredX = lapply(1:D, function(i) scale(X_list[[i]], scale = coef[[i]], center = T))

    # generate augmented data matrices X' and Y'
    # generate augmented C with given alpha
    bigx = xPrime(centeredX, alpha = alpha)
    # generate augmented Y
    transy = transformYCpp(Z)
    bigy = matrix(0, ncol = ncol(transy), nrow = nrow(bigx))
    bigy[1:(n * D), ] = sqrt(alpha) * do.call(rbind, replicate(D, transy, simplify = FALSE))/sqrt(n * D)
  }
  return(list(bigx = bigx, bigy = bigy, coef = coef))
}
