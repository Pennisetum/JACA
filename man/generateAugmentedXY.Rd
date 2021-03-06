% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Augmented.R
\name{generateAugmentedXY}
\alias{generateAugmentedXY}
\title{Generate augmented data matrices Y' and X' for JACA based on supplied Z, X_list and alha}
\usage{
generateAugmentedXY(Z, X_list, alpha = 0.5, missing = FALSE)
}
\arguments{
\item{Z}{An N by K class indicator matrix; rows are samples and columns are class indicator vectors with z_k = 1 if observation belongs to class k.}

\item{X_list}{A list of input data matrices; in each sublist, rows are samples and columns are features.}

\item{alpha}{The parameter to control the weight between optimal scoring and CCA part. Default is 0.5.}

\item{missing}{Logical. If False, input data \code{X_list} must be complete and have no missing values.
If True, input data \code{X_list} should contain missing values.}
}
\value{
A list with
\item{bigx}{Augmented matrix X'.}
\item{bigy}{Augmented matrix Y'.}
\item{coef}{A list of length D of scaling coefficients from standardization of X_list when forming X'.}
}
\description{
Generate augmented data matrices Y' and X' for JACA based on supplied Z, X_list and alha
}
\examples{
set.seed(1)
# Generate class indicator matrix Z
n = 10
Z=matrix(c(rep(1, n),rep(0, 2 * n)), byrow = FALSE, nrow = n)
for(i in 1:n){
  Z[i, ] = sample(Z[i, ])
}

# Generate input data X_list
d = 2; p = 5
X_list = sapply(1:d, function(i) list(matrix(rnorm(n * p), n, p)))

# Generate augmented X' and Y'
out = generateAugmentedXY(Z, X_list)
out$bigx
out$bigy

}
