# JACA: Joint association and classification analysis of multi-view data
The R package 'JACA' implements a joint statistical framework that  performs association and classification analysis for multi-view data, where Multi-view data refers to matched sets of measurements on the same subjects. The corresponding reference is

["Joint association and classification analysis of multi-view data" by Zhang and Gaynanova](https://arxiv.org/abs/1811.08511) (2018+).

## Installation

To install the latest version from Github, use
```s
library(devtools)
devtools::install_github("Pennisetum/JACA")
```
## Usage

```s
library(JACA)

# Example
set.seed(1)
# Generate class indicator matrix Z
n = 100
Z=matrix(c(rep(1, n),rep(0, 2 * n)), byrow = FALSE, nrow = n)
for(i in 1:n){
 Z[i, ] = sample(Z[i, ])
}

# Generate input data X_list
d = 2
X_list = sapply(1:d, function(i) list(matrix(rnorm(n * 20), n, 20)))

# Train JACA model
W = jacaTrain(Z, X_list, lambda = rep(0.05, 2), verbose = FALSE, alpha= 0.5, rho = 0.2)

# Show the number of non-zero rows of each matrix of discriminant vectors
sapply(W, function(x) sum(rowSums(x) != 0))

# Test semi supervised learning
# Set certain class labels and subsets of views as missing 
Z[90:100, ] = rep(NA, 3)
X_list[[1]][1:10, ] = NA
X_list[[2]][11:20, ] = NA
W = jacaTrain(Z, X_list, kmax = 200, eps = 1e-06, lambda = rep(0.05, 2),alpha = 0.5, rho = 0.2, missing = TRUE)

# Show the number of non-zero rows of each matrix of discriminant vectors
sapply(W, function(x) sum(rowSums(x) != 0))
