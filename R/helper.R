# a functioin to generate the augmented data matrix X'
xPrime <- function(centeredX, alpha) {
  ### a functioin to generate the augmented data matrix X' number of measurements
  D = length(centeredX)
  # sample size
  n = nrow(centeredX[[1]])

  if (D == 1) {
    return(sqrt(alpha) * centeredX[[1]]/sqrt(n * D))
  } else {
    # generate diagonal part of X'
    diagX = do.call(Matrix::bdiag, centeredX)/sqrt(n * D) * sqrt(alpha)
    lower1 = matrix(0, ncol = D, nrow = (D - 1) * D/2)
    rou_count = 1
    for (i in 1:(D - 1)) {
      for (j in i:(D - 1)) {
        lower1[rou_count, i] = 1
        lower1[rou_count, j + 1] = -1
        rou_count = rou_count + 1
      }
    }

    # generate the lower part of X'
    lowerX_temp = list()
    for (i in 1:nrow(lower1)) {
      temp = lapply(1:D, function(idx) lower1[i, idx] * centeredX[[idx]])
      lowerX_temp[[i]] = do.call(cbind, temp)
    }

    lowerX = do.call(rbind, lowerX_temp)/sqrt(n * D * (D - 1)) * sqrt(1 - alpha)
    xPrime = matrix(0, ncol = ncol(lowerX), nrow = nrow(lowerX) + nrow(diagX))
    xPrime[1:nrow(diagX), ] = as.matrix(diagX)
    xPrime[-(1:nrow(diagX)), ] = lowerX
    return(xPrime)
  }
}

# a functioin to generate the augmented data matrix X' when input data has missing value
xprimeMissing <- function(centeredX, inter_ID_YX, inter_ID_XX, alpha) {
  ### a functioin to generate the augmented data matrix X' number of measurements
  D = length(centeredX)
  n = nrow(centeredX[[1]])
  if (D == 1) {
    return(sqrt(alpha) * centeredX[[1]][inter_ID_YX[[1]], ]/sqrt(n * D))
  } else {
    # generate diagonal part of X'
    centeredX_diag = lapply(1:D, function(idx) centeredX[[idx]][inter_ID_YX[[idx]], , drop = F]/sqrt(n *
                                                                                                       D))
    diagX = do.call(Matrix::bdiag, centeredX_diag) * sqrt(alpha)
    lower1 = matrix(0, ncol = D, nrow = (D - 1) * D/2)
    rou_count = 1
    for (i in 1:(D - 1)) {
      for (j in i:(D - 1)) {
        lower1[rou_count, i] = 1
        lower1[rou_count, j + 1] = -1
        rou_count = rou_count + 1
      }
    }

    # generate the lower part of X'
    lowerX_temp = list()
    for (i in 1:nrow(lower1)) {
      temp = lapply(1:D, function(idx) lower1[i, idx] * centeredX[[idx]][inter_ID_XX[[i]], , drop = F]/sqrt(n *
                                                                                                              D * (D - 1)))

      lowerX_temp[[i]] = do.call(cbind, temp)
    }

    lowerX = do.call(rbind, lowerX_temp) * sqrt(1 - alpha)
    # Xprime = as.matrix(rbind(diagX, lowerX))
    Xprime = matrix(0, ncol = ncol(lowerX), nrow = nrow(lowerX) + nrow(diagX))
    Xprime[1:nrow(diagX), ] = as.matrix(diagX)
    Xprime[-(1:nrow(diagX)), ] = lowerX
    Xprime[is.na(Xprime)] <- 0
    return(Xprime)
  }
}

# define the evaluation function for cross validation
objectiveJACA = function(W_list, test_z, train_x, test_x, D, theta, alpha) {
  ## W_list: the results we get from JACA. D: number of data sets.
  Ytilde = test_z %*% theta
  n = nrow(test_z)

  if (D != 1) {
    # Xbar: the mean values of train_x
    Xbar = lapply(train_x, function(mat) colMeans(mat, na.rm = T))
    centeredX = lapply(1:D, function(num) test_x[[num]] - matrix(rep(Xbar[[num]], nrow(test_x[[num]])),
                                                                 byrow = T, nrow = nrow(test_x[[num]])))

    # Use correlation as cv criterion the part1 of the objective function
    part1 = sum(sapply(1:D, function(num) -sum(cv_cor(Ytilde, centeredX[[num]] %*% W_list[[num]]))))
    # the part2 of the objective function
    part2 = sapply(combn(1:D, 2, simplify = F), function(com) -sum(cv_cor(centeredX[[com[1]]] %*%
                                                                            W_list[[com[1]]], centeredX[[com[2]]] %*% W_list[[com[2]]])))
    part2 = sum(part2)

    return(alpha * part1 + (1 - alpha) * part2/(D - 1))
  } else {
    # Xbar: the mean values of train_x
    Xbar = lapply(train_x, function(mat) colMeans(mat))
    centeredX = lapply(1:D, function(num) test_x[[num]] - matrix(rep(Xbar[[num]], nrow(test_x[[num]])),
                                                                 byrow = T, nrow = nrow(test_x[[num]])))

    # the part1 of the objective function
    part1 = sum(sapply(1:D, function(num) -sum(cv_cor(Ytilde, centeredX[[num]] %*% W_list[[num]]))))
    # the part2 of the objective function
    part2 = 0

    return(alpha * part1 + (1 - alpha) * part2/(D - 1))
  }
}

# Given the the class indicator matrix Z, generate \tilde Y
transformY <- function(Z) {
  # calculate the size of each class and the cumulative of sizes of classes (for later use)
  nclass = colSums(Z, na.rm = T)
  cumnclass = cumsum(nclass)

  # number of observations and number of class minus 1
  n = nrow(Z)
  r = ncol(Z) - 1

  # calculate the Theta matrix
  theta = matrix(0, nrow = r + 1, ncol = r)
  for (l in 1:r) {
    theta[, l] = c(rep(sqrt(n * nclass[l + 1]/(cumnclass[l] * cumnclass[l + 1])), l), -sqrt(n * cumnclass[l]/(nclass[l +
                                                                                                                       1] * cumnclass[l + 1])), rep(0, r - l))
  }

  # return Ytilde and theta matrix
  return(list(Ytilde = Z %*% theta, theta = theta))
}

# the function to compare the difference from two subspace
cv_cor = function(u, v) {
  idx = apply(u, 1, function(x) sum(is.na(x)) == 0) & apply(v, 1, function(x) sum(is.na(x)) == 0)

  u = scale(u[idx, , drop = F], center = T, scale = F)
  v = scale(v[idx, , drop = F], center = T, scale = F)

  Su = tcrossprod(u)
  Sv = tcrossprod(v)
  if (sum(sum(diag(Su)) != 0) == 0 | sum(sum(diag(Sv)) != 0) == 0)
    return(cor = 0)
  crossterm = sum(diag(Su %*% Sv))
  denominator = sqrt(sum(diag(Su %*% Su)) * sum(diag(Sv %*% Sv)))
  return(cor = sqrt(crossterm/denominator))
}

# calculate lambda such that parameters exceed these values will generate trivial model. For imput
# data with missing value
lambda_max <- function(Z, X_list, D, Z_nonNA_id, X_nonNA_id) {
  n = nrow(Z)
  Ytilde = matrix(NA, ncol = ncol(Z) - 1, nrow = nrow(Z))
  Ytilde[Z_nonNA_id, ] = transformYCpp(Z[Z_nonNA_id, ])

  # compute the union ID (Y and X)
  inter_ID_YX = lapply(X_nonNA_id, function(x) sort(intersect(Z_nonNA_id, x)))

  # center and scale x
  coef = lapply(X_list, function(x) apply(x, 2, sd, na.rm = TRUE) * sqrt((nrow(x) - 1)/nrow(x)))
  centeredX = lapply(1:length(X_list), function(i) scale(X_list[[i]], scale = coef[[i]], center = T))
  # calculate lambda_{max,d}
  lambda_t = rep(0, D)
  for (i in 1:D) {
    # n = length(inter_ID_YX[[i]])
    xSubset = centeredX[[i]][inter_ID_YX[[i]], , drop = F]
    ySubset = Ytilde[inter_ID_YX[[i]], , drop = F]
    lambda_t[i] = max(sqrt(rowSums(abs(crossprod(xSubset, ySubset))^2))/(n * D))
  }
  lambda_t
}

# assign cross-validation ID
assignID <- function(Z, X_list, D, nfolds) {
  n = nrow(Z)
  miss_pattern = matrix(0, ncol = D + 1, nrow = n)

  # record missing pattern for Z
  miss_pattern[, 1] = as.numeric(apply(Z, 1, function(x) sum(is.na(x)) == 0))
  for (i in 1:D) {
    miss_pattern[, i + 1] = as.numeric(apply(X_list[[i]], 1, function(x) sum(is.na(x)) == 0))
  }
  # find unique missing pattern
  unique_pattern = unique(miss_pattern)

  # generate the fold index
  id <- rep(-1, n)
  for (p_i in 1:nrow(unique_pattern)) {
    # if Z is not missing
    if (unique_pattern[p_i, 1] == 1) {
      for (i in 1:ncol(Z)) {
        # for each class, randomely assign id to each observation
        is_Z = (Z[, i] == 1)
        is_Z[is.na(is_Z)] = FALSE
        zNotMissing = apply(miss_pattern, 1, function(u) all(u == unique_pattern[p_i, ]))
        id[is_Z & zNotMissing] <- sample(rep(seq_len(nfolds), length.out = sum(is_Z & zNotMissing)))
      }
    } else {
      # if Z is missing, randomely assign id to each observation
      zNotMissing = apply(miss_pattern, 1, function(u) all(u == unique_pattern[p_i, ]))
      id[zNotMissing] <- sample(rep(seq_len(nfolds), length.out = sum(zNotMissing)))
    }
  }
  id
}


# define the evaluation function for cross validation for JACA with missing values
objectiveJACAMissing = function(W_list, test_z, train_x, test_x, D, theta, dis, alpha) {
  ## W_list: the results we get from JACA. D: number of data sets.
  Ytilde = test_z %*% theta
  n = nrow(test_z)

  if (D != 1) {
    # Xbar: the mean values of train_x
    Xbar = lapply(train_x, function(mat) colMeans(mat, na.rm = T))
    centeredX = lapply(1:D, function(num) test_x[[num]] - matrix(rep(Xbar[[num]], nrow(test_x[[num]])),
                                                                 byrow = T, nrow = nrow(test_x[[num]])))

    # the part1 of the objective function
    part1 = sum(sapply(1:D, function(num) -sum(cv_cor(Ytilde, centeredX[[num]] %*% W_list[[num]]))))
    # the part2 of the objective function
    part2 = sapply(combn(1:D, 2, simplify = F), function(com) -sum(cv_cor(centeredX[[com[1]]] %*%
                                                                            W_list[[com[1]]], centeredX[[com[2]]] %*% W_list[[com[2]]])))
    part2 = sum(part2)

    return((alpha * part1 + (1 - alpha) * part2/(D - 1)))
  } else {
    # Xbar: the mean values of train_x
    Xbar = lapply(train_x, function(mat) colMeans(mat))
    centeredX = lapply(1:D, function(num) test_x[[num]] - matrix(rep(Xbar[[num]], nrow(test_x[[num]])),
                                                                 byrow = T, nrow = nrow(test_x[[num]])))


    # the part1 of the objective function
    part1 = sum(sapply(1:D, function(num) -sum(cv_cor(Ytilde, centeredX[[num]] %*% W_list[[num]]))))
    # the part2 of the objective function
    part2 = 0


    return(alpha * part1 + (1 - alpha) * part2/(D - 1))
  }
}

jacaTrain_augmented = function(bigy, bigx, coef, D, p_n, lambda, rho, missing = F, alpha = 0.5, W_list = NULL, kmax = 500, eps = 1e-06,
                          verbose = F) {
  # check the number of lambda provided. lambda should be a non-negtive vector
  if (any(lambda < -eps))
    stop("lambda must be nonnegtive!")
  if (length(lambda) != D)
    stop("Check the number of lambda!")


  if (!is.null(W_list)) {
    W_list[[1]] = do.call(rbind, W_list)
  }

  lambda_vec = rep(lambda, p_n)

  if (nrow(bigx) != nrow(bigy))
    stop("Dimensions of X_list and Z don't match!")

  # fit model
  result = jacaCpp(Y = bigy, X_list = bigx, lambda = lambda_vec, rho = rho, W_list = W_list, kmax = kmax,
                   eps = eps, verbose = verbose)
  order_idx = c(0, cumsum(p_n))
  W_d = lapply(1:D, function(idx) diag(1/coef[[idx]], length(coef[[idx]])) %*% result[(order_idx[idx] + 1):(order_idx[idx + 1]), , drop = F])

  for (i in 1:D) {
    W_d[[i]][abs(W_d[[i]]) < eps] = 0
  }

  return(W_d)
}
