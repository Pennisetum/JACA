#' Solve Joint Association and Classification Analysis problem for multi-view data.
#'
#' Given Class indicator matrix Z and a list of matrices X_list,
#' \code{jacaTrain} estimates view-specific matrices of discriminant vectors by solving
#' Joint Association and Classification Analysis problem. It can also be used to perform semi-supervised learning,
#' that is to use information from both labeled and unlabeled subjects to construct classification rules.
#'
#' @param Z An N by K class indicator matrix; rows are samples and columns are class indicator vectors with z_k = 1 if observation belongs to class k.
#' @param X_list A list of input data matrices; in each sublist, rows are samples and columns are features.
#' @param lambda A vector of L1 penalty parameters; if there are D input data matrices, lambda should also contain D elements. \code{lambda} controls the sparsity of the solutions and must be between 0 and 1 (small L1 bound corresponds to less penalization).
#' @param rho  Scaler, l2 regularization penulty parameter on X'X/n. Should also  be between 0 and 1 (small L2 bound corresponds to less penalization).
#' @param missing Logical. If False, input data \code{X_list} must be complete and have no missing values.
#' If True, input data \code{X_list} should contain missing values.
#' @param alpha The parameter to control the weight between optimal scoring and CCA part. Default is 0.5.
#' @param W_list A list of inital guess of matrices of discriminant vectors (can be served as a warm start).
#' @param kmax Max iteration number.
#' @param eps  Threashold value used to determine the convergence of the optimization algorithm; the default value is 1e-06.
#' @param verbose Logical. If False, the algorithm will stay silent. If True, it will print out fitting progress.
#'
#' @return \item{W_d}{A list of view-specific matrices of discriminant vectors.}
#'
#' @examples
#' set.seed(1)
#' # Generate class indicator matrix Z
#' n = 100
#' Z=matrix(c(rep(1, n),rep(0, 2 * n)), byrow = FALSE, nrow = n)
#' for(i in 1:n){
#'   Z[i, ] = sample(Z[i, ])
#' }
#'
#' # Generate input data X_list
#' d = 2
#' X_list = sapply(1:d, function(i) list(matrix(rnorm(n * 20), n, 20)))
#'
#' # Train JACA model
#' W = jacaTrain(Z, X_list, lambda = rep(0.05, 2), verbose = FALSE, alpha= 0.5, rho = 0.2)
#'
#' # Show the number of non-zero rows of each matrix of discriminant vectors
#' sapply(W, function(x) sum(rowSums(x) != 0))
#'
#' # Test semi supervised learning
#' # Set certain class labels and subsets of views as missing
#' Z[90:100, ] = rep(NA, 3)
#' X_list[[1]][1:10, ] = NA
#' X_list[[2]][11:20, ] = NA
#' W = jacaTrain(Z, X_list, kmax = 200, eps = 1e-06, lambda = rep(0.05, 2),
#'               alpha = 0.5, rho = 0.2, missing = TRUE)
#'
#' # Show the number of non-zero rows of each matrix of discriminant vectors
#' sapply(W, function(x) sum(rowSums(x) != 0))
#'
#' @export
jacaTrain = function(Z, X_list, lambda, rho, missing = F, alpha = 0.5, W_list = NULL, kmax = 500, eps = 1e-06,
                      verbose = F) {

  # record the number of views
  if (is.list(X_list)) {
    D = length(X_list)
  } else {
    stop("Check your X!")
  }

  # check the number of lambda provided. lambda should be a non-negtive vector
  if (any(lambda < -eps))
    stop("lambda must be nonnegtive!")
  if (length(lambda) != D)
    stop("Check the number of lambda!")

  # record sample size
  n = nrow(X_list[[1]])

  # check the input W_list
  if (!is.null(W_list)) {
    for (i in 1:D) {
      if (nrow(W_list[[i]]) != ncol(X_list[[i]]))
        stop("Check W_list!")
    }
    W_list[[1]] = do.call(rbind, W_list)
  }

  # Generate augmented X and Y
  out = generateAugmentedXY(Z, X_list, alpha = alpha, missing = missing)

  # compute number of features in each data set
  p_n = sapply(X_list, ncol)
  # formulate lambda vector
  lambda_vec = rep(lambda, p_n)

  # fit model
  result = jacaCpp(Y = out$bigy, X_list = out$bigx, lambda = lambda_vec, rho = rho, W_list = W_list, kmax = kmax,
                   eps = eps, verbose = verbose)
  order_idx = c(0, cumsum(p_n))
  W_d = lapply(1:D, function(idx) diag(1/out$coef[[idx]], length(out$coef[[idx]])) %*% result[(order_idx[idx] + 1):(order_idx[idx + 1]), , drop = F])

  for (i in 1:D) {
    W_d[[i]][abs(W_d[[i]]) < eps] = 0
  }

  return(W_d)
}

#' Cross-validation function for JACA
#'
#' Chooses optimal tuning parameters lambda and rho for function \code{jacaTrain} using cross-validation.
#'
#' @param Z An N by K class indicator matrix; rows are samples and columns are class indicator vectors with z_k = 1 if observation belongs to class k.
#' @param X_list A list of input data matrices; in each sublist, rows are samples and columns are features.
#' @param nfolds Number of cross-validation folds.
#' @param lambda_seq The set of L1 penalty parameters to be considered. Should be chosen from 0 to 1.
#' The default value is NULL and \code{jacaCV} generates its own sequence.
#' @param n_lambda The number of \code{lambda_seq} considered. Used only if \code{lambda_seq = NULL}.
#' @param rho_seq The set of L2 penalty parameters to be considered. Should be chosen from 0 to 1.
#' @param n_rho The number of \code{rho_seq} considered. Used only if \code{rho_seq = NULL}.
#' @param missing Logical. If False, input data \code{X_list} must be complete and have no missing values.
#' If True, input data \code{X_list} should contain missing values.
#' @param alpha The parameter to control the weight between optimal scoring and CCA part. Default is 0.5.
#' @param W_list A list of inital guess of matrices of discriminant vectors (can be served as a warm start).
#' @param kmax Max iteration number.
#' @param eps  Threshold value used to determine the convergence of the optimization algorithm; the default value is 1e-06.
#' @param verbose Logical. If False, the algorithm will stay silent. If True, it will print out fitting progress.
#' @param foldID User-supplied fold seperation. The default is NULL, and \code{jacaCV} generates its own folds.

#' @return \item{W_min}{A list of view-specific matrices of discriminant vectors. Generated using the parameters chosen by cross-validation method.}
#' @return \item{lambda_min}{The value of L1 penalty parameter that resulted in the minimal mean cross-validation error.}
#' @return \item{rho_min}{The value of L2 penalty parameter that resulted in the minimal mean cross-validation error.}
#' @return \item{grid_seq}{The matrix of tuning parameters used.}
#' @return \item{error_mean}{The mean cross-validated error of each combination of the tuning parameters.}
#' @return \item{error_se}{The standard error of cross-validated error of each combination of the tuning parameters.}
#'
#' @examples
#' set.seed(1)
#' # Generate class indicator matrix Z
#' n = 500
#' Z = matrix(c(rep(1, n),rep(0, 2 * n)), byrow = FALSE, nrow = n)
#' for(i in 1:n){
#'   Z[i, ] = sample(Z[i, ])
#' }
#'
#' # Generate input data X_list
#' d = 2
#' X_list = sapply(1:d, function(i) list(matrix(rnorm(n * 20), n, 20)))
#' id <- 1:nrow(Z)
#'
#' # Train JACA model using cross validation
#' result = jacaCV(Z, X_list, nfolds = 3, lambda_seq = c(0.02, 0.04), rho_seq = c(0.3, 0.6))$W_min
#'
#' # Test semi supervised learning
#' # Set certain class labels and subsets of views as missing
#' Z[90:100, ] = rep(NA, 3)
#' X_list[[1]][1:10, ] = NA
#' X_list[[2]][11:20, ] = NA
#'
#' # Train JACA model using cross validation
#' result = jacaCV(Z, X_list, nfolds = 3, lambda_seq = c(0.02, 0.04),
#'                 rho_seq = c(0.3, 0.6), missing = TRUE)$W_min
#'
#' @export

jacaCV <- function(Z, X_list, nfolds = 5, lambda_seq = NULL, n_lambda = 50, rho_seq = seq(0.01, 1, length = 20),
                   n_rho = 5, missing = F, alpha = 0.5, W_list = NULL, kmax = 500, eps = 1e-06, verbose = F, foldID = NULL) {
  # record the number of views
  if (is.list(X_list)) {
    D = length(X_list)
  } else {
    stop("Check your X!")
  }

  # compute maximum lambda allowed
  if (missing) {
    # Record the missing ID
    Z_nonNA_id = which(apply(Z, 1, function(x) sum(is.na(x)) == 0))
    X_nonNA_id = lapply(X_list, function(sublist) {
      which(apply(sublist, 1, function(x) sum(is.na(x)) == 0))
    })
    # compute the union ID (Y and X)
    inter_ID_YX = lapply(X_nonNA_id, function(x) sort(intersect(Z_nonNA_id, x)))
    # Compute maximum lambda
    lambda_t = lambda_max(Z = Z, X_list = X_list, D = D, Z_nonNA_id = Z_nonNA_id, X_nonNA_id = X_nonNA_id) * alpha
  } else {
    if (any(anyNA(Z), sapply(X_list, anyNA))) {
      stop("Datasets contain missing value! Please set missing = T")
    }
    # caculate transformed Ytilde and scale X
    Ytilde = transformYCpp(Z)
    centeredX = list()
    n = nrow(Z)
    coef = lapply(X_list, function(x) apply(x, 2, sd, na.rm = TRUE) * sqrt((n - 1)/n))

    # calculate maximum of lambda
    lambda_t = rep(0, D)
    for (i in 1:D) {
      centeredX[[i]] = scale(X_list[[i]], scale = coef[[i]], center = T)
      lambda_t[i] = max(sqrt(rowSums(abs(crossprod(centeredX[[i]], Ytilde))^2))/(n * D)) * alpha
    }
  }

  # If lambda_seq is supplied, the function should only keep values that satisfy 0 < lambda_i <
  # lambda_max. Otherwise, generate lambda_seq
  if (!is.null(lambda_seq)) {
    lambda_seq = lambda_seq[(lambda_seq >= 0) & (lambda_seq <= 1)]
    if (length(lambda_seq) == 0) {
      warning("Check input lambda_seq")
      # generate lambda_seq
      lambda_seq <- exp(seq(log(1), log(1e-01), length.out = n_lambda))
    }
  } else lambda_seq <- exp(seq(log(1), log(1e-01), length.out = n_lambda))
  n_lambda = length(lambda_seq)

  # If rho is supplied, the function should only keep values that satisfy 0<=rho<=1
  if (!is.null(rho_seq)) {
    rho_seq = rho_seq[(rho_seq >= 0) & (rho_seq <= 1)]
    rho_seq <- sort(rho_seq)
  } else {
    stop("Check rho_seq")
  }
  n_rho = length(rho_seq)
  p_n = sapply(X_list, ncol)

  if (missing) {
    ## generate the fold index
    if (is.null(foldID)) {
      id <- assignID(Z = Z, X_list = X_list, D = D, nfolds = nfolds)
    } else {
      id <- foldID
    }
  } else {
    ## generate the fold index
    if (is.null(foldID)) {
      id <- 1:nrow(Z)
      for (i in 1:ncol(Z)) {
        id[Z[, i] == 1] <- sample(rep(seq_len(nfolds), length.out = sum(Z[, i] == 1)))
      }
    } else {
      id <- foldID
    }
  }

  #### calculate the normal JACA model.  generate the grid to search the parameters.
  init_grid = expand.grid(lambda_seq, rho_seq)
  cv = matrix(0, nfolds, nrow(init_grid))

  fit_model_old = W_list
  # cross validation part
  for (i in 1:nfolds) {
    ## during each iteration, select one as test, the rest (K-1) as train
    test_idx = which(id == i)
    test_x = lapply(X_list, function(X) X[test_idx, ])
    test_z = Z[test_idx, ]
    train_x = lapply(X_list, function(X) X[-test_idx, ])
    train_z = Z[-test_idx, ]
    theta = transformY(train_z)$theta

    # generate augmented x and augmented y
    out = generateAugmentedXY(Z = train_z, X_list = train_x, alpha = alpha, missing = missing)

    # train JACA model
    for (cv_ind in 1:nrow(init_grid)){
      fit_model = jacaTrain_augmented(bigy = out$bigy, bigx = out$bigx, coef = out$coef, D = D, p_n = p_n,  lambda = init_grid[cv_ind, 1] * lambda_t, rho = init_grid[cv_ind, 2], missing = missing, alpha = alpha, W_list = fit_model_old, eps = eps, kmax = kmax, verbose = F)
      cv[i, cv_ind] = objectiveJACA(W_list = fit_model, test_z = test_z, train_x = train_x, test_x = test_x, D = D,
                                    theta = theta, alpha = alpha)
      if (verbose){
        print(paste0("Parameters: ", round(init_grid[cv_ind, 1] * lambda_t, 3), ", ", init_grid[cv_ind, 2], ". Objective (sum corr): ", -cv[i, cv_ind]))
      }
      fit_model_old = fit_model
    }
    cat("Complete", i, "\n")
  }

  # compute mean cv errors and choose parameters
  cvm = colSums(cv)/nfolds
  cvse = apply(cv, 2, sd)/sqrt(nfolds)
  min_id = which.min(cvm)
  grid_min = init_grid[min_id, ]

  # compute final W_list
  W_min = jacaTrain(Z = Z, X_list = X_list, lambda = as.numeric(grid_min[1]) * lambda_t, rho = as.numeric(grid_min[2]),
                    missing = missing, alpha = alpha, verbose = F, W_list = W_list, kmax = kmax, eps = eps)


  # selected parameters
  lambda_min = as.numeric(grid_min[1])
  rho_min = as.numeric(grid_min[2])


  return(list(W_min = W_min, lambda_min = lambda_min, rho_min = rho_min, grid_seq = init_grid, error_mean = cvm,
              error_se = cvse))
}

#' Classification method for JACA model
#'
#' Classify observations in the test set based on the supplied matrices of canonical vectors \code{W_list} and the training datasets.
#' This function first fits a linear discriminant analysis model based on \code{W_list}, \code{trainx} and \code{trainz}, then makes predictions
#' by \code{W_list}, \code{testx}.
#'
#' @param W_list A list of view-specific matrices of discriminant vectors. Should be the results of \code{jacaTrain} or \code{jacaCV}.
#' @param trainx A list of input data matrices that are used to generate the model: in each sublist, samples are rows and columns are features.
#' @param trainz An N by K class indicator matrix that are used to generate the model; rows are samples and columns are class indicator vectors with z_k = 1 if observation belongs to class k.
#' @param testx A list of input data matrices Predictions will be made using it.
#' @param posterior Logical. If TRUE, the function outputs the posterior probabilities for the classes, otherwise it will outputs the class assignments results.
#'
#' @return \item{prediction}{An n by D matrix with predicted group labels for the test set if posterior = FALSE. Columns are predictions made by each dataset seperately. Otherwise it is a list of posterior probabilities for the test set.}
#'
#' @examples
#' set.seed(1)
#' # Generate class indicator matrix Z
#' n = 100
#' Z=matrix(c(rep(1, n),rep(0, 2 * n)), byrow = FALSE, nrow = n)
#' for(i in 1:n){
#'   Z[i, ] = sample(Z[i, ])
#' }
#'
#' # Generate input data X_list
#' d = 2
#' X_list = sapply(1:d, function(i) list(matrix(rnorm(n * 20), n, 20)))
#'
#' # Train the model
#' result = jacaTrain(Z, X_list, lambda = rep(0.05, 2), verbose = FALSE, alpha= 0.5, rho = 0.2)
#'
#' # Make predictions by each dataset seperately
#' zTest = predictJACA(result, X_list, Z, X_list)
#'
#' # Make predictions by the concatenated dataset
#' wJoint = list()
#' wJoint[[1]] = do.call(rbind, result)
#' xTrain = list()
#' xTrain[[1]] = do.call(cbind, X_list)
#' zTestConcatenated = predictJACA(wJoint, xTrain, Z, xTrain)
#'
#' @export
predictJACA = function(W_list, trainx, trainz, testx, posterior = F) {
  Y.train = apply(trainz, 1, which.max) - 1
  n = nrow(trainz)
  D = length(testx)
  Xbar = lapply(trainx, function(mat) colMeans(mat))
  centeredX = lapply(1:D, function(num) testx[[num]] - matrix(rep(Xbar[[num]], nrow(testx[[num]])),
                                                              byrow = T, nrow = nrow(testx[[num]])))
  Y.pred = list()
  for (i in 1:D) {
    pred.data = scale(trainx[[i]], scale = F) %*% W_list[[i]]
    nonzero_col = apply(pred.data, 2, var) > 1e-06
    if (sum(nonzero_col) != 0) {
      pred.data = data.frame(Y.train = Y.train, pred.data[, nonzero_col])
      colnames(pred.data)[-1] <- paste("X", 1:(ncol(pred.data) - 1), sep = "")
      lda.model = MASS::lda(Y.train ~ ., data = pred.data)
      temp = (centeredX[[i]] %*% W_list[[i]])[, nonzero_col, drop = F]
      colnames(temp) <- paste("X", 1:(ncol(pred.data) - 1), sep = "")
      if(posterior)
        Y.pred[[i]] = predict(lda.model, data.frame(temp))$posterior
      else
        Y.pred[[i]] = as.numeric(as.character(predict(lda.model, data.frame(temp))$class))
    } else {
      Y.pred[[i]] = sample(0:max(Y.train), nrow(testx[[1]]), replace = T)
      Y.pred = do.call(cbind, Y.pred)
      warning("Constant canonical vectors!")
    }
  }


  return(prediction = Y.pred)
}



# The evaluation function for cross validation when alpha = 1 (error rate only)
objectiveJACAone = function(W_list, test_z, train_x, test_x, D, theta) {
  ## W_list: the results we get from JACA. D: number of data sets.
  Ytilde = test_z %*% theta
  n = nrow(test_z)

  # Xbar: the mean values of train_x
  Xbar = lapply(train_x, function(mat) colMeans(mat, na.rm = T))
  centeredX = lapply(1:D, function(num) test_x[[num]] - matrix(rep(Xbar[[num]], nrow(test_x[[num]])), byrow = T, nrow = nrow(test_x[[num]])))

  # Use correlation as cv criterion the part1 of the objective function
  part1 = sum(sapply(1:D, function(num) -sum(cv_cor(Ytilde, centeredX[[num]] %*% W_list[[num]]))))
  # the part2 of the objective function

  return(part1)
}


#' Alternative cross-validation function for JACA tailored towards minimization of error rates (allows to also cross-validate for alpha)
#'
#' Chooses optimal tuning parameters lambda, rho and alpha for function \code{jacaTrain} based on cross-validation to minimize out-of-sample misclassifiation error rate (measured by objective function in JACA with alpha = 1).
#'
#' @inheritParams jacaCV
#' @param alpha_seq A sequence of alpha values to consider, where alpha must be strictly greater than 0, and less or equal to 1.
#'
#' @return \item{W_min}{A list of view-specific matrices of discriminant vectors. Generated using the parameters chosen by cross-validation method.}
#' @return \item{lambda_min}{The value of L1 penalty parameter that resulted in the minimal mean cross-validation error.}
#' @return \item{rho_min}{The value of L2 penalty parameter that resulted in the minimal mean cross-validation error.}
#' @return \item{alpha_min}{The value of alpha parameter that resulted in the minimal mean cross-validation error.}
#' @return \item{grid_seq}{The matrix of tuning parameters used.}
#' @return \item{error_mean}{The mean cross-validated error of each combination of the tuning parameters.}
#' @return \item{error_se}{The standard error of cross-validated error of each combination of the tuning parameters.}
#' @export
#'
#' @examples
#' #' set.seed(1)
#' # Generate class indicator matrix Z
#' n = 20
#' Z=matrix(c(rep(1, n),rep(0, 2 * n)), byrow = FALSE, nrow = n)
#' for(i in 1:n){
#'   Z[i, ] = sample(Z[i, ])
#' }
#'
#' # Generate input data X_list
#' d = 2; p=10
#' X_list = sapply(1:d, function(i) list(matrix(rnorm(n * p), n, p)))
#'
#' # Train JACA model using cross validation to minimize the error rate
#' result = jacaCVerror(Z, X_list, nfolds = 3, lambda_seq = c(0.02, 0.04), rho_seq = c(0.3, 0.6))
#'
jacaCVerror <- function(Z, X_list, nfolds = 5, lambda_seq = NULL, n_lambda = 50, rho_seq = seq(0.01, 1, length = 5),
                        n_rho = 5, missing = F, alpha_seq = c(0.5, 0.7), W_list = NULL, kmax = 500, eps = 1e-06, verbose = F, foldID = NULL) {
  # record the number of views
  D = length(X_list)
  p_n = sapply(X_list, ncol)

  if (missing) {
    # Record the missing ID
    Z_nonNA_id = which(apply(Z, 1, function(x) sum(is.na(x)) == 0))
    X_nonNA_id = lapply(X_list, function(sublist) {
      which(apply(sublist, 1, function(x) sum(is.na(x)) == 0))
    })
    # compute the union ID (Y and X)
    inter_ID_YX = lapply(X_nonNA_id, function(x) sort(intersect(Z_nonNA_id, x)))

    # Compute maximum lambda
    lambda_t = lambda_max(Z = Z, X_list = X_list, D = D, Z_nonNA_id = Z_nonNA_id, X_nonNA_id = X_nonNA_id)
  } else {
    if (any(anyNA(Z), sapply(X_list, anyNA))) {
      stop("Datasets contain missing value! Please set missing = T")
    }
    # caculate transformed Ytilde and scale X
    Ytilde = transformY(Z)$Ytilde
    centeredX = list()
    n = nrow(Z)
    coef = lapply(X_list, function(x) apply(x, 2, sd, na.rm = TRUE) * sqrt((n - 1)/n))

    # calculate maximum of lambda
    lambda_t = rep(0, D)
    for (i in 1:D) {
      centeredX[[i]] = scale(X_list[[i]], scale = coef[[i]], center = T)
      lambda_t[i] = max(sqrt(rowSums(abs(crossprod(centeredX[[i]], Ytilde))^2))/(n * D))
    }
  }

  # If lambda_seq is supplied, the function should only keep values that satisfy 0 < lambda_i <
  # lambda_max. Otherwise, generate lambda_seq
  if (is.null(lambda_seq)) {
    lambda_seq <- exp(seq(log(1), log(1e-01), length.out = n_lambda))
  }
  n_lambda = length(lambda_seq)

  # If rho is supplied, the function should only keep values that satisfy 0<=rho<=1
  if (!is.null(rho_seq)) {
    rho_seq = rho_seq[(rho_seq >= 0) & (rho_seq <= 1)]
    rho_seq <- sort(rho_seq)
  } else {
    stop("Check rho_seq")
  }
  n_rho = length(rho_seq)

  # Make sure that supplied alphas are not too close to 0, and are not above 1
  if (any(alpha_seq <= 0.005) & any(alpha_seq > 1)){
    stop("Check alpha_seq")
  }
  n_alpha = length(alpha_seq)

  if (missing) {
    ## generate the fold index
    if (is.null(foldID)) {
      id <- assignID(Z = Z, X_list = X_list, D = D, nfolds = nfolds)
    } else {
      id <- foldID
    }
  } else {
    ## generate the fold index
    if (is.null(foldID)) {
      id <- 1:nrow(Z)
      for (i in 1:ncol(Z)) {
        id[Z[, i] == 1] <- sample(rep(seq_len(nfolds), length.out = sum(Z[, i] == 1)))
      }
    } else {
      id <- foldID
    }
  }

  #### calculate the normal JACA model.  generate the grid to search the parameters.
  init_grid = expand.grid(lambda_seq, rho_seq, alpha_seq) #first go through all lambdas, then through rhos, then through alphas

  # Figure out in advance the indexes that correspond to a change in alpha_seq
  index_alpha = (0:(n_alpha - 1)) * (n_lambda * n_rho) + 1

  cv = matrix(0, nfolds, nrow(init_grid))
  fit_model_old = W_list

  # cross validation part
  for (i in 1:nfolds) {
    ## during each iteration, select one as test, the rest (K-1) as train
    test_idx = which(id == i)
    test_z = Z[test_idx, ]
    test_x = lapply(X_list, function(X) X[test_idx, ])


    # Some of these will be NA so for misclassification error rate remove them. Also if check for error rate, can only do on complete cases test data
    # isnotna_test = (rowSums(is.na(test_z)) == 0) & (rowSums(is.na(do.call("cbind", test_x))) == 0)
    # test_z = test_z[isnotna_test, ]
    # test_x = lapply(test_x, function(x) x[isnotna_test, ])

    # While I can fit the model on full data, I can only do full classification on complete data
    train_x = lapply(X_list, function(X) X[-test_idx, ])
    train_z = Z[-test_idx, ]

    # Get training theta, only need to construct once
    theta = transformY(train_z)$theta

    for (cv_ind in 1:nrow(init_grid)){
      # Generate augmented data matrix with that particular alpha at each cv_ind where alpha value changes
      if (cv_ind %in% index_alpha){
        out = generateAugmentedXY(Z = train_z, X_list = train_x, alpha = init_grid[cv_ind, 3], missing = missing)
      }

      # Update the model
      fit_model = jacaTrain_augmented(bigy = out$bigy, bigx = out$bigx, coef = out$coef, D = D, p_n = p_n,  lambda = init_grid[cv_ind, 1] * lambda_t * init_grid[cv_ind, 3], rho = init_grid[cv_ind, 2], missing = missing, alpha = init_grid[cv_ind, 3], W_list = fit_model_old, eps = eps, kmax = kmax, verbose = F)

      # Get the error
      cv[i, cv_ind] = objectiveJACAone(W_list = fit_model, test_z = test_z, train_x = train_x, test_x = test_x, D = D, theta = theta)

      # Save model fit
      fit_model_old = fit_model
    }

    # Finish this fold
    cat("Complete", i, "\n")
  }

  # compute mean cv errors and choose parameters
  cvm = colSums(cv)/nfolds
  cvse = apply(cv, 2, sd)/sqrt(nfolds)


  # Average over all minimums
  grid_select = init_grid[cvm == min(cvm), ]
  grid_min = colMeans(grid_select)

  # compute final W_list
  W_min = jacaTrain(Z = Z, X_list = X_list, lambda = as.numeric(grid_min[1]) * lambda_t * as.numeric(grid_min[3]), rho = as.numeric(grid_min[2]),
                    missing = missing, alpha = as.numeric(grid_min[3]), verbose = F, W_list = W_list, kmax = kmax, eps = eps)


  # selected parameters
  lambda_min = as.numeric(grid_min[1])
  rho_min = as.numeric(grid_min[2])
  alpha_min = as.numeric(grid_min[3])


  return(list(W_min = W_min, lambda_min = lambda_min, rho_min = rho_min, grid_seq = init_grid, error_mean = cvm, error_se = cvse, alpha_min = alpha_min))
}

