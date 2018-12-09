#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

/*
##########################################################
Notice!!! The following codes can only be called in c++.
##########################################################
*/

// Generate a vector that repeat u n times.
// Input:
//  u  ------ the value to be repeated
//  n  ------ an integer-valued (non-negative) number. Repeat times
colvec repc(const double& u, const int& n) {
  vec v(n);
  v.fill(u);
  return(v);
}

// Calculate the matrix multiplication. Return Z^\top \times X
mat crossprodc(const mat& Z,const mat& X){
  mat crossproduct(Z.n_cols, X.n_cols);
  for (int i=0; i<Z.n_cols; i++){
    for (int j=0; j<X.n_cols; j++){
      crossproduct(i, j) = dot(Z.col(i), X.col(j));
    }
  }
  return crossproduct;
}

// objective function of JACA
double objectiveJACA(const mat& Ytilder, const mat& Xtilder, const mat& W, vec lambda, double rho){
  double value = accu(square(Ytilder - Xtilder*W))/(2) + accu(sqrt(lambda % sum(square(W),1))) +
    rho/2*accu(square(W)) - rho/(2)*accu(square(Xtilder*W));
  return value;
}

// main function of dLDA for cleaned data
// Input:
//  Ytilder --- Transformed reponse matrix
//  Xtilder --- Multi-view data sets
//  lambda ---- penalty parameter vector
//  rhor ------ l2 penalty parameter
//  W0   ------ Initial guess for the solution
//  kmax ------ Maximum iteration count
//  eps  ------ Threashold when compare the difference
//  verbose --- logical. Print the information of performance
mat blockCoorDescent(const mat& Ytilder, const mat& Xtilder, vec lambda, double rhor,
            Nullable<NumericMatrix> W0 = R_NilValue, int kmax = 500,
            double eps = 0.000001, bool verbose = false){

  // Declare variables. Count the number of observation and number of views
  int nXcol = Xtilder.n_cols;
  int nYcol = Ytilder.n_cols;
  mat W0r(nXcol,nYcol), W1r, r;
  double newv, oldv;
  rowvec Temp;
  // Initialize the iteration counter.
  int k = 0;
  // check if W0 is given
  // initialize W_0r
  if (W0.isNotNull()){
    W0r = as<mat>(W0);
  } else {
    W0r.fill(0.0);
  }
  // W1r is the updated version of W0r
  W1r = W0r;

  if (rhor == 1){
    // speciel function when rhor == 1, since we don't need to loop
    r = crossprodc(Xtilder, Ytilder);
    // update each coordinate, follows the pseudocode in the proposal with rhor=1. Then some terms disappear.
    for (int j=0; j < nXcol;j++){
      // This is the same in r:T = r[j,]
      Temp = r.row(j);
      W1r.row(j) = max(0.0,1.0-lambda[j]/(sqrt(sum(square(Temp)))))*Temp;
    }
    return W1r;
  }

  // compute r (residure) in the proposal
  r = Ytilder - (1-rhor)*Xtilder*W0r;
  // block coordinate descent
  newv = objectiveJACA(Ytilder, Xtilder, W1r, lambda, rhor);
  oldv = newv-1.0;
  while (abs(oldv-newv)>eps){
    // update the counter
    oldv = newv;
    k=k+1;
    W0r=W1r;

    // update each row of W1r
    for(int j=0; j<nXcol;j++){
      // Temp = x_j^\top r + \|x_j\|_2^2 w^{(k-1)}_j
      Temp = crossprodc(Xtilder.col(j), r) + (1-rhor)*accu(square(Xtilder.col(j)))*W0r.row(j);
      if(accu(square(Temp))==0){
        W1r.row(j) = Temp;
      }else {
        W1r.row(j) = max(0.0, 1.0-lambda[j]/(sqrt(accu(square(Temp))))) *
          Temp/((1-rhor)*accu(square(Xtilder.col(j)))+rhor);}
      r=r+(1-rhor)*(Xtilder.col(j)*(W0r.row(j)-W1r.row(j)));
    }
    newv=objectiveJACA(Ytilder, Xtilder, W1r, lambda, rhor);

    // output the trace
    if (verbose){
      Rcout<<max(max(abs(W1r-W0r)))<<", "<< oldv-newv<<endl;
    }

    // stop the function if exceed the max num of loops
    if (k>kmax){
      Rcout<<"Inner loop fails to converge"<<endl;
      break;
    }
  }
  return W1r;
}


/*
###############################################################################################################
###############################################################################################################
Notice!!! The following functions can be called in R.                                           ###############
###############################################################################################################
###############################################################################################################
*/

/* This function use dLDA procedure to solve the JACA function
* It should inherit the parameters from JACA_Cpp function
*
*/

// [[Rcpp::export]] //
NumericMatrix transformYCpp(NumericMatrix Zr){
  //calculate the size of each class and the cumulative of sizes of classes (for later use)
  mat Z = as<mat>(Zr);
  rowvec nclass = sum(Z);
  rowvec cumnclass = cumsum(nclass);

  //number of observations and number of class minus 1
  int n = Z.n_rows;
  int r = Z.n_cols-1;

  //calculate the Theta matrix
  mat theta(r+1,r);
  theta.fill(0.0);
  for (int i = 0; i < r; i++ ) {
    theta.submat(0, i, i, i) = repc(sqrt(n*nclass(i+1)/(cumnclass(i)*cumnclass(i+1))), i+1);
    theta(i+1, i)= -sqrt(n*cumnclass(i)/(nclass(i+1)*cumnclass(i+1)));
  }
  mat Ytilde = Z*theta;
  return wrap(Ytilde);
}

// [[Rcpp::export]] //
NumericMatrix jacaCpp(NumericMatrix Y, NumericMatrix X_list, NumericVector lambda, double rho,
             Nullable<List> W_list = R_NilValue, int kmax = 500, double eps = 0.000001, bool verbose = false){
  // caculate transformed Ytilde
  mat Ytilder   = as<mat>(Y);
  mat centeredX = as<mat>(X_list);
  vec lambdar   = as<vec>(lambda);
  List  W_listr;
  mat W_d;

  // call main function to compute solutions
  if (W_list.isNotNull()){
    W_listr = clone(W_list);
    NumericMatrix W_0r = as<NumericMatrix>(W_listr[0]);
    W_d= blockCoorDescent(Ytilder, centeredX, lambdar, rho, W_0r,
                    kmax, eps, verbose);
  } else {
    W_d = blockCoorDescent(Ytilder, centeredX, lambdar, rho, R_NilValue,
                     kmax, eps, verbose);
  }

  return wrap(W_d);
}
