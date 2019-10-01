#include <RcppDist.h>
#include <Rcpp/Benchmark/Timer.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
// [[Rcpp::export]]
bool stab_test(const mat B_mat, const int N, const int P) {
  bool stab; // return value
  int size = N * P;
  int size_wo_lag1 = N * (P - 1);
  mat Bcomp = zeros(size, size); 
  //Rcout << Bcomp << std::endl << B
  Bcomp.rows(0, N - 1) = B_mat.submat(0, 0, N - 1, size - 1); //deleting the intercept: size - 1
  Bcomp.submat(N, 0, size - 1, size_wo_lag1 - 1) = eye(size_wo_lag1, size_wo_lag1);
  cx_vec eigenval = eig_gen(Bcomp);
  if (max(abs(eigenval)) < 1) {
    stab = true;
  } else {
    stab = false;
  }
  return stab;
}

/*
// [[Rcpp::export]]
cx_vec eigenv(const mat B_mat) {
  cx_vec eigenval = eig_gen(B_mat);
  return eigenval;
}
 */

// [[Rcpp::export]]
List carterkohn_tva_tvb(mat Y, mat X, mat a_elements, mat hlast, mat Q, vec b00, mat p00, vec mu, mat Fmat) {
  int size_b = b00.size(); //a_elemnts and hlast will form time-varying R
  int N = Y.n_cols; //number of variables
  int P = (size_b - N) / pow(N, 2); //number of lags for stability test and Bmat ( - N in numerator to correct for N intercepts in vec(b0))
  int nobs = Y.n_rows;
  mat beta_fw(nobs, size_b); // Kalman filter result
  cube p_fw(size_b, size_b, nobs);
  mat beta_bw(nobs, size_b); // Carter Kohn: to be returned
  mat errors(nobs, N); // Carter Kohn: to be returned
  vec b11 = b00; // starting values
  mat p11 = p00; // starting values
  mat diag_mat = eye(N, N); //for creating Hmat with Kronecker Product
  mat Amat = eye(N, N); //covariance terms will be lower triangular
  mat F_transp = trans(Fmat);
  //Kalman filter
  for (int i = 0; i < nobs; i++) {
    mat Hmat = kron(diag_mat, X.row(i)); //X expanded to allow for multiplying with vec(b)
    mat H_transp = trans(Hmat);
    //create R matrix: first create A, then R with hlast
    int a_counter = 0;
    for (int i_row = 0; i_row < N; i_row++) { // creating Amatrix: lower triangular
      for (int i_col = 0; i_col < i_row; i_col++) {
        Amat(i_row, i_col) = a_elements(i, a_counter);
        a_counter++;
      }
    }
    Amat = inv(Amat);
    mat R = Amat * diagmat(hlast.row(i + 1)) * trans(Amat); // i + 1 because of starting value of hlast which has to be ignored
    // Kalman steps
    vec b10 = mu + Fmat * b11;
    mat p10 = Fmat * p11 * F_transp + Q;
    vec eta10 = trans(Y.row(i)) - Hmat * b10;
    mat f10 = Hmat * p10 * H_transp + R;
    mat Kgain = p10 * H_transp * inv(f10);
    b11 = b10 + Kgain * eta10;
    p11 = p10 - Kgain * Hmat * p10;
    beta_fw.row(i) = trans(b11);
    p_fw.slice(i) = p11;
  }
  // Carter Kohn
  bool chk = false; //checking for stability
  while (!chk) {
    int stab = 1; // stability result
    rowvec B_sample = rmvnorm(1, trans(beta_fw.row(nobs - 1)), p_fw.slice(nobs - 1));
    mat B_mat = trans(reshape(B_sample, N * P + 1, N)); //for stability test and error/residual computation
    stab = stab * stab_test(B_mat, N, P); // if it's stable than function returns true
    beta_bw.row(nobs - 1) = B_sample;
    errors.row(nobs - 1) = Y.row(nobs - 1) - X.row(nobs - 1) * trans(B_mat); //residuals
    for (int i = nobs - 2; i >= 0; i--) {
      rowvec next_b = beta_bw.row(i + 1);
      rowvec curr_b = beta_fw.row(i);
      mat curr_p = p_fw.slice(i);
      //Rcout << trans(mu - Fmat * trans(curr_b)) << std::endl << next_b;
      vec eta_bw = trans(next_b) - mu - Fmat * trans(curr_b);
      mat f_bw = Fmat * curr_p * F_transp + Q;
      mat Kgain = curr_p * F_transp * inv(f_bw);
      vec b_update = trans(curr_b) + Kgain * eta_bw; //mean of cond distribution
      mat p_update = curr_p - Kgain * Fmat * curr_p; //variance of cond distribution
      rowvec B_sample = rmvnorm(1, b_update, p_update); //drawing the beta value
      mat B_mat = trans(reshape(B_sample, N * P + 1, N)); 
      stab = stab * stab_test(B_mat, N, P);
      if (stab == 0)
        break; //some matrix is not stable -> we start all over again
      beta_bw.row(i) = B_sample;
      errors.row(i) =  Y.row(i) - X.row(i) * trans(B_mat);
    }
    if (stab == 1)
      chk = true;
  }
  return List::create(Named("beta_bw") = beta_bw, Named("errors") = errors);
}


// [[Rcpp::export]]
List carterkohn_ca_cb(mat Y, mat X, mat Amat, mat hlast, mat Q, vec b00, mat p00, vec mu, mat Fmat, int maxtry = 1e3) {
  int size_b = b00.size(); //a_elemnts and hlast will form time-varying R
  int N = Y.n_cols; //number of variables
  int P = (size_b - N) / pow(N, 2); //number of lags for stability test and Bmat ( - N in numerator to correct for N intercepts in vec(b0))
  int nobs = Y.n_rows;
  int size = N * P; //next 4 lines for stability testing
  int size_wo_lag1 = N * (P - 1);
  int attempts; //for return only
  mat A_eps = trans(Amat); //needed for ortho errors; calculated here before inv() of A matrix
  mat Bcomp = zeros(size, size); 
  Bcomp.submat(N, 0, size - 1, size_wo_lag1 - 1) = eye(size_wo_lag1, size_wo_lag1);
  // mat beta_fw(nobs, size_b); // Kalman filter result
  // cube p_fw(size_b, size_b, nobs);
  // mat beta_bw(nobs, size_b); // Carter Kohn: to be returned
  // mat errors(nobs, N); // Carter Kohn: to be returned
  vec b11 = b00; // starting values
  mat p11 = p00; // starting values
  mat diag_mat = eye(N, N); //for creating Hmat with Kronecker Product
  mat F_transp = trans(Fmat);
  /*
  mat Amat = eye(N, N); //covariance terms will be lower triangular
  int a_counter = 0;
  for (int i_row = 0; i_row < N; i_row++) { // creating Amatrix: lower triangular
    for (int i_col = 0; i_col < i_row; i_col++) {
      Amat(i_row, i_col) = a_elements(a_counter);
      a_counter++;
    }
  }
  */
  Amat = inv(Amat);
  mat Amat_transp = trans(Amat);
  //Kalman filter
  for (int i = 0; i < nobs; i++) {
    mat Hmat = kron(diag_mat, X.row(i)); //X expanded to allow for multiplying with vec(b)
    mat H_transp = trans(Hmat);
    //create R matrix: first create A, then R with hlast
    //SHOULD BE i+1!!!!!
    mat R = Amat * diagmat(hlast.row(i + 1)) * Amat_transp; // i + 1 because of starting value of hlast which has to be ignored
    // Kalman steps
    vec b10 = mu + Fmat * b11;
    mat p10 = Fmat * p11 * F_transp + Q;
    vec eta10 = trans(Y.row(i)) - Hmat * b10;
    mat f10 = Hmat * p10 * H_transp + R;
    mat Kgain = p10 * H_transp * inv(f10);
    b11 = b10 + Kgain * eta10;
    p11 = p10 - Kgain * Hmat * p10;
    // beta_fw.row(i) = trans(b11);
    // p_fw.slice(i) = p11;
  }
  // Carter Kohn step
  // bool chk = false; //checking for stability
  //bool stab;
  bool stab = false;
  mat B_mat;
  rowvec B_sample;
  cx_vec eigenval;
  //int tries = 1;
  //while (!chk) {
  for (int i = 0; i < maxtry; i++) {
    // int stab = 1; // stability result
    B_sample = rmvnorm(1, b11, p11);
    B_mat = trans(reshape(B_sample, N * P + 1, N)); //for stability test and error/residual computation
    Bcomp.rows(0, N - 1) = B_mat.submat(0, 0, N - 1, size - 1); //deleting the intercept: size - 1
    eigenval = eig_gen(Bcomp);
    if (max(abs(eigenval)) < 1) {
      // chk = true;
      stab = true;
      attempts = i + 1;
      // Rcout << i;
      break;
    }
    /*
    else if (tries == maxtry) {
      chk = true;
      stab = false;
    }
    */
    //tries++;
  }
  //stab = stab_test(B_mat, N, P); //final stability test
  mat errors = Y - X * trans(B_mat);
  mat ortho_errors = errors * A_eps;
  return List::create(Named("beta") = B_sample, 
                      Named("B_mat") = B_mat,
                      Named("errors") = errors, 
                      Named("ortho_errors") = ortho_errors,
                      Named("attempts") = attempts,
                      Named("stability") = stab,
                      Named("b") = b11,
                      Named("p") = p11);//,   
                            //Named("beta_tt") = b11);
}

// [[Rcpp::export]]
cube comp_Hmat(const mat X, const int N) {
  int X_rows = X.n_rows;
  int X_cols = X.n_cols;
  cube H_cube(N, N * X_cols, X_rows);
  // cube H_cube_t(N * X_cols, N, X_rows);
  mat diag_mat = eye(N, N);
  for (int i = 0; i < X_rows; i++) {
    mat curr_element = kron(diag_mat, X.row(i));
    H_cube.slice(i) = curr_element;
    // H_cube_t.slice(i) = curr_element.t();
  }
  // return List::create(Named("H_cube") = H_cube, 
  //                     Named("H_cube_t") = H_cube_t);
  return H_cube;
} //WORKS


// [[Rcpp::export]]
List carterkohn_ca_cb_fast2(const mat Y, const mat X, mat Amat, const mat hlast, const vec b00, const mat p00, mat Bcomp, const cube H_cube, const int maxtry = 1e3) {
  // Timer timer;
  // timer.step("start"); 
  // TESTED: same b11 and p11 as carterkohn_ca_cb!!!!!!!!!
  int nobs = Y.n_rows;
  int N = Y.n_cols; //number of variables
  // int P = (size_b - N) / pow(N, 2); //number of lags for stability test and Bmat ( - N in numerator to correct for N intercepts in vec(b0))
  int P = (X.n_cols - 1) / N;
  int size = N * P; 
  int size_b = (size + 1) * N; //a_elemnts and hlast will form time-varying R
  int attempts; //for return only
  mat A_eps = trans(Amat); //needed for ortho errors; calculated here before inv() of A matrix
  vec b11 = b00; // starting values
  mat p11 = p00; // starting values
  // mat diag_mat = eye(N, N); //for creating Hmat with Kronecker Product
  Amat = inv(Amat);
  mat Amat_transp = trans(Amat);
  mat Y_transp = trans(Y); //don't have to transpose all the time in the loop
  // cube H_cube = H_list["H_cube"];
  // cube H_cube_t = H_list["H_cube_t"];
  // timer.step("prep");
  //Kalman filter
  for (int i = 0; i < nobs; i++) {
    // mat Hmat = kron(diag_mat, X.row(i)); //X expanded to allow for multiplying with vec(b)
    // mat H_transp = trans(Hmat);
    mat Hmat = H_cube.slice(i);
    // mat H_transp = H_cube_t.slice(i);
    mat H_transp = Hmat.t();
    //create R matrix: first create A, then R with hlast
    //SHOULD BE i+1 (in Matlab only "i")!!!!!
    mat R = Amat * diagmat(hlast.row(i + 1)) * Amat_transp; // i + 1 because of starting value of hlast which has to be ignored
    // Kalman steps
    // vec b10 = b11; //NEW
    // mat p10 = p11; //NEW
    vec eta10 = Y_transp.col(i) - Hmat * b11; //NEW should be b10 
    mat f10 = Hmat * p11 * H_transp + R; //NEW should be p10 rhs 
    mat Kgain = p11 * H_transp * inv(f10); //NEW should be p10 rhs
    b11 = b11 + Kgain * eta10; //NEW should be b10
    p11 = p11 - Kgain * Hmat * p11; //NEW should be p10 rhs
  }
  //timer.step("Kalman");
  // Carter Kohn step
  bool stab = false;
  mat B_mat;
  rowvec B_sample;
  cx_vec eigenval;
  //int tries = 1;
  //while (!chk) {
  B_sample = rmvnorm(1, b11, p11);
  B_mat = trans(reshape(B_sample, N * P + 1, N)); //for stability test and error/residual computation
  for (int i = 0; i < maxtry; i++) {
    // Rcout << "Inside Loop!!!";
    // int stab = 1; // stability result
    // Bcomp.rows(0, N - 1) = B_mat.submat(0, 0, N - 1, size - 1); //deleting the intercept: size - 1
    Bcomp.rows(0, N - 1) = B_mat.cols(0, size - 1);
    eigenval = eig_gen(Bcomp);
    if (max(abs(eigenval)) < 1) {
      // chk = true;
      stab = true;
      // attempts = i + 1;
      // Rcout << i;
      break;
    }
    B_sample = rmvnorm(1, b11, p11);
    B_mat = trans(reshape(B_sample, N * P + 1, N)); //for stability test and error/residual computation
  }
  //timer.step("loop");
  //stab = stab_test(B_mat, N, P); //final stability test
  mat errors = Y - X * B_mat.t();
  mat ortho_errors = errors * A_eps;
  return List::create(Named("beta") = B_sample, 
                      Named("B_mat") = B_mat,
                      Named("errors") = errors, 
                      Named("ortho_errors") = ortho_errors,
                      Named("stability") = stab);
                      // Named("b") = b11,
                      // Named("p") = p11);
                      //Named("timer") = timer);
}


// [[Rcpp::export]]
List carterkohn_ca_cb_fast(const mat Y, const mat X, mat Amat, const mat hlast, const vec b00, const mat p00, mat Bcomp, const int maxtry = 1e3) {
  int size_b = b00.size(); //a_elemnts and hlast will form time-varying R
  int nobs = Y.n_rows;
  int N = Y.n_cols; //number of variables
  // int P = (size_b - N) / pow(N, 2); //number of lags for stability test and Bmat ( - N in numerator to correct for N intercepts in vec(b0))
  int P = (X.n_cols - 1) / N;
  int size = N * P; //next 4 lines for stability testing
  int attempts; //for return only
  // int size_wo_lag1 = N * (P - 1);
  // mat Bcomp = zeros(size, size); 
  // Bcomp.submat(N, 0, size - 1, size_wo_lag1 - 1) = eye(size_wo_lag1, size_wo_lag1);
  mat A_eps = trans(Amat); //needed for ortho errors; calculated here before inv() of A matrix
  vec b11 = b00; // starting values
  mat p11 = p00; // starting values
  mat diag_mat = eye(N, N); //for creating Hmat with Kronecker Product
  Amat = inv(Amat);
  mat Amat_transp = trans(Amat);
  //Kalman filter
  for (int i = 0; i < nobs; i++) {
    mat Hmat = kron(diag_mat, X.row(i)); //X expanded to allow for multiplying with vec(b)
    mat H_transp = trans(Hmat);
    //create R matrix: first create A, then R with hlast
    //SHOULD BE i+1 (in Matlab only "i")!!!!!
    mat R = Amat * diagmat(hlast.row(i + 1)) * Amat_transp; // i + 1 because of starting value of hlast which has to be ignored
    // Kalman steps
    vec b10 = b11;
    mat p10 = p11;
    vec eta10 = trans(Y.row(i)) - Hmat * b10;
    mat f10 = Hmat * p10 * H_transp + R;
    mat Kgain = p10 * H_transp * inv(f10);
    b11 = b10 + Kgain * eta10;
    p11 = p10 - Kgain * Hmat * p10;
  }
  // Carter Kohn step
  bool stab = false;
  mat B_mat;
  rowvec B_sample;
  cx_vec eigenval;
  //int tries = 1;
  //while (!chk) {
  B_sample = rmvnorm(1, b11, p11);
  B_mat = trans(reshape(B_sample, N * P + 1, N)); //for stability test and error/residual computation
  for (int i = 0; i < maxtry; i++) {
    // Rcout << "Inside Loop!";
    // int stab = 1; // stability result
    // Bcomp.rows(0, N - 1) = B_mat.submat(0, 0, N - 1, size - 1); //deleting the intercept: size - 1
    Bcomp.rows(0, N - 1) = B_mat.cols(0, size - 1);
    eigenval = eig_gen(Bcomp);
    if (max(abs(eigenval)) < 1) {
      // chk = true;
      stab = true;
      attempts = i + 1;
      // Rcout << i;
      break;
    }
    B_sample = rmvnorm(1, b11, p11);
    B_mat = trans(reshape(B_sample, N * P + 1, N)); //for stability test and error/residual computation
  }
  //stab = stab_test(B_mat, N, P); //final stability test
  mat errors = Y - X * trans(B_mat);
  mat ortho_errors = errors * A_eps;
  return List::create(Named("beta") = B_sample, 
                      Named("B_mat") = B_mat,
                      Named("errors") = errors, 
                      Named("ortho_errors") = ortho_errors,
                      Named("attempts") = attempts,
                      Named("stability") = stab);
}

// [[Rcpp::export]]
List carterkohn_univ(mat Y, mat X, vec R, mat Q, vec b00, mat p00, vec mu, mat Fmat) {
  int nobs = Y.n_rows;
  int size_b = b00.size();
  mat b_fw(nobs, size_b); //Kalman filter
  cube p_fw(size_b, size_b, nobs);
  mat b_bw(nobs, size_b); // Carter Kohn return
  mat errors(nobs, 1);
  vec b11 = b00; //starting values
  mat p11 = p00;
  mat F_transp = trans(Fmat);
  for (int i = 0; i < nobs; i++) {
    vec b10 = mu + Fmat * b11;
    mat p10 = Fmat * p11 * F_transp + Q;
    vec eta10 = Y.row(i) - X.row(i) * b10;
    mat f10 = X.row(i) * p10 * trans(X.row(i)) + R(i + 1);
    mat Kgain = p10 * trans(X.row(i)) * inv(f10);
    b11 = b10 + Kgain * eta10;
    p11 = p10 - Kgain * X.row(i) * p10;
    b_fw.row(i) = trans(b11);
    p_fw.slice(i) = p11;
  }
  b_bw.row(nobs - 1) = rmvnorm(1, trans(b_fw.row(nobs - 1)), p_fw.slice(nobs - 1));
  errors.row(nobs - 1) = Y.row(nobs - 1) - X.row(nobs - 1) * trans(b_bw.row(nobs - 1));
  //Rcout <<  errors.row(nobs - 1) << std::endl;
  for (int i = nobs - 2; i >= 0; i--) {
    rowvec next_b = b_bw.row(i + 1);
    rowvec curr_b = b_fw.row(i);
    mat curr_p = p_fw.slice(i);
    vec eta_bw = trans(next_b) - mu - Fmat * trans(curr_b);
    mat f_bw = Fmat * curr_p * F_transp + Q;
    mat Kgain = curr_p * F_transp * inv(f_bw);
    vec b_update = trans(curr_b) + Kgain * eta_bw;
    mat p_update = curr_p - Kgain * Fmat * curr_p;
    b_bw.row(i) = rmvnorm(1, b_update, p_update);
    errors.row(i) = Y.row(i) - X.row(i) * trans(b_bw.row(i));
  }
  return List::create(Named("beta_bw") = b_bw, Named("errors") = errors);
}


// [[Rcpp::export]]
mat svol_mh(mat hlast_mat, mat eps_mat, vec g_vec, vec mubar_vec, double sigbar) { // for all variables -> first loop over variables
  int N = hlast_mat.n_cols; 
  int nobs = hlast_mat.n_rows; //nobs = "nobs + 1" because of starting value
  mat hnew_mat(nobs, N);
  for (int var_iter = 0; var_iter < N; var_iter++) { //iterating "h_i": one column at a time
    //prep: taking out the columns according to var_iter
    vec hlast = hlast_mat.col(var_iter);
    vec hnew(nobs); //results for this "var_iter"
    vec eps = eps_mat.col(var_iter);
    double g = g_vec(var_iter);
    double mubar = mubar_vec(var_iter);
    //first period: starting value
    int i = 0;
    double hlead = hlast(i + 1);
    double sig = sigbar * g / (sigbar + g);
    double mu = sig * (mubar / sigbar + log(hlead) / g);
    hnew(i) = exp(mu + pow(sig, .5) * R::rnorm(0, 1));
    //inbetween periods
    for (int i = 1; i < nobs - 1; i++) { //i=1 because i=0 was done before
      double hlead = hlast(i + 1);
      double curr_eps = eps(i - 1);
      double hlag = hnew(i - 1);
      double mu = (log(hlead) + log(hlag)) / 2;
      double sig = g / 2;
      double hprop = exp(mu + pow(sig, .5) * R::rnorm(0, 1));
      double pi_new = -.5 * log(hprop) - pow(curr_eps, 2) / (2 * hprop);
      double pi_old = -.5 * log(hlast(i)) - pow(curr_eps, 2) / (2 * hlast(i));
      if (exp(pi_new - pi_old) > R::runif(0, 1)) {
        hnew(i) = hprop;
      } else {
        hnew(i) = hlast(i);
      }
    }
    // last period
    i = nobs - 1;
    double curr_eps = eps(i - 1);
    mu = log(hnew(i - 1));
    sig = g;
    double hprop = exp(mu + pow(sig, .5) * R::rnorm(0, 1));
    //double hprop = exp(mu + pow(sig, .5) * rv(i));
    double pi_new = -.5 * log(hprop) - pow(curr_eps, 2) / (2 * hprop);
    double pi_old = -.5 * log(hlast(i)) - pow(curr_eps, 2) / (2 * hlast(i));
    if (exp(pi_new - pi_old) > R::runif(0, 1)) {
    //if (exp(pi_new - pi_old) > ru(i - 1)) {
      hnew(i) = hprop;
    } else {
      hnew(i) = hlast(i);
    }
    
    // save results for this iteration
    hnew_mat.col(var_iter) = hnew;
  }
  return hnew_mat;
}

// [[Rcpp::export]]
vec sample_g(mat hlast, int Tg0, int g0) {
  hlast.shed_row(0); //starting value deleted
  int nrow = hlast.n_rows;
  int N = hlast.n_cols;
  vec g(N);
  mat hlast_log = log(hlast);
  mat gerr = hlast_log.rows(1, nrow - 1) - hlast_log.rows(0, nrow - 2);
  mat vpv_mat = gerr.t() * gerr; //we only need the diagonal terms in loop
  // Rcout << nrow << std::endl << Tg0 << std::endl << vpv_mat(0,0) << std::endl << g0;
  for (int i = 0; i < N; i++) {
    //double vpv = as_scalar(dot(gerr.col(i), gerr.col(i)));
    // if (i == 0)
    // Rcout << (nrow - 1 + Tg0) / 2.0 << std::endl << (vpv_mat(i, i)  + g0) / 2.0;
    g(i) = 1 / as<double>(rgamma(1, (nrow - 1 + Tg0) / 2.0 , 2.0 / (vpv_mat(i, i) + g0))); //2 or 2.0: no diff
  } 
  return g;
}

/*
gerr <- diff(log(hlast)) #we have sampled ht but the state equation works with ln(ht)
#all.equal(gerr, mat_loop$all.g.errors) #works
  for (i in 1:N) {
    vpv <- as.numeric(crossprod(gerr[, i]))
    g[i] <- rinvgamma(1, (eff_nobs + Tg0) / 2, (vpv + g0) / 2) #excatley nobs + Tg0 because h has nobs + 1 obs and one obs is lost
  }
*/
 
// [[Rcpp::export]]
List sample_A(mat res, mat hlast, int V0) { //returns A not A^-1!; V0 is the var of the normal prior
  int N = res.n_cols; //residuals from VAR
  int nobs = res.n_rows;
  rowvec a_draw(N * (N - 1) / 2);
  int a_counter = 0;
  for (int i = 1; i < N; i++) {
    vec y_res = res.col(i);
    mat x_res = res.submat(0, 0, nobs - 1, i - 1);
    vec curr_h = sqrt(hlast.submat(1, i, nobs, i)); //nobs is correct because of additional row in h; starting value excluded
    mat curr_h_mat(nobs, i);
    for (int j = 0; j < i; j++) 
      curr_h_mat.col(j) = curr_h;
    y_res = y_res / curr_h;
    x_res = x_res / curr_h_mat * (-1);
    // Rcout << y_res << std::endl;
    mat x_res_trans = trans(x_res);
    mat Sigma0 = eye(i, i) * V0; //prior variance
    vec A0 = zeros(i); //prior mean hard coded; zero vector
    mat Vstar = inv(inv(Sigma0) + x_res_trans * x_res); // scalar for i = 1
    vec Mstar = Vstar * (inv(Sigma0) * A0 + x_res_trans * y_res);
    // Rcout << Vstar << std::endl << Mstar;
    rowvec curr_a_draw = rmvnorm(1, Mstar, Vstar);
    int end_pos = a_counter + curr_a_draw.size() - 1;
    a_draw.subvec(a_counter, end_pos) = curr_a_draw;
    a_counter = end_pos + 1;
  }
  //Rcout << a_draw;
  mat Amat = eye(N, N); //covariance terms will be lower triangular
  a_counter = 0;
  for (int i_row = 0; i_row < N; i_row++) { // creating Amatrix: lower triangular
    for (int i_col = 0; i_col < i_row; i_col++) {
      Amat(i_row, i_col) = a_draw(a_counter);
      a_counter++;
    }
  }
  return List::create(Named("a_draw") = a_draw, Named("Amat") = Amat);
}

// [[Rcpp::export]]
mat comp_eps_cpp(mat errors, mat a_elem_comb) { //for time varying A models
  int nobs = errors.n_rows;
  int N = errors.n_cols;
  mat eps(nobs, N); //ortho residuals: will be returned
  mat Amat = eye(N, N);
  for (int i = 0; i < nobs; i++) {
    int a_counter = 0;
    for (int i_row = 0; i_row < N; i_row++) { // creating Amatrix: lower triangular
      for (int i_col = 0; i_col < i_row; i_col++) {
        Amat(i_row, i_col) = a_elem_comb(i, a_counter);
        a_counter++;
      }
    }
    eps.row(i) = errors.row(i) * trans(Amat);
  }
  return eps;
}

// [[Rcpp::export]]
mat svol_ca_forecast_cpp(mat forecast_start, int P, int N, int forecast_horizon, mat B_mat, mat hlast, vec g, mat A) {
  vec v_zeros = zeros(N);
  mat y_hat = zeros(P + forecast_horizon, N);
  mat h_hat = zeros(P + forecast_horizon, N); //"too big" initalization makes it easier later on
  mat gmat = diagmat(g);
  colvec y_hat_multiplier(N * P + 1);
  y_hat_multiplier(N * P) = 1; // last element is the intercept
  y_hat.rows(0, P - 1) = forecast_start;
  h_hat.row(P - 1) = log(hlast.tail_rows(1)); // last row
  mat invA = inv(A);
  mat invA_trans = trans(invA);
  for (int i = P; i < forecast_horizon + P; i++) {
    //stochastic volatility first
    h_hat.row(i) = rmvnorm(1, trans(h_hat.row(i - 1)), gmat); 
    // Rcout << h_hat.row(i);
    mat Sigma = invA * diagmat(exp(h_hat.row(i))) * invA_trans;
    mat curr_y = y_hat.rows(i - P, i - 1);
    y_hat_multiplier.subvec(0, N * P - 1) = trans(vectorise(reverse(curr_y), 1));
    y_hat.row(i) = trans(B_mat * y_hat_multiplier) + rmvnorm(1, v_zeros, Sigma);
  }
  return y_hat;
}

// [[Rcpp::export]]
mat svol_ca_ar_forecast_cpp(mat forecast_start, int P, int N, int forecast_horizon, mat B_mat, mat hlast, vec curr_mu, vec curr_phi, vec curr_sig, mat A) {
  vec v_zeros = zeros(N);
  mat y_hat = zeros(P + forecast_horizon, N);
  mat h_hat = zeros(P + forecast_horizon, N); //"too big" initalization makes it easier later on
  // mat gmat = diagmat(g);
  colvec y_hat_multiplier(N * P + 1);
  y_hat_multiplier(N * P) = 1; // last element is the intercept
  y_hat.rows(0, P - 1) = forecast_start;
  h_hat.row(P - 1) = log(hlast.tail_rows(1)); // last row
  mat invA = inv(A);
  mat invA_trans = trans(invA);
  for (int i = P; i < forecast_horizon + P; i++) {
    //stochastic volatility first
    for (int j = 0; j < N; j++) {
      h_hat(i, j) = curr_mu(j) + curr_phi(j) * (h_hat(i - 1, j) - curr_mu(j)) + 
        R::rnorm(0, curr_sig(j));
    }
    // h_hat.row(i) = rmvnorm(1, trans(h_hat.row(i - 1)), gmat); 
    // Rcout << h_hat.row(i);
    mat Sigma = invA * diagmat(exp(h_hat.row(i))) * invA_trans;
    mat curr_y = y_hat.rows(i - P, i - 1);
    y_hat_multiplier.subvec(0, N * P - 1) = trans(vectorise(reverse(curr_y), 1));
    y_hat.row(i) = trans(B_mat * y_hat_multiplier) + rmvnorm(1, v_zeros, Sigma);
  }
  return y_hat;
}

// [[Rcpp::export]]
List test_restriction(mat Amat, mat restr) { //restriction holds the restrictions for each column for given row
  int nrow_A = Amat.n_rows;
  int n_restr = restr.n_rows;
  vec rest_col = restr.col(1); //
  vec rest_ineq = restr.col(2);
  for (int irow = 0; irow < nrow_A; irow++) {
    int counter = 0;
    int counter_neg = 0;
    for (int irest = 0; irest < n_restr; irest++) {
      if (Amat(irow, rest_col(irest) - 1) * rest_ineq(irest) > 0) //-1 correct?
        counter++;
      if (Amat(irow, rest_col(irest) - 1) * rest_ineq(irest) < 0)
        counter_neg++;
    }
    if (counter == n_restr)
      return List::create(Named("row") = irow, Named("mult") = 1);
    if (counter_neg == n_restr)
      return List::create(Named("row") = irow, Named("mult") = -1);
  }
  return List::create(Named("row") = 0, Named("mult") = 0); //no column found
}

// [[Rcpp::export]]
mat qr_decomp(mat K) {
  int size_Q = K.n_rows;
  mat Q, R;
  qr_econ(Q, R, K);
  for (int i = 0; i < size_Q; i++) {
    if (R(i, i) < 0)
      Q.col(i) = -Q.col(i);
  }
  return(Q);
}

/*
// [[Rcpp::export]]
arma::field<arma::cube> irf_return(int nobs, int horizon, int N, int samples) {
  arma::field<arma::cube> irf(samples);
  irf.fill(arma::cube(nobs, horizon, N, arma::fill::zeros));
  return irf;
}
*/

/*
cube sign_restr_irf(List est_res, mat restrictions, vec shock, int A_iter, int horizon = 40) {
  cube beta = est_res["out_beta"];
  cube Amat = est_res["out_A"];
  cube Hmat = est_res["out_h"];
  int nobs = beta.n_rows;
  int N = Amat.n_cols;
  int samples = Amat.n_slices;
  field irf_res = irf_return(nobs, horizon, N, samples);
  return irf_res;
}
 */

/*
 * for i=1:size(out1,1);
 
 for j=1:size(out1,2)
 
 H=diag(squeeze(out2(i,j,:)));
 a=squeeze(out4(i,j,:));
 A=chofac(N,a);
 sigma=inv(A)*H*inv(A)';  %covariance matrix
 %sign restrictions
 chck=-1;
 while chck<0
 K=randn(N,N);
 QQ=getQR(K);
 A0hat=chol(sigma);
 A0hat1=(QQ*A0hat);  %candidate draw
 for m=1:N
 %check signs in each row
 e1=A0hat1(m,1)<0;  %Response of Y
 e2=A0hat1(m,2)<0;  %Response of P
 e3=A0hat1(m,3)>0;  %Response of R
 
 if e1+e2+e3==3
 MP=A0hat1(m,:);
 chck=10;
 else
 %check signs but reverse them
 e1=-A0hat1(m,1)<0;  %Response of Y
 e2=-A0hat1(m,2)<0;  %Response of P
 e3=-A0hat1(m,3)>0;  %Response of R
 
 if e1+e2+e3==3
 MP=-A0hat1(m,:);
 chck=10;
 end
 end
 end
 end
 %re-shuffle rows of A0hat1 and insert MP in the first row
 A0x=[]; %will hold rows of A0hat1 not equal to MP
 for m=1:N
 ee=sum(abs(A0hat1(m,:))==abs(MP));
 if ee==0
 A0x=[A0x;A0hat1(m,:)];
 end
 end
 A0new=[A0x;MP]; %A0 matrix to be used for impulse response
 shock=[0 0 1]; 
 btemp=squeeze(out1(i,j,:));
 btemp=reshape(btemp,N*L+1,N);
 zz=irfsim(btemp,N,L,A0new,shock,horz+L);
 zz=zz./repmat(zz(1,3),horz,N);
 irfmat(i,j,:,:)=zz;
 end
 end
 */