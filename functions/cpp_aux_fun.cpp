#include <RcppDist.h>
#include <Rcpp/Benchmark/Timer.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
// [[Rcpp::export]]
double lik_cpp(const mat x, const vec mu, const mat sigma, bool loglik = true) { //WORKS
  int nobs = x.n_rows; //number of rows -> observatiosn
  int N = x.n_cols; //number of cols -> number of variables
  double constant = -(nobs * N) / 2.0 * log(2.0 * M_PI);
  mat C_inv = inv(trimatu(chol(sigma)));
  double log_det_cov_sqrt = sum(log(C_inv.diag())); //sqrt of det logged // log det div by 2
  colvec x_i;
  rowvec x_C_inv;
  double quadform;
  
  double ll = 0; //loglik
  for (int i = 0; i < nobs; i++) {
    x_i = trans(x.row(i)) - mu;
    x_C_inv = trans(x_i) * C_inv;
    quadform = dot(x_C_inv, x_C_inv);
    ll += -0.5 * quadform;
  }
  ll = ll + nobs * log_det_cov_sqrt + constant;
  if(!loglik){
    ll = exp(ll);
  }
  return ll;
}

// [[Rcpp::export]]
mat riwish_cpp(const int v, const mat S) {
  return riwish(v, S);
}

// [[Rcpp::export]]
mat mvrnorm_cpp(const int n, const vec mu, const mat sigma, bool column_vec = false) { // needed for STVAR
  mat x = rmvnorm(n, mu, sigma);
  if (column_vec) {
    x = trans(x);
  }
  return x;
}

/*
// [[Rcpp::export]]
mat mvrnorm_simple(vec mu, mat sigma) { // needed for STVAR
  mat x = rmvnorm(1, mu, sigma);
  return x;
}

// [[Rcpp::export]]
mat mvrnorm_simple2(vec mu, mat sigma) { // needed for STVAR
  mat x = rmvnorm(1, mu, sigma);
  return x;
}

 
// [[Rcpp::export]]
mat mvrnorm_colvec(const vec mu, const mat S) {
  uword m = S.n_cols, j;
  rowvec result(m);
  //rowvec Mu = mu.t();
  for (j = 0; j < m; j++) {
     result(j) = R::rnorm(0.0, 1.0);
  }
  return trans(result * chol(S)) + mu;
}
*/

// [[Rcpp::export]]
mat xpx_inv_offset_cpp(const mat x, mat offset) { //WORKS
  return inv(trans(x) * x + offset);
}

// [[Rcpp::export]]
mat xpx_inv_cpp(const mat x) { //WORKS
  return inv(trans(x) * x);
}

// [[Rcpp::export]]
mat resids_cpp(const mat Y, const mat X, const mat B, bool crossprod = false) { //WORKS
  mat resid = Y - X * B.t(); //B mat is of dim N x (N*P) + 1 -> .t()
  if (crossprod) {
    resid = resid.t() * resid;
  }
  return resid;
}

// [[Rcpp::export]]
List sample_B_cpp(const vec B_star, const mat Sigma_sample, const mat xpx_inv, const int N, const int P, const int max_attempts) {
  int attempt = 1; //if no stability test -> max_attempts = 0 -> no loop -> no eigenvalue eval
  int size = N * P;
  int size_wo_lag1 = N * (P - 1);
  bool chk = false;
  mat Sig = kron(Sigma_sample, xpx_inv);
  cx_vec eigenval;
  mat Bcomp = zeros(size, size);
  colvec B_sample = trans(rmvnorm(1, B_star, Sig));
  mat B_mat = trans(reshape(B_sample, size + 1, N)); //size + 1 because of intercept
  Bcomp.rows(0, N - 1) = B_mat.submat(0, 0, N - 1, size - 1);
  Bcomp.submat(N, 0, size - 1, size_wo_lag1 - 1) = eye(size_wo_lag1, size_wo_lag1);
  while (attempt <= max_attempts) {
    eigenval = eig_gen(Bcomp);
    // Rcout << eigenval << std::endl;
    if (max(abs(eigenval)) < 1) {
      chk = true;
      break;
    } else {
      B_sample = trans(rmvnorm(1, B_star, Sig)); //nothing to gain in speed
      B_mat = trans(reshape(B_sample, size + 1, N));
      Bcomp.rows(0, N - 1) = B_mat.submat(0, 0, N - 1, size - 1);
      attempt++;
    }
  }
  return List::create(Named("B_sample") = B_sample,
                      Named("B_mat") = B_mat,
                      Named("chk") = chk,
                      Named("attempt") = attempt);
                      //Named("Eigenvalue") = eigenval);
}

// [[Rcpp::export]]
List sample_B_fast(const vec B_star, const mat Sigma_sample, const mat xpx_inv, mat Bcomp, const int size, const int N, const int P, const int max_attempts) {
  // Timer timer;
  // timer.step("start");
  int attempt = 1; //if no stability test -> max_attempts = 0 -> no loop -> no eigenvalue eval
  bool chk = false;
  mat Sig = kron(Sigma_sample, xpx_inv);
  cx_vec eigenval;
  //timer.step("sample1");
  colvec B_sample = trans(rmvnorm(1, B_star, Sig)); //first attempt/draw
  // colvec B_sample = mvrnorm_colvec(B_star, Sig);
  // timer.step("transf1");
  mat B_mat = trans(reshape(B_sample, size + 1, N)); //size + 1 because of intercept
  // timer.step("loop");
  while (attempt <= max_attempts) {
    Bcomp.rows(0, N - 1) = B_mat.submat(0, 0, N - 1, size - 1); //size - 1 because no intercept and indexing starts with zero
    eigenval = eig_gen(Bcomp);
    // Rcout << eigenval << std::endl;
    if (max(abs(eigenval)) < 1) {
      chk = true;
      break;
    } // next part should work wihout ELSE!!
    B_sample = trans(rmvnorm(1, B_star, Sig)); //nothing to gain in speed: rmvnorm from Dist package is fast
    // B_sample = mvrnorm_colvec(B_star, Sig);
    B_mat = trans(reshape(B_sample, size + 1, N));
    // Bcomp.rows(0, N - 1) = B_mat.submat(0, 0, N - 1, size - 1);
    attempt++;
  }
  // Rcout << attempt << std::endl;
  // timer.step("return_Val");
  return List::create(Named("B_sample") = B_sample,
                      Named("B_mat") = B_mat,
                      Named("chk") = chk);
                      // Named("attempt") = attempt,
                      // Named("timer") = timer);
  //Named("Eigenvalue") = eigenval);
}

// [[Rcpp::export]]
NumericVector eval_delay_thresh_cpp(const mat ZZ, const double tar_value, const mat YY, const mat XX, mat B1_sample_mat, mat B2_sample_mat, mat Sigma1_sample, mat Sigma2_sample) {
  int n_delays = ZZ.n_cols;
  vec zero_vec = zeros(YY.n_cols);
  NumericVector lh(n_delays);
  for (int i = 0; i < n_delays; i++) {
    vec ZZ_star = ZZ.col(i);
    mat Y1_new = YY.rows(find(ZZ_star <= tar_value));
    mat X1_new = XX.rows(find(ZZ_star <= tar_value));
    mat Y2_new = YY.rows(find(ZZ_star > tar_value));
    mat X2_new = XX.rows(find(ZZ_star > tar_value));
    mat resid1 = resids_cpp(Y1_new, X1_new, B1_sample_mat);
    mat resid2 = resids_cpp(Y2_new, X2_new, B2_sample_mat);
    lh(i) = lik_cpp(resid1, zero_vec, Sigma1_sample) + lik_cpp(resid2, zero_vec, Sigma2_sample);
  }
  return exp(lh - max(lh));
} //WORKS

// [[Rcpp::export]]
NumericVector sm_fun_cpp(const vec x, const vec param, const String sm_fun) {
  String fn = "logistic";
  if (sm_fun == fn) {
    vec res = pow(1 + exp(-exp(param(1)) * (x - param(0))), -1);
    return(NumericVector(res.begin(), res.end()));
    // return(wrap(pow(1 + exp(-exp(param(1)) * (x - param(0))), -1)));
  } else {
    vec res = 1 - exp(-exp(param(1))  * pow(x - param(0), 2));
    return(NumericVector(res.begin(), res.end()));
    // return(wrap(1 - exp(-exp(param(1))  * pow(x - param(0), 2))));
  }
}

// [[Rcpp::export]]
NumericVector eval_delay_smooth_cpp(const mat ZZ, vec mh_param, const mat YY, const mat XX, const mat B1_sample_mat, const mat B2_sample_mat, mat Sigma_sample, String sm_fun) {
  int n_delays = ZZ.n_cols;
  int n_rows_X = XX.n_rows;
  int n_cols_X = XX.n_cols;
  mat X1 = XX;
  mat X2(n_rows_X, n_cols_X);
  vec zero_vec = zeros(YY.n_cols);
  NumericVector lh(n_delays);
  for (int i = 0; i < n_delays; i++) {
    vec ZZ_star = ZZ.col(i);
    NumericVector st_value_new = sm_fun_cpp(ZZ_star, mh_param, sm_fun);
    for (int j = 0; j < n_rows_X; j++) {
      X2.row(j) = XX.row(j) * st_value_new(j);
    }
    mat resids = YY - X1 * trans(B1_sample_mat) - X2 * trans(B2_sample_mat);
    lh(i) = lik_cpp(resids, zero_vec, Sigma_sample);
  }
  return exp(lh - max(lh));
} //NOT SURE IF WORKS


// eval_delay_smooth <- function(delay, ZZ, mh_param, YY, XX, B1_sample_mat, B2_sample_mat,
//                               Sigma_sample, zero_vec, sm_fun) {
//   ZZ_star <- ZZ[, delay]
//   st_value_new <- sm_fun(ZZ_star, mh_param)
// #X1 <- XX * st_value_new; X2 <- XX * (1 - st_value_new)
//   X1 <- XX; X2 <- XX * st_value_new
//     resids <- YY - X1 %*% t(B1_sample_mat) - X2 %*% t(B2_sample_mat)
//     lh <- lik_cpp(resids, zero_vec, Sigma_sample, loglik = T)
//     return(lh)
// }

// [[Rcpp::export]]
mat lin_forecast_cpp(const mat forecast_start, const int P, const int N, const int forecast_horizon, const mat B_mat, const mat Sigma) {
  vec v_zeros = zeros(N);
  mat y_hat = zeros(P + forecast_horizon, N);
  colvec y_hat_multiplier(N * P + 1);
  y_hat_multiplier(N * P) = 1; // last element is the intercept
  y_hat.rows(0, P - 1) = forecast_start;
  for (int i = P; i < forecast_horizon + P; i++) {
    mat curr_y = y_hat.rows(i - P, i - 1);
    y_hat_multiplier.subvec(0, N * P - 1) = trans(vectorise(reverse(curr_y), 1));
    y_hat.row(i) = trans(B_mat * y_hat_multiplier) + rmvnorm(1, v_zeros, Sigma);
  }
  return y_hat;
}

// [[Rcpp::export]]
vec tar_fn_cpp(vec x, List tar_transform) {
  // Rcout << x << std::endl;
  String fn = tar_transform["fun"];
  int len_x = x.size();
  vec ret(len_x);
  if (fn == "growth") {
    int lag = tar_transform["inp"];
    vec x_lag(len_x);
    std::fill(x_lag.begin(), x_lag.end(), NumericVector::get_na());
    x_lag.subvec(lag, len_x - 1) = x.subvec(0, len_x - lag - 1);
    ret = log(x) - log(x_lag);
  } else if (fn == "ma" || fn == "sum") {
    int window = tar_transform["inp"];
    double summed = 0.0;
    for (int i = 0; i <= window - 1; i++) { //before first window values
      summed += x[i];
      ret[i] = NA_REAL;
    }
    ret[window - 1] = summed;
    for (int i = window; i < len_x; i++) {
      summed += x[i] - x[i - window];
      ret[i] = summed;
    }
    if (fn == "ma") {
      ret = ret / window;
    }
  } else { // ident function
    ret = x;
  }
  return ret;
} //WORKS

// [[Rcpp::export]]
mat tvar_forecast_cpp(mat forecast_start, int tar_variable, int P, int N, int forecast_horizon, int d_sample, int T_tardata, double tar_value, vec tardata, mat B1_mat, mat B2_mat, mat Sigma1, mat Sigma2, List tar_transform, mat data) {
  vec data_tar_col = data.col(tar_variable - 1); //tar data column
  // vec data_tar_col = data.col(tar_variable - 1); // tar data column
  int size_zz = P + forecast_horizon + d_sample;
  vec v_zeros = zeros(N); //for drawing random errors as mean
  mat y_hat = zeros(P + forecast_horizon, N);
  // colvec y_hat_multiplier = ones(N * P + 1);
  colvec y_hat_multiplier(N * P + 1);
  y_hat_multiplier(N * P) = 1; // last element is the intercept
  y_hat.rows(0, P - 1) = forecast_start;
  // Rcout << y_hat << std::endl;
  vec ZZ_elements(size_zz);
  std::fill(ZZ_elements.begin(), ZZ_elements.end(), NumericVector::get_na());
  ZZ_elements.subvec(0, P - 1) = zeros(P);
  // STILL TO CHECK
  ZZ_elements.subvec(P, P + d_sample - 1) = tardata.subvec(T_tardata - d_sample, T_tardata - 1);
  // Rcout << ZZ_elements << std::endl;
  int first_na = P + d_sample; //index for element which will be updated in loop
  // Rcout << first_na << std::endl;
  for (int i = P; i < (P + forecast_horizon); i++) {
    mat curr_y = y_hat.rows(i - P, i - 1);
    y_hat_multiplier.subvec(0, N * P - 1) = trans(vectorise(reverse(curr_y), 1));
    bool reg_indicator = ZZ_elements(i) > tar_value;
    // if (i == P || i == (P + 1) || i == (P + 2))
    //   Rcout << curr_y << std::endl << trans(y_hat_multiplier) << ZZ_elements(i) << std::endl << tar_value << std::endl << std::endl;
    // if (i == P || i == (P + 1))
    //   Rcout << ZZ_elements(i) << std::endl << tar_value << std::endl << reg_indicator;
    if (reg_indicator == 0) {
      y_hat.row(i) = trans(B1_mat * y_hat_multiplier) + rmvnorm(1, v_zeros, Sigma1);
      // Rcout << "reg low" << std::endl;
    } else {
      y_hat.row(i) = trans(B2_mat * y_hat_multiplier) + rmvnorm(1, v_zeros, Sigma2);
      // Rcout << "reg high" << std::endl;
    }
    vec tar_temp = join_vert(data_tar_col, y_hat.submat(P, tar_variable - 1, i, tar_variable - 1));
    tar_temp = tar_fn_cpp(tar_temp, tar_transform);
    ZZ_elements(first_na) = tar_temp(tar_temp.size() - 1); //first_na or first_na - 1 ?
    first_na++;
  }
  return y_hat;
}  // WORKS


// [[Rcpp::export]]
mat stvar_forecast_cpp(mat forecast_start, int tar_variable, int P, int N, int forecast_horizon, int d_sample, int T_tardata, vec param_mh, String sm_fun, vec tardata, mat B1_mat, mat B2_mat, mat Sigma, List tar_transform, mat data) {
  String fn = "logistic";
  double curr_st_value;
  vec data_tar_col = data.cols(tar_variable - 1, tar_variable - 1);
  int size_zz = P + forecast_horizon + d_sample;
  vec v_zeros = zeros(N);
  mat y_hat = zeros(P + forecast_horizon, N);
  // colvec y_hat_multiplier = ones(N * P + 1);
  colvec y_hat_multiplier(N * P + 1);
  y_hat_multiplier(N * P) = 1; // last element is the intercept
  // colvec y_hat_multiplier_smooth = y_hat_multiplier;
  y_hat.rows(0, P - 1) = forecast_start;
  vec ZZ_elements(size_zz);
  std::fill(ZZ_elements.begin(), ZZ_elements.end(), NumericVector::get_na());
  ZZ_elements.subvec(0, P - 1) = zeros(P);
  // STILL TO CHECK
  ZZ_elements.subvec(P, P + d_sample - 1) = tardata.subvec(T_tardata - d_sample, T_tardata - 1);
  int first_na = P + d_sample;
  for (int i = P; i < (P + forecast_horizon); i++) {
    if (sm_fun == fn) {
      curr_st_value = pow(1 + exp(-exp(param_mh(1)) * (ZZ_elements(i) - param_mh(0))), -1);
    } else {
      curr_st_value = 1 - exp(-exp(param_mh(1))  * pow(ZZ_elements(i) - param_mh(0), 2));
    }
    // NumericVector curr_st_value = sm_fun_cpp(ZZ_elements(i), param_mh, sm_fun);
    mat curr_y = y_hat.rows(i - P, i - 1);
    y_hat_multiplier.subvec(0, N * P - 1) = trans(vectorise(reverse(curr_y), 1));
    colvec y_hat_multiplier_smooth = y_hat_multiplier * curr_st_value;
    y_hat.row(i) = trans(B1_mat * y_hat_multiplier + B2_mat * y_hat_multiplier_smooth) + 
      rmvnorm(1, v_zeros, Sigma);
    vec tar_temp = join_vert(data_tar_col, y_hat.submat(P, tar_variable - 1, i, tar_variable - 1));
    tar_temp = tar_fn_cpp(tar_temp, tar_transform);
    ZZ_elements(first_na) = tar_temp(tar_temp.size() - 1); //first_na or first_na - 1 ?
    first_na++;
  }
  return y_hat;
} 

// y_hat <- matrix(NA, P + forecast_horizon, N)
//   y_hat[1:P, ] <- start_forecast
// #forecast  
//   ZZ_elements <- rep(NA, nrow(y_hat) + d_sample)
//   ZZ_elements[1:P] <- 0
// ZZ_elements[(P + 1):(P + d_sample)] <- tardata[(T_tardata - d_sample + 1):T_tardata, 1] #have another look if correct
//   first_na <- min(which(is.na(ZZ_elements)))
//   for (i in (P + 1):(P + forecast_horizon)) {
//     curr_st_value <- sm_fun(ZZ_elements[i], param_mh)
//     y_hat[i, ] <- (B1_sample_mat %*% c(c(t(y_hat[(i - 1):(i - P), ])), 1)) +
//       (B2_sample_mat %*% c(c(t(y_hat[(i - 1):(i - P), ])), 1)) * curr_st_value +
//       mvrnorm_cpp(1, zero_vec, Sigma_sample, column_vec = TRUE) 
// #has to change after transform
// #------------------------------
//     tar_temp <- tar_fn_cpp(c(data[, tar_variable], y_hat[(P + 1):i, tar_variable]), tar_transform)
//     ZZ_elements[first_na] <- tar_temp[length(tar_temp)]
//     first_na <- first_na + 1
// #------------------------------
//   }
//   out_yhat[iter - burn,, ] <- y_hat

// GIRF functions follow

// [[Rcpp::export]]
mat apply_girf_cpp(cube girf) { //WORKS
  int n_slices = girf.n_slices;
  int n_rows = girf.n_rows;
  int n_cols = girf.n_cols;
  mat girf_mean = zeros(n_rows, n_cols);
  for (int i = 0; i < n_slices; i++) {
    girf_mean += girf.slice(i);
  }
  girf_mean = girf_mean / n_slices;
  return girf_mean;
}

// [[Rcpp::export]]
mat struc_covmat_cpp(mat X, String type = "Chol") { //WORKS
  mat X_chol = trans(chol(X));
  String triang = "Triangular";
  if (triang == type) {
    mat gam_mat = diagmat(1 / X_chol.diag());
    X_chol = X_chol * gam_mat;
  }
  return X_chol;
}

// [[Rcpp::export]]
mat struc_error_cpp(mat Sigma, mat resids) { //WORKS
  return inv(Sigma) * trans(resids);
}

// [[Rcpp::export]]
mat simul_diff_cpp(cube curr_B_array, int N, vec history, mat str_err_wo, cube struc_Sigma, int horizon, int curr_reg, vec tardata, int curr_pos_in_ds, int P, int max_del, double curr_tar, int tar_variable, int shock_var, double shock_size, List tar_transform) {
  mat sim_res_wo = zeros(N, horizon);  //save simulation results
  mat sim_res_add = zeros(N, horizon);
  mat str_err_add = str_err_wo; //only structural errors without shock are passed to function
  str_err_add(shock_var - 1, 0) = shock_size; //creating structural errors with inital shock
  mat reg_struc_Sigma = struc_Sigma.slice((curr_reg - 1)); //cov matrix for transforming back to reduced form errors
  mat shock_wo = reg_struc_Sigma * str_err_wo; //building reduced form errors
  mat shock_add = reg_struc_Sigma * str_err_add;
  tar_variable = tar_variable - 1; //column pos of tar variable in ds
  int curr_reg_wo = curr_reg - 1; //evaluating which regime we are in during simulation
  int curr_reg_add = curr_reg - 1;
  mat curr_B_wo = curr_B_array.slice((curr_reg - 1)); //Betas for regimes
  mat curr_B_add = curr_B_wo;
  int tar_inp = tar_transform["inp"]; // for finding correct ZZ value in loop
  // Rcout << tar_inp << std::endl << std::endl;
  int len_history = history.size();  //N*P + 1
  vec history_wo = history; //histories for regimes
  vec history_add = history;
  vec temp_history = zeros(len_history); //used for updating histories
  int tard_counter = max_del - 1; //used for CPP indexing tard_used_wo and tard_used_add (therefore - 1)
  //int tard_counter = max_del;
  int len_tard_used = tard_counter + horizon + 1;
  vec tard_used_wo = zeros(len_tard_used);
  // CHECK AGAIN!!!!!!!!!!!!!!!!!!!!!
  tard_used_wo.subvec(0, tard_counter) = tardata.subvec((curr_pos_in_ds + P - max_del), 
                      (curr_pos_in_ds + P - 1));
  //tard_used_wo.subvec(0, tard_counter) = tardata.subvec((curr_pos_in_ds + P - max_del - 1), (curr_pos_in_ds + P - 1));
  // Rcout << trans(tard_used_wo) << std::endl << curr_pos_in_ds + P - max_del << std::endl;
  vec tard_used_add = tard_used_wo;
  vec ZZ_used_wo; //transformation of tar data
  vec ZZ_used_add;
  for (int i = 0; i < horizon; i++) {
    sim_res_wo.col(i) = curr_B_wo * history_wo + shock_wo.col(i);
    sim_res_add.col(i) = curr_B_add * history_add + shock_add.col(i);
    // if (i == 0 || i == 1)
    //   Rcout << sim_res_wo << std::endl << std::endl;
    //tard_counter += 1;
    tard_counter++;
    tard_used_wo[tard_counter] = sim_res_wo(tar_variable, i);
    tard_used_add[tard_counter] = sim_res_add(tar_variable, i);
    ZZ_used_wo = tar_fn_cpp(tard_used_wo.subvec(0, tard_counter), tar_transform);  //transform tardata
    ZZ_used_add = tar_fn_cpp(tard_used_add.subvec(0, tard_counter), tar_transform);
    curr_reg_wo = (ZZ_used_wo[i + tar_inp] > curr_tar);
    curr_reg_add = (ZZ_used_add[i + tar_inp] > curr_tar);
    // if (i == 0 || i == 1)
    //   Rcout << trans(tard_used_wo) << std::endl << trans(ZZ_used_wo) << std::endl << ZZ_used_wo[i + tar_inp] << std::endl << curr_tar << std::endl << curr_reg_wo << std::endl;
    //curr_reg_wo = (ZZ_used_wo[i + 1] > curr_tar);
    //curr_reg_add = (ZZ_used_add[i + 1] > curr_tar);
    curr_B_wo = curr_B_array.slice(curr_reg_wo);
    curr_B_add = curr_B_array.slice(curr_reg_add);
    temp_history.subvec(0, (N - 1)) = sim_res_wo.col(i); //update history_wo
    temp_history.subvec(N, (len_history - 2)) = history_wo.subvec(0, (len_history - N - 2)); //-2 because of intercept and indexing starting at zero
    temp_history.tail(1) = 1;
    history_wo = temp_history;
    temp_history.subvec(0, (N - 1)) = sim_res_add.col(i); //update history_add
    temp_history.subvec(N, (len_history - 2)) = history_add.subvec(0, (len_history - N - 2));
    temp_history.tail(1) = 1; //can be taken out
    history_add = temp_history;
    // if (i == 0 || i == 10 || i == 20)
    //   Rcout << curr_reg_wo << std::endl << curr_B_wo << std::endl << std::endl;
  }
  return trans(sim_res_add - sim_res_wo);
}

// EXPERIMENTS
//List sim_sample_cpp(List histories, List struc_errors, List pos_in_ds, int horizon, int N, int P, int max_delay_needed, int n_histories, int n_boot, cube curr_B_array, vec history, cube struc_Sigma, vec tardata, double curr_tar, int tar_variable, int shock_var, double shock_size, List tar_transform) {

// [[Rcpp::export]]
List sim_sample_cpp(List histories, List struc_errors, List pos_in_ds, int horizon, int N, int P, int max_delay_needed, cube curr_B_array, cube struc_Sigma, vec tardata, double curr_tar, int tar_variable, int shock_var, double shock_size, List tar_transform, int n_histories, int n_boot) {
  List final_girf_all_sample = List::create(Named("reg1") = mat(N, horizon),
                                            Named("reg2") = mat(N, horizon));
  mat curr_regime;
  mat regime_struc_errors;
  for (int reg_iter = 0; reg_iter < 2; reg_iter++) {
    curr_regime = as<mat>(histories[reg_iter]);
    regime_struc_errors = as<mat>(struc_errors[reg_iter]);
    int nobs_reg = curr_regime.n_rows;
    IntegerVector sample_seq = seq(1, nobs_reg);
    //Rcout << sample_seq << std::endl;
    cube girf_curr_reg(horizon, N, n_histories);
    for (int hist_iter = 0; hist_iter < n_histories; hist_iter++) {
      int curr_hist_pos = as<int>(sample(sample_seq, 1));
      vec curr_hist_pos_in_ds_temp = pos_in_ds[reg_iter];
      int curr_hist_pos_in_ds = curr_hist_pos_in_ds_temp[curr_hist_pos - 1];
      rowvec curr_hist = curr_regime.row(curr_hist_pos - 1);
      //Rcout << curr_hist_pos << std::endl << curr_hist_pos_in_ds << curr_hist << std::endl;
      cube girf_boot(horizon, N, n_boot);
      for (int boot_iter = 0; boot_iter < n_boot; boot_iter++) {
        //vec curr_error_pos = as<vec>(sample(sample_seq, horizon));
        uvec curr_error_pos = as<uvec>(sample(sample_seq, horizon)) - 1;
        //Rcout << curr_error_pos << std::endl << max(curr_error_pos) << std::endl << nobs_reg;
        //mat curr_str_errors(N, horizon);
        mat curr_str_errors = regime_struc_errors.cols(curr_error_pos);
        // for (int i = 0; i < curr_error_pos.size(); i++ ) {
        //   int curr_index = curr_error_pos(i);
        //   curr_str_errors.col(i) = regime_struc_errors.col(curr_index - 1);
        // } 
        // Rcout << curr_str_errors << std::endl;
        //reg_iter + 1 for the moment to make comparison with R code possible, we could pass reg_iter and change code in slice in simul_cpp
        girf_boot.slice(boot_iter) = simul_diff_cpp(curr_B_array, N, trans(curr_hist), curr_str_errors, 
                        struc_Sigma, horizon, reg_iter + 1, tardata, curr_hist_pos_in_ds, P, max_delay_needed, 
                        curr_tar, tar_variable, shock_var, shock_size, tar_transform);
        //if (boot_iter == 0)
        //  Rcout << girf_boot.slice(0) << std::endl;
        //Rcout << boot_iter << std::endl;
      }
      girf_curr_reg.slice(hist_iter) = apply_girf_cpp(girf_boot);
    }
    //mat temp_girf = apply_girf_cpp(girf_curr_reg);
    final_girf_all_sample[reg_iter] = apply_girf_cpp(girf_curr_reg);
  }
  return final_girf_all_sample;
}

/*
// [[Rcpp::export]]
mat tt(mat x) {
  IntegerVector sample_seq = seq(1, x.n_cols);
  IntegerVector index = sample(sample_seq, 10);
  Rcout << index << std::endl;
  uvec index1 = as<uvec>(index) - 1;
  Rcout << index1 << std::endl;
  //mat res(x.n_rows, 10);
  mat res = x.cols(index1);
  return res;
}
*/