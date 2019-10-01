#dims = over which dim to apply: (time_dim, variable_dim); not needed: sample_dim
#cols = which variables to apply to?

#what we should do: (1+g)^2 = (1+g1)(1+g2) => g = sqrt((1+g1)(1+g2)) - 1
#what we do: g = (g1+g2) / 2
cumulate_data <- function(ds, dims, cols) { 
  n_cols <- length(cols) 
  dim_ds <- dim(ds) #matrix or array?
  time_dim <- dims[1]
  var_dim <- dims[2]
  stopifnot(length(dim_ds) >= max(dims))
  if (!(time_dim %in% c(1, 2)))
    stop("Time has to be the first or second dimension.")
  if (!(var_dim %in% c(2, 3)) || var_dim == time_dim)
    stop("Columns have to be the second or third dimension.")
  nobs <- dim_ds[time_dim]
  if (length(dims) != 2) #time and variable dimension needed
    stop("2 dimensions have to be provided for aggregation: time and variable.")
  if (!(length(dim_ds) %in% c(2, 3))) #matrix or array has to be supplied
    stop("Function 'cumulate_data' only works with matrices and 3-dim arrays.")
  if (length(dim_ds) == 3) {
    apply_dim <- 1:3
    n_samples <- dim(ds)[apply_dim[!(apply_dim %in% dims)]] #find the sample dimension and number of samples
    apply_dim <- apply_dim[apply_dim != time_dim] #finding the dimensions for apply
    ds_cumul <- apply(ds[,, cols], apply_dim, cumsum) #changes over 2dim -> MARGIN = c(1, 3)
    div_block <- array(rep(1:nobs, n_cols * n_samples), dim = dim(ds_cumul))
    ds_cumul <- ds_cumul / div_block
    if (time_dim != 1) {#time not in first position? apply reorder -> we have to restore correct order
      ds_cumul <- R.utils::wrap.array(ds_cumul, list(2, 1, 3))
      ds[,, cols] <- ds_cumul
    } else {#time dim is first dim; samples ordered last
      ds[, cols, ] <- ds_cumul
    }
  } else {#2dim
    apply_dim <- 2  #changs over first dimension -> time dim
    ds_cumul <- apply(ds[, cols, drop = F], apply_dim, cumsum)
    div_block <- matrix(rep(1:nobs, n_cols), ncol = n_cols)
    ds_cumul <- ds_cumul / div_block
    ds[, cols] <- ds_cumul
  }
  return(ds)
}

#the Matlab version is equivalent to curr_mse <- colMeans(resids[1:i,, drop = F]^2) => averaged over time
#Matlab verion is also transposed
get_mse <- function(fc, obs, horizons, sqrt = T, cumulative = F) {
  resids <- as.matrix(fc - obs)
  #all_mse <- matrix(NA, ncol(resids), length(horizons))
  all_mse <- matrix(NA, length(horizons), ncol(resids), dimnames = 
                       list(paste0("h", horizons), colnames(obs)))
  write_iter <- 1
  for (i in horizons) {
    if (cumulative) {
      curr_mse <- colMeans(resids[1:i,, drop = F]^2)
    } else {
      curr_mse <- resids[i,, drop = F]^2
    }
    if (sqrt)
      curr_mse <- sqrt(curr_mse)
    #all_mse[, write_iter] <- curr_mse
    all_mse[write_iter, ] <- curr_mse
    write_iter <- write_iter + 1
  }
  return(all_mse)
}

density_eval <- function(fc, obs, variables, horizons, fn = c("crps_sample", "logs_sample", 
                                                              "es_sample", "vs_sample")) {
  fn_string <- match.arg(fn)
  fn <- match.fun(fn)
  fc <- fc[, horizons,, drop = F] #retain only those horizons we are going to use
  obs <- obs[horizons,, drop = F]
  if (fn_string %in% c("crps_sample", "logs_sample")) { #univariate measures => loop over variables
    writer_index <- 1
    score_res <- matrix(NA, length(horizons), length(variables),
                        dimnames = list(paste0("h", horizons),
                                        paste0("v", variables)))
    for (i_var in variables) {
      score_res[, writer_index] <- fn(obs[, i_var], t(fc[,, i_var])) #evaluates all horizons
      writer_index <- writer_index + 1
    }
    #score_res <- score_res[horizons, ]
  } else {
    score_res <- rep(NA, length(horizons))
    names(score_res) <- paste0("h", horizons)
    #multivariate scores work only with one time period to keep it simple => loop over time/horizon
    for (i_h in 1:length(horizons))  
      score_res[i_h] <- fn(obs[i_h, variables], t(fc[, i_h, variables]))
  }
  return(score_res)
}

#really use this one? 
get_rprob_density <- function(data, horizons, probs, variables = NULL) {
  n_variables <- dim(data)[3]
  if (is.null(variables)) #all variables evaluated
    variables <- 1:n_variables
  stopifnot(all(variables %in% 1:n_variables))
  #array of size probs_to_evlaute \times n_variables \times horizon
  rec_prob <- array(NA, dim = c(length(probs), n_variables, length(horizons)),
                     dimnames = list(paste0("p", probs), colnames(data), paste0("h", horizons)))
  write_horizon <- 1 #just for matrix indexing
  for (i_hor in horizons) {
    for (i_var in variables) {
      curr_ds <- data[, i_hor, i_var]
      density_var <- density(curr_ds)
      lower_bound <- min(curr_ds) - 2 * density_var$bw
      rec_prob[, i_var, write_horizon] <- try(sapply(probs, psample, density_var, lower_bound),
                                              silent = T)
    }
    write_horizon <- write_horizon + 1
  }
  return(rec_prob)
}

psample <- function(upper_b, sample_dens, lower_b) {
  return(integrate(approxfun(sample_dens), lower_b, upper_b)$value)
}

#recession prob based on samples: makes more sense
get_rprob_sample <- function(data, probs, horizons, variables = NULL) {
  if (is.null(variables)) #variables = NULL => use all variables
    variables <- 1:dim(data)[3]
  n_samples <- dim(data)[1] #should even work if I have to delete some samples because of NA
  rprob <- array(NA, dim = c(length(probs), length(variables), max(horizons)),
                 dimnames = list(paste0("p", probs), paste0("v", variables),
                                 paste0("h", 1:max(horizons))))
  for (i_var in 1:length(variables)) {
    for (i_prob in 1:length(probs)) {
      curr_variable <- variables[i_var]
      rprob[i_prob, i_var, ] <- apply(data[,, curr_variable, drop = F], 2, 
                                      function(x) sum(x < probs[i_prob])) #apply over horizons
    }
  }
  rprob <- rprob[,, horizons] #deleting unused horizons; nice side effect: simplifies to matrix if possible
  return(rprob / n_samples)
}

#used for PIT calc
get_pit_bin <- function(fc, obs, variables, horizons) {
  n_variables <- length(variables)
  #quant_seq <- seq(0, .9, .1)
  quant_seq <- seq(.1, .9, .1)
  res <- matrix(NA, n_variables, max(horizons), 
                dimnames = list(paste0("v", 1:n_variables),
                                paste0("h", 1:max(horizons))))
  for (i_var in 1:n_variables) {
    curr_variable <- variables[i_var] #find variable column
    curr_quan <- apply(fc[,, curr_variable], 2, quantile, quant_seq) #quantiles of dim 10 * max(horizon)
    for (i_hor in 1:max(horizons))  { #loop over horizons
      #obs falls outside obs -> too small by construction of quant_seq: 0 bin
      #curr_bin <- sum(obs[i_hor, curr_variable] >= curr_quan[, i_hor]) 
      #res[i_var, i_hor] <- ifelse(curr_bin == 0, 1, curr_bin)
      res[i_var, i_hor] <- sum(obs[i_hor, curr_variable] >= curr_quan[, i_hor]) + 1
    }
  }
  return(t(res)) #return all horizons until max
}

adjust_output <- function(arr_input) {#delets nas in columns and reorganizes
  temp_arr <- array(NA, dim = dim(arr_input), dimnames = dimnames(arr_input))
  if (length(dim(temp_arr)) == 3) {#it's an array
    temp_arr <- temp_arr[1:sum(!is.na(arr_input[, 1, 1])),, ]
    for (hor_index in 1:dim(temp_arr)[2]) {
      curr_mat <- arr_input[, hor_index, ]
      to_del <- unique(which(is.na(curr_mat), arr.ind = T)[, 1])
      temp_arr[, hor_index, ] <- curr_mat[-to_del, ]
    }
  } else {#it's a matrix
    temp_arr <- temp_arr[1:sum(!is.na(arr_input[, 1])), ]
    for (hor_index in 1:dim(temp_arr)[2]) {
      curr_mat <- arr_input[, hor_index]
      to_del <- unique(which(is.na(curr_mat)))
      temp_arr[, hor_index] <- curr_mat[-to_del]
    }
  }
  return(temp_arr)
}

# psample <- function(pr, sample_dens) { #area under the curve until pr
#   return(integrate(sample_dens, -Inf, pr)$value)
# }

comp_growthrate <- function(tseries, freq = 4, log = T) {
  freq <- freq * 100
  if (log) {
    return(freq * (log(tseries) - dplyr::lag(log(tseries))))
  } else {
    return(freq * (tseries - dplyr::lag(tseries)))
  }
}

aggregate_data <- function(x, freq = c("monthly", "quarterly", "yearly"), 
                           method = c("average", "sum", "max", "end_of_period"), verbose = T) {
  freq <- match.arg(freq)
  method <- match.arg(method)
  quart_flag <- FALSE #in case the data is already quarterly: avoid month = month - 2
  if (freq == "monthly") {
    agg_freq <- xts::apply.monthly
    xts_ident <- "months" #only used for end_of_period
  } else if (freq == "quarterly") {
    agg_freq <- xts::apply.quarterly
    xts_ident <- "quarters"
  } else {
    agg_freq <- xts::apply.yearly
    xts_ident <- "years"
  }
  if (method == "average") {
    agg_fn <- mean
  } else if (method == "sum") {
    agg_fn <- sum
  } else if (method == "max") {
    agg_fn <- max
  }
  if (verbose) {
    x_freq <- min(diff(x$date))
    if (x_freq <= 3) {
      cat("Daily data detected.\n")
    } else if (x_freq <= 7) {
      cat("Weekly data detected.\n")
    } else if (x_freq <= 31) {
      cat("Monthly data detcted.\n")
    } else if (x_freq <= 91) {
      cat("Quarterly data detected.\n")
      quart_flag <- TRUE
    } else {
      cat("Unknown data frequency.\n")
    }
  } 
  x_zoo <- zoo::zoo(x[, 2], x$date)
  if (method != "end_of_period") {
    x_agg <- agg_freq(x_zoo, agg_fn, na.rm = TRUE)
  } else {
    x_agg <- x_zoo[xts:::endof(x_zoo, xts_ident)]
  }    
  day(time(x_agg)) <- 1
  x_agg <- data.frame(as.Date(time(x_agg)), as.numeric(x_agg))
  colnames(x_agg) <- colnames(x)
  if (freq == "quarterly" && !quart_flag)
    month(x_agg$date) <- month(x_agg$date) - 2
  return(x_agg)
}

transform_fred <- function(fred_db, transf_data, columns = "all", scale = F, date_colum = T) {
  stopifnot(ncol(fred_db) == ncol(transf_data))
  N <- ncol(fred_db)
  #input checking and finding col_to_pick
  if (columns[1] == "all") {
    col_to_pick <- 1:N
  } else if (is.character(columns)) {
    col_to_pick <- numeric()
    for (i in columns) {
      pos_var <- which(colnames(transf_data) == i)
      if (length(pos_var) == 0) {
        warning(paste("Variable", i, "does not exist."))
      } else {
        col_to_pick <- c(col_to_pick, pos_var)
      }
    }
  } else if (is.numeric(columns)) {
    not_available <- columns[columns > N || columns < 0]
    if (length(not_available) != 0) 
      warning(paste("Columns", not_available, "are ignored."))
    col_to_pick <- columns[columns <= N && columns > 0]
  }
  stopifnot(length(col_to_pick) > 0)
  coln <- colnames(fred_db)[col_to_pick]
  #get data
  #reduced_ds <- fred_db[, col_to_pick]
  #transf_ds <- transf_data[, col_to_pick]
  final_db <- data.frame(row.names = 1:nrow(fred_db))
  for (i in col_to_pick) {
    t_code <- transf_data[, i]
    if (t_code == 1) {
      curr_col <- fred_db[, i]
    } else if (t_code == 2) {
      curr_col <- c(NA, diff(fred_db[, i]))
    } else if (t_code == 3) {
      curr_col <- c(NA, NA, diff(fred_db[, i], differences = 2))
    } else if (t_code == 4) {
      curr_col <- log(fred_db[, i])
    } else if (t_code == 5) {
      #curr_col <- log(fred_db[, i]) - log(dplyr::lag(fred_db[, i]))
      curr_col <- c(NA, diff(log(fred_db[, i])))
    } else if (t_code == 6) {
      # curr_col <- log(fred_db[, i]) - log(dplyr::lag(fred_db[, i]))
      # curr_col <- c(NA, diff(curr_col))
      curr_col <- c(NA, NA, diff(log(fred_db[, i]), differences = 2))
    } else {#case 7
      curr_col <- c(NA, diff(fred_db[, i] / dplyr::lag(fred_db[, i] - 1)))
    }
    final_db <- cbind(final_db, curr_col)
  }
  if (scale)
    final_db <- data.frame(scale(final_db))
  colnames(final_db) <- coln
  if (date_colum)
    final_db <- cbind(date = fred_db[, 1], final_db)
  return(final_db)
} 

transform_data <- function(db, transf_data, columns = "all", scale = F, date_colum = T) {
  stopifnot(ncol(db) == ncol(transf_data))
  N <- ncol(db)
  pos_date_col <- which(sapply(db, is.Date))
  if (length(pos_date_col) > 1) {
    warning("More than one date column provided. First date column is used.")
    pos_date_col <- pos_date_col[1]
  }
  #input checking and finding col_to_pick
  if (columns[1] == "all") {
    col_to_pick <- which(sapply(db, is.numeric))
    # transf_data <- col_to_pick
  } else if (is.character(columns)) {
    col_to_pick <- numeric()
    for (i in columns) {
      pos_var <- which(colnames(transf_data) == i)
      if (length(pos_var) == 0) {
        warning(paste("Variable", i, "does not exist."))
      } else {
        col_to_pick <- c(col_to_pick, pos_var)
      }
    }
  } else if (is.numeric(columns)) {
    not_available <- columns[columns > N || columns < 0]
    if (length(not_available) != 0) 
      warning(paste("Columns", not_available, "are ignored."))
    col_to_pick <- columns[columns <= N && columns > 0]
  }
  stopifnot(length(col_to_pick) > 0)
  coln <- colnames(db)[col_to_pick]
  #get data
  #reduced_ds <- db[, col_to_pick]
  #transf_ds <- transf_data[, col_to_pick]
  final_db <- data.frame(row.names = 1:nrow(db))
  for (i in col_to_pick) {
    t_code <- transf_data[, i]
    if (t_code == 1) {
      curr_col <- db[, i]
    } else if (t_code == 2) {
      curr_col <- c(NA, diff(db[, i]))
    } else if (t_code == 3) {
      curr_col <- c(NA, NA, diff(db[, i], differences = 2))
    } else if (t_code == 4) {
      curr_col <- log(db[, i])
    } else if (t_code == 5) {
      curr_col <- c(NA, diff(log(db[, i])))
    } else if (t_code == 6) {
      curr_col <- c(NA, NA, diff(log(db[, i]), differences = 2))
    } else {#case 7
      curr_col <- c(NA, diff(db[, i] / dplyr::lag(db[, i] - 1)))
    }
    final_db <- cbind(final_db, curr_col)
  }
  if (scale)
    final_db <- data.frame(scale(final_db))
  colnames(final_db) <- coln
  if (date_colum) {
    if (length(pos_date_col) == 0) {
      warning("No date column provided")
    } else {
      final_db <- cbind(date = db[, pos_date_col], final_db)
    }
  }
  return(final_db)
} 

#Recession dates
rec_date <- function(freq = c("monthly", "quarterly", "daily"), start_date = NULL, end_date = NULL) {
  freq <- match.arg(freq)
  if (freq == "monthly") {
    indic <- "USREC"
  } else if (freq == "quarterly") {
    indic <- "USRECQ"
  } else {
    indic <- "USRECD"
  }
  recession <- alfred::get_fred_series(indic, "US_REC")
  if (is.null(start_date))
    start_date <- min(recession$date)
  if (is.null(end_date))
    end_date <- max(recession$date)
  recession <- dplyr::filter(recession, date >= start_date & date <= end_date) %>%
    dplyr::mutate(first_diff = c(NA, diff(US_REC)), 
                  start = ifelse(first_diff == 1, 1, 0),
                  end = ifelse(first_diff == -1, 1, 0))
  rec_start <- which(recession$start == 1)
  rec_end <- which(recession$end == 1)
  if (rec_start[1] > rec_end[1]) #we start with recession
    rec_start <- c(1, rec_start)
  if (length(rec_start) > length(rec_end)) #we end with recession
    rec_end <- c(rec_end, nrow(recession))
  rec_db <- data.frame(start_date = recession$date[rec_start], 
                       end_date = recession$date[rec_end])
  # rec_db$end <- rec_db$end - 1
  return(rec_db)
}
