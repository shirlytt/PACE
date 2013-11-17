require(caTools)
require(pracma)

# options(warn = 2)
options(error = recover)

# source local weighted least square functions
source('lwls.R')

# done
isregular = function(t) {
  ans = 'sparse'
  tt = unlist(t)
  f = length(tt) / length(unique(tt)) / length(t)
  if (f > 0.75) {
    ans = 'missing'
  }
  if (f == 1) {
    ans = 'regular'
  }
  ans
}

# need discuss
bin_number = function(t) {
  mm = max(sapply(t, length))
  m = median(sapply(t, length))
  n = length(t)
  # for test
  # if (mm > 50) {
  #   return(50)
  # }
  # else {
  #   return(20)
  # }
  if (n <= 5000) {
    if (mm <= 400 | m <= 20) {
      # ans = mm
      ans = 50
    }
    else {
      ans = 400
    }
  }
  else {
    m_star = max(20, round(5000 - n) * 19 / 2250 + 400)
    if (m > m_star) {
      ans = m_star
    }
    else {
      # ans = mm
      ans = 50
    }
  }
}

# will port to c for sparse and missing
reformat = function(t, y, regular) {
  if (regular == 'sparse') {
    num_bins = bin_number(t)
    ans = bin_sparse(t, y, num_bins)
  }
  else {
    tt = unique(unlist(t))
    if (regular == 'regular') {
      ans = matrix(unlist(y), nrow = length(y), ncol = length(tt), byrow = TRUE)
    }
    else {
      ans = matrix(NA, nrow = length(y), ncol = length(tt), byrow = TRUE)
      for (i in 1:length(y)) {
        ans[i, match(t[[i]], tt)] = y[[i]]
      }
    }
    colnames(ans) = tt
  }
  ans
}

bin_sparse = function(t, y, num_bins) {
  tt = unlist(t)
  yy = unlist(y)
  n = length(yy)
  # index for subject
  # -1 for easy to use in C
  ii = unlist(lapply(1:length(y), function(i) rep(i, length(y[[i]])))) - 1
  mint = min(tt)
  maxt = max(tt)
  grids = seq(mint, maxt, length.out = num_bins + 1)
  d = grids[2] - grids[1]
  data = matrix(NA, ncol = num_bins, nrow = length(y))
  subject_sum = rep(0, num_bins)
  subject_count = rep(0, num_bins)
  # will rewrite in C
  for (i in 1:length(y)) {
    index = as.integer(cut(t[[i]], grids, include.lowest = TRUE))
    junk = tapply(y[[i]], index, mean)
    data[i, match(names(junk), 1:num_bins)] = junk
  }
  colnames(data) = (grids[1:(length(grids)-1)] + grids[2:length(grids)]) / 2
  data
  # has problem, so just don't do sparse for now, however it will
  # convert to missing case
  # ans = .C('bin_sparse', tt = as.double(tt), yy = as.double(yy), ii = as.integer(ii),
  #          num_bins = as.integer(num_bins), grids = as.double(grids), d = as.double(d),
  #          data = as.double(data), n = as.integer(n),
  #          subject_sum = as.double(subject_sum), subject_count = as.double(subject_count),
  #          mint = as.double(mint), NAOK = TRUE)$data
  # rownames(ans) = grids
  # t(ans)
}

# no of grids for cov-matrix won't be too large, just keep this part in R
get_raw_cov = function(data, mu, regular) {
  ans = list()
  n = nrow(data)
  p = ncol(data)
  data = as.matrix(sweep(data, 2, mu, '-'))
  if (regular == 'regular') {
    ans$raw_cov = t(data) %*% data / n
    ans$w2_in = matrix(n, ncol = ncol(data), nrow = ncol(data))
  }
  else {
    ans$raw_cov = matrix(NA, p, p)
    ans$w2_in = matrix(NA, p, p)
    for (i in 1:p) {
      for (j in i:p) {
        ans$raw_cov[i,j] = ans$raw_cov[j,i] = mean(data[,i] * data[,j], na.rm = TRUE)
        ans$w2_in[i,j] = ans$w2_in[j,i] = sum(!(is.na(data[,i]) | is.na(data[,j])))
      }
    }
    colnames(ans$raw_cov) = colnames(data)
    rownames(ans$raw_cov) = colnames(data)
    colnames(ans$w2_in) = colnames(data)
    rownames(ans$w2_in) = colnames(data)
  }
  ans
}

sigma_est = function(raw_cov, w2_in, bandwidth, kernel, userrange) {
  # x_in
  x_in = as.double(colnames(raw_cov))
  if (!is.null(userrange)) {
    junk_index = x_in >= userrange[1] & x_in <= userrange[2]
    x_in = x_in[junk_index]
    raw_cov = raw_cov[junk_index, junk_index]
    w2_in = w2_in[junk_index, junk_index]
  }
  # y_in
  var_y = as.double(diag(raw_cov))
  # variance of y
  var_y = lwls(bandwidth, kernel, x_in, var_y)$output

  # variance by x
  raw_cov_x_only = raw_cov
  diag(raw_cov_x_only) = NA
  w2_in[is.na(raw_cov_x_only) | is.nan(raw_cov_x_only)] = 0
  # observably need rewrite
  adhoc_raw_cov_only_x = adhoc_regular_cov(raw_cov_x_only, w2_in)
  var_x = diag(lwls_2d_vanilla(c(bandwidth[1], bandwidth[2]), kernel, adhoc_raw_cov_only_x))
  # lwls C version need to handle nan case

  # variance by yself
  sigma = trapz(x_in, var_y - var_x) / (max(x_in) - min(x_in))

  if (sigma < 0) {
    cat('Estimated sigma is negative, reset to zero now!\n')
  }

  max(0, sigma)
}

pcs_est = function(data, mu, hash, unique_hash, lambda,
                   eigen_func_data, cov_surface_fitted, sigma_old,
                   # switch
                   # no switch for xi_est and y_pred since those are always needed
                   sigma_new_switch = FALSE, xi_var_switch = FALSE) {
  sigma_new = 0
  xi_est = matrix(NA, ncol = ncol(eigen_func_data), nrow = nrow(data))
  y_pred = matrix(NA, ncol = ncol(data), nrow = nrow(data))
  xi_var = list()
  # hash table acceleration does not work for very different missing pattern
  for (hash_iterator in unique_hash) {
    index = hash == hash_iterator
    flag = unlist(strsplit(hash_iterator, " ")) == 1
    part_cov_surface_fitted = cov_surface_fitted[flag, flag] + diag(sigma_old, sum(flag))
    p = ncol(part_cov_surface_fitted)
    if (sum(flag) > 1) {
      junk = .Internal(La_svd('A', part_cov_surface_fitted, double(p),
                              matrix(double(p^2),p), matrix(double(p^2),p)))
      # eval this middle_matrix since it's need for both xi_est and xi_var
      middle_matrix = diag(lambda) %*% t(eigen_func_data[flag,]) %*% t(junk$vt) %*% diag(1 / junk$d) %*% t(junk$u)
      if (sum(index) > 1) {
        junk_y = as.matrix(sweep(data[index,], 2, mu, '-')[,flag])
      }
      else {
        junk_y = matrix((data[index,] - mu)[flag], nrow = 1)
      }
      xi_est[index,] = t(middle_matrix %*% t(junk_y))
    }
    else {
      junk = 1 / part_cov_surface_fitted[1,1]
      middle_matrix = diag(lambda) %*% matrix(eigen_func_data[flag,], ncol = 1) * junk
      if (sum(index) > 1) {
        junk_y = as.matrix(sweep(data[index,], 2, mu, '-')[,flag])
      }
      else {
        junk_y = matrix((data[index,] - mu)[flag], nrow = 1)
      }
      xi_est[index,] = sapply(junk_y, function(junk_yy) middle_matrix * junk_yy)
    }
    # k: number of components
    # p: grids of output without missing (=sum(flag))

    # diag(lambda):                                  k *          k
    # t(eigen_func_data[flag,]):                     k *          p
    # pinv:                                          p *          p
    # t(junk_y):                                     p * sum(index)
    # xi_est[index,]:                                sum(index) * k

    # be careful here
    if (sum(flag) > 1) {
      y_pred[index, flag] = xi_est[index,] %*% t(eigen_func_data[flag,])
    }
    else {
      y_pred[index, flag] = xi_est[index,] %*% matrix(eigen_func_data[flag,], ncol = 1)
    }
    # xi_est[index,]:                     sum(index) * k
    # t(eigen_func_data[flag,]):          k          * p
    # y_pred[index, flag]:                sum(index) * p
    if (sigma_new_switch) {
      sigma_new = sigma_new + sum(apply(junk_y - y_pred[index, flag], 1, function(x) mean(x^2)))
    }
    if (xi_var_switch) {
      for (ii in which(index)) {
        if (sum(flag) > 1) {
          xi_var[[ii]] = diag(lambda) - middle_matrix %*% t(diag(lambda) %*% t(eigen_func_data[flag,]))
        }
        else {
          xi_var[[ii]] = diag(lambda) - middle_matrix %*% t(diag(lambda) %*% matrix(eigen_func_data[flag,], ncol = 1))
        }
      }
    }
  }
  if (!sigma_new_switch) {
    sigma_new = NA
  }
  if (!xi_var_switch) {
    xi_var = NA
  }
  list(xi_est = xi_est,
       y_pred = y_pred,
       sigma_new = sigma_new / nrow(data),
       xi_var = xi_var)
}

cv_rho = function(data, mu_data, lambda, eigen_func_data, cov_surface_fitted, alpha_range = c(0.01, 0.225)) {
  clean = function(one, full) {
    full[full != one]
  }
  all_time = as.numeric(colnames(data))
  time_range = range(all_time)
  gamma = sqrt(trapz(all_time, mu_data^2) + sum(lambda)) / sqrt(time_range[2] - time_range[1])
  exist_matrix = !is.na(data)
  exist_list = apply(exist_matrix, 1, which)
  if (class(exist_list) == "list") {
    more_than_one = sapply(exist_list, length) > 1
  }
  else {
    more_than_one = apply(exist_matrix, 1, sum) > 1
  }
  data = data[more_than_one,]
  exist_list = exist_list[more_than_one]
  if (class(exist_list) == "list") {
    leave_one = sapply(exist_list, sample ,1)
  }
  else {
    leave_one = apply(exist_matrix, 1, function(x) sample(which(x), 1))
  }
  # keep = mapply(clean, leave_one, exist_list)
  y_true = rep(0, nrow(data))
  data_cv = data
  # may improve in case n is too large
  for (i in 1:nrow(data_cv)) {
    y_true[i] = data_cv[i, leave_one[i]]
    data_cv[i, leave_one[i]] = NA
  }
  #
  hash = apply(data_cv, 1, function(x) do.call(paste, as.list(as.integer(!is.na(x)))))
  unique_hash = unique(hash)
  ad_hoc_cv_func = function(rho) {
    xi_est = pcs_est(data_cv, mu_data, hash, unique_hash, lambda, eigen_func_data, cov_surface_fitted, rho)$xi_est
    pred_y = mu_data[leave_one]
    for (i in 1:length(lambda)) {
      pred_y = pred_y + xi_est * eigen_func_data[leave_one, i]
    }
    sum((y_true - pred_y)^2)
  }
  ans = optimize(ad_hoc_cv_func, alpha_range * gamma)$minimum
  cat("rho for ridge regression is:", ans, "\n")
  ans
}

FPCA = function(
  # raw data input format: 2 choice
  # 1:
  # data: a data.frame, ith row is the vector of measurements
  #       for the ith subject, i=1,...,n.
  #       colnames is the value for time points.
  #       (default)
  # 2:
  # y: a list, y[[i]] is the vector of measurements for the
  #    ith subject, i=1,...,n.
  # x: a list, y[[i]] is the vector of time points for the
  #    ith subject, i=1,...,n.
  data, y, x,
  # mean related
  bandwidth_mean, bandwidth_mean_cv_mode = c("geo_mean", "gcv", "cv"),
  mu_data = NULL, presmooth = TRUE,
  # covariance related
  bandwidth_cov, bandwidth_cov_cv_mode = c("gcv", "geo_mean", "cv"),
  # choosing number of components
  no_of_components = -1,
  # method of choosing
  selection = c("FVE", "BIC", "AIC"), FVE_threshold = 0.95,
  # maximum components
  maxk = 20,
  # measurements error or not
  error = TRUE,
  # grid of smoothed covariance surface
  cov_smooth_grid = 50,
  # used for computing random effects \xi_{ik}
  # 'CE': conditional expectation method
  # 'IN': classical integration method
  method_pcs = c("CE", "IN"),
  # applying shrinkage to estimates of random
  # coefficients (for regular data only)
  shrink = FALSE,
  # kernel
  kernel = c("Auto", "Epanechnikov", "Rectangular", "Gaussian", "Quartic", "Gausvar"),

  # number of bins
  # Don't let user set this part? It's just for speed up.
  # num_bins = NULL,

  # truncation step agreement
  # use cross validation for rho or user define its value
  # using new method in FPCscore.pdf to obtain ridge
  new_ridge = TRUE,
  rho_cv = FALSE, rho = 0,
  # method to estimate mu
  method_mu = c("PACE", "RARE"),
  # abandon rate
  out_percent = NULL, userrange = NULL) {
  ## if vanilla is true, use 2d smooth loop in R
  vanilla = TRUE
  # check format 2 data validity
  if (missing(data)) {
    # it's ok since we have missing or other cases
    # if (sum(sapply(y, function(yy) any(is.na(yy))))) {
    #   stop('FPCA is aborted because y contain NA(s)!')
    # }
    if (sum(sapply(x, function(xx) any(is.na(xx))))) {
      stop('FPCA is aborted because x contain NA(s)!')
    }
    if ((length(y) != length(x)) | (any(sapply(y, length) != sapply(x, length)))){
      stop('FPCA is aborted because y and x don\'t match!')
    }
    ni = sapply(y, length)
    if (all(ni == 1)) {
      stop('FPCA is aborted because the data do not contain repeated measurements!')
    }
    rm(ni)
  }
  kernel = match.arg(kernel)
  bandwidth_mean_cv_mode = match.arg(bandwidth_mean_cv_mode)
  bandwidth_mean_cv_mode = switch(bandwidth_mean_cv_mode,
                                  "cv" = 1,
                                  "gcv" = 2,
                                  "geo_mean" = 3)

  bandwidth_cov_cv_mode = match.arg(bandwidth_cov_cv_mode)
  bandwidth_cov_cv_mode = switch(bandwidth_cov_cv_mode,
                                 "cv" = 1,
                                 "gcv" = 2,
                                 "geo_mean" = 3)

  # These three should move to part 3
  # However check them first to avoid unexpected exit in future and weast time
  selection = match.arg(selection)
  # method_pcs = match.arg(method_pcs)
  method_mu = match.arg(method_mu)

  sigma = NULL
  mu_cov = NULL

  # if (missing(method_pcs) && shrink) {
  #   shrink = FALSE
  #   cat('Shrinkage method is avaiable when method_pcs is "IN" and error is TRUE!\n',
  #       'Reset to shrinkage method = FALSE!\n', sep = "")
  # }

  # check and reset number of components
  if (!missing(no_of_components)) {
    if (is.integer(no_of_components)) {
      if (no_of_components > (cov_smooth_grid - 2)) {
        no_of_components = (cov_smooth_grid - 2)
        cat('no_of_components can only be less than or equal to cov_smooth_grid - 2!\n',
            'Reset it to be cov_smooth_grid - 2 now!\n', sep = "")
      }
      if (no_of_components < 0) {
        no_of_components = -1
        cat('no_of_components must be a positive integer!\n',
            'Will choose automatically (Default method: FVE >= 0.95)!\n', sep = "")
      }
    }
    else {
      no_of_components = -1
      cat('no_of_components input is unsupported!\n',
          'Will choose automatically (Default method: FVE >= 0.95)!\n', sep = "")
    }
  }

  # check regular and reformat y and x input into data.frame
  sparse_flag = FALSE
  if (missing(data)) {
    regular = isregular(x)
    # the real problem is not whether sparse but procedure irregular
    # however, as long as data is irregular, previous procedure
    # will judge it as sparse and do binning

    # not clean here
    # however not a big deal
    if (regular == 'sparse' & presmooth) {
      sparse_flag = TRUE
      if (is.null(mu_data)) {
        junkx = as.double(unlist(x))
        junky = as.double(unlist(y))
        junk_index = order(junkx)
        junkx = junkx[junk_index]
        junky = junky[junk_index]
        num_bins = bin_number(junkx)
        minx = min(junkx)
        maxx = max(junkx)
        grids = seq(minx, maxx, length.out = num_bins + 1)
        x_data = (grids[1:(length(grids)-1)] + grids[2:length(grids)]) / 2
        x_cov = seq(min(x_data), max(x_data), length.out = cov_smooth_grid + 1)
        x_cov = (x_cov[1:cov_smooth_grid] + x_cov[2:(cov_smooth_grid+1)]) / 2
        kernel = switch(kernel,
                        "Epanechnikov" = 0,
                        "Rectangular" = 1,
                        "Gaussian" = 2,
                        "Quartic" = 3,
                        "Gausvar" = 4,
                        2)

        # add kernel choose here
        if (missing(bandwidth_mean)) {
          bandwidth_mean = bandwidth_choice_1d(kernel = as.integer(kernel),
                                               x_in = junkx, y_in = junky, w_in = rep(1, length(junkx)),
                                               cv_mode = bandwidth_mean_cv_mode)
        }
        # cat('Bandwidth for mean function is: ', bandwidth_mean, '.\n', sep = "")
        mu_data = lwls(bandwidth_mean, kernel = as.integer(kernel), x_in = junkx, y_in = junky,
                       count_in = rep(1, length(junkx)), x_out = x_data)$output
        mu_cov = lwls(bandwidth_mean, kernel = as.integer(kernel), x_in = junkx, y_in = junky,
                      count_in = rep(1, length(junkx)), x_out = x_cov)$output
        rm(junkx, junky, junk_index, num_bins, minx, maxx, grids)
      }
    }
    # for sparse data, do binning same time
    data = reformat(x, y, regular)
  }
  # from here we no longer have sparse case
  if (any(is.na(data))) {
    regular = 'missing'
  }
  else {
    regular = 'regular'
  }

  if (sparse_flag) {
    cat('The design is ', regular, ' (convert from sparse).\n', sep = "")
  }
  else {
    cat('The design is ', regular, '.\n', sep = "")
  }

  # if need, add binning part here

  if (kernel == "Auto") {
    if(regular == 'regular' & ncol(data) >= 20) {
      kernel = "Epanechnikov"
    }
    else {
      kernel = "Gaussian"
    }
  }
  kernel = switch(kernel,
                  "Epanechnikov" = 0,
                  "Rectangular" = 1,
                  "Gaussian" = 2,
                  "Quartic" = 3,
                  "Gausvar" = 4)

  # Part I: Obtain smoothed mean curve.
  cat('Part I: Obtain smoothed mean curve.\n')
  x_data = as.double(colnames(data))
  x_cov = as.double(seq(min(x_data), max(x_data),
                        length.out = cov_smooth_grid + 1))
  x_cov = (x_cov[1:(length(x_cov)-1)] + x_cov[2:length(x_cov)]) / 2

  y_data = as.double(colMeans(data, na.rm = TRUE))
  w_data = as.double(colSums(!is.na(data)))

  if (missing(bandwidth_mean)) {
    bandwidth_mean = bandwidth_choice_1d(kernel = as.integer(kernel),
                                         x_in = x_data, y_in = y_data, w_in = w_data,
                                         cv_mode = bandwidth_mean_cv_mode)
  }
  cat('Bandwidth for mean function is: ', bandwidth_mean, '.\n', sep = "")
  if (is.null(mu_data)) {
    mu_data = lwls(bandwidth_mean, kernel, x_data, y_data, w_data, x_out = x_data)$output
  }
  if (is.null(mu_cov)) {
    mu_cov = lwls(bandwidth_mean, kernel, x_data, y_data, w_data, x_out = x_cov)$output
  }
  rm(y_data, w_data)

  # Part II: Obtain smoothed covariance surface.
  cat('Part II: Obtain smoothed covariance surface.\n')
  # add weight
  junk = get_raw_cov(data, mu_data, regular)
  raw_cov = junk$raw_cov
  w2_cov = junk$w2_in
  rm(junk)

  if (vanilla) {
    adhoc_raw_cov = adhoc_regular_cov(raw_cov, w2_cov)
    diag(adhoc_raw_cov$newy) = NA
    diag(adhoc_raw_cov$count) = 0
  }

  # need fix after replace lwls_2d_vanilla()
  if (missing(bandwidth_cov)) {
    # bandwidth_cov = bandwidth_choice_2d(kernel, raw_cov, w2_cov,
    #                                     cv_mode = bandwidth_cov_cv_mode,
    #                                     symmetric = TRUE)
    if (vanilla) {
      bandwidth_cov = bandwidth_choice_2d(kernel, adhoc_raw_cov,
                                          cv_mode = bandwidth_cov_cv_mode,
                                          symmetric = TRUE)
    }
  }

  cat('Bandwidth for covariance surface are: ', bandwidth_cov[1],
      ' ', bandwidth_cov[2], '.\n', sep = "")

  # cov_surface_smoothed = lwls_2d(bandwidth_cov, kernel, raw_cov, w2_cov,
  #                                x_out1 = x_cov, x_out2 = x_cov)
  if (vanilla) {
    cov_surface_smoothed = lwls_2d_vanilla(bandwidth_cov, kernel, adhoc_raw_cov,
                                           x_out1 = x_cov, x_out2 = x_cov)
  }

  cov_surface_smoothed = (cov_surface_smoothed + t(cov_surface_smoothed)) / 2

  # Part III: Choose number of principal components functions.
  cat('Part III: Choose number of principal components functions.\n')

  if (is.numeric(out_percent)) {
    if (out_percent < 0) {
      out_percent = NULL
      cat('out_percent must between 0 to 1!\n Do not leave out boundary!\n')
    }
    if (out_percent > 0.25) {
      out_percent = 0.25
      cat('leaving out ', out_percent,
          ' percent of the data in the boundary is too much, reset out_percent to 0.25 !\n', sep = "")
    }
  }
  else {
    if (!is.null(out_percent)) {
      cat('out_percent input is unsupported!\n Do not leave out boundary!\n')
    }
    out_percent = NULL
  }

  # not clean here, combin userrange with no cut off as userrange = full in future
  if (!is.null(out_percent)) {
    cat('Use', out_percent, 'for leave out boundary.\n')
    userrange = as.numeric(quantile(x_data, c(out_percent, 1 - out_percent)))
  }
  if (!is.null(userrange)) {
    cat('Use userrange from', userrange[1],
        'to', userrange[2], 'for leave out boundary.\n')
    junk_index = x_data >= userrange[1] & x_data <= userrange[2]
    junk_x_data = x_data
    x_data = x_data[junk_index]
    junk_x_cov = x_cov
    x_cov = seq(min(x_data), max(x_data), length.out = cov_smooth_grid + 1)
    mu_data= mu_data[junk_index]
    mu_cov = approx(junk_x_cov, mu_cov, xout = x_cov)$y

    new_coordinate = expand.grid(x_cov, x_cov)
    cov_surface_smoothed = matrix(interp2(junk_x_cov, junk_x_cov, cov_surface_smoothed,
                                          new_coordinate[,2],  new_coordinate[,1],  method = "linear"),
                                  nrow = length(x_cov), byrow = FALSE)
    rownames(cov_surface_smoothed) = x_cov
    colnames(cov_surface_smoothed) = x_cov
    data = data[,junk_index]
    data = data[!apply(data, 1, function(x) all(is.na(x))),]
    rm(junk_index, junk_x_data, junk_x_cov)
  }

  if (error) {
    # should use quadratic form on diagonal to estimate Var(x(t))
    sigma = sigma_est(raw_cov, w2_cov, bandwidth_cov, kernel, userrange)
  }
  else {
    sigma = 0
  }

  h = x_cov[2] - x_cov[1]
  junk = eigen(cov_surface_smoothed, symmetric = TRUE)
  index = Im(junk$values) == 0 & Re(junk$values) > 0
  eigen_values = junk$values[index]
  eigen_func_cov = junk$vectors[,index]
  rm(junk)

  maxk = min(maxk, sum(index))
  if (no_of_components > maxk) {
    cat('At most ', maxk, ' of PC can be selected!\n', sep = "")
    no_of_components = maxk
  }

  lambda = h * eigen_values
  FVE = cumsum(lambda / sum(lambda))
  cat("FVE calculated from", length(lambda), "possible eigenvalues:\n ", FVE[1:min(5, maxk)], "...\n")
  if (missing(no_of_components) | no_of_components == -1) {
    if (selection == "FVE") {
      no_of_components = min(which(FVE >= FVE_threshold), maxk)
      if (no_of_components < 2 & maxk > 1) {
        no_of_components = 2
        cat("Use less than 2 components, reset to 2.\n")
      }
      cat("Best number of principal components selected by FVE: ", no_of_components, ".\n",
          "It accounts for ", FVE[no_of_components] * 100, "% of total variation (threshold = ", FVE_threshold, ").\n", sep = "")
    }
    if (selection == "BIC") {
      no_of_components = 2
      cat("No implement for BIC yet.")
    }
    if (selection == "AIC") {
      no_of_components = 2
      cat("No implement for AIC yet.")
    }
  }

  eigen_values = eigen_values[1:no_of_components]
  lambda = h * eigen_values
  eigen_func_cov = eigen_func_cov[,1:no_of_components]

  for (i in 1:ncol(eigen_func_cov)) {
    eigen_func_cov[,i] = eigen_func_cov[,i] / sqrt(trapz(x_cov, eigen_func_cov[,i] ^ 2))
    ## no idea why do so...
    if (eigen_func_cov[1,i] > eigen_func_cov[2,i]) {
      eigen_func_cov[,i] = -eigen_func_cov[,i]
    }
  }

  eigen_func_data = matrix(as.double(NA), ncol = ncol(eigen_func_cov), nrow = length(x_data))
  for (i in 1:ncol(eigen_func_data)) {
    eigen_func_data[,i] = spline(x = x_cov, y = eigen_func_cov[,i], xout = x_data)$y
    eigen_func_data[,i] = eigen_func_data[,i] / sqrt(trapz(x_data, eigen_func_data[,i] ^ 2))
  }
  cov_surface_fitted = matrix(as.double(0), ncol = length(x_data), nrow = length(x_data))
  for (i in 1:ncol(eigen_func_data)) {
    cov_surface_fitted = cov_surface_fitted +
      lambda[i] * eigen_func_data[,i] %*% t(eigen_func_data[,i])
  }
  cor_surface_fitted = cov2cor(cov_surface_fitted)
  rownames(cov_surface_fitted) = x_data
  colnames(cov_surface_fitted) = x_data
  rownames(cor_surface_fitted) = x_data
  colnames(cor_surface_fitted) = x_data

  # Part IV: Perform principal components analysis.
  cat('Part IV: Perform principal components analysis.\n')

  hash = apply(data, 1, function(x) do.call(paste, as.list(as.integer(!is.na(x)))))
  unique_hash = unique(hash)
  # original getOriCurves.m

  # set up hash table to identify missing situation
  # combin same missing situation subjects to improve speed
  # not meanful when subjects are all different

  # always update sigma, means no case rho = -1 in paper

  # if needed
  # if (rho != -1)
  #   ...
  # else
  # sigma_new_2 = sigma

  if (error) {
    # \hat{\sigma}^2_{new,1} in FPCscore.pdf Step 1
    sigma_new_1 = pcs_est(data, mu_data, hash, unique_hash, lambda, eigen_func_data,
                          cov_surface_fitted, sigma, sigma_new_switch = TRUE)$sigma_new

    # \hat{\sigma}^2_{new,2} in FPCscore.pdf Step 2
    sigma_new_2 = pcs_est(data, mu_data, hash, unique_hash, lambda, eigen_func_data,
                          cov_surface_fitted, sigma_new_1, sigma_new_switch = TRUE)$sigma_new
    # Step 3
    if (rho_cv) {
      rho = cv_rho(data, mu_data, lambda, eigen_func_data, cov_surface_fitted)
    }
    sigma_new_2 = max(sigma_new_2, rho)
  }
  else {
    sigma_new_2 = 0
  }
  # original getScores1.m
  junk = pcs_est(data, mu_data, hash, unique_hash, lambda, eigen_func_data,
                 cov_surface_fitted, sigma_new_2, xi_var_switch = TRUE)
  if (!error) {
    sigma = NA
    rho_opt = NA
    sigmanew = NA
  }
  if (sparse_flag) {
    regular = paste(regular, '(convert from sparse)')
  }
  ans = list(# mean related
    x_data = x_data, mu_data = mu_data,
    x_cov = x_cov, mu_cov = mu_cov,
    bandwidth_mean = bandwidth_mean,
    # covariance(correlation) related
    cov_surface_smoothed = cov_surface_smoothed,
    cov_surface_fitted = cov_surface_fitted,
    bandwidth_cov = bandwidth_cov,
    cor_surface_fitted = cor_surface_fitted,
    # PCA related
    no_of_components = no_of_components, lambda = lambda,
    # variables below are formatted by subjects
    # eigen function value at input points
    eigen_func_data = eigen_func_data,
    # eigen function value at output points
    eigen_func_cov = eigen_func_cov,
    # Principal component scores and so on
    xi_est = junk$xi_est,
    y_pred = junk$y_pred,
    xi_var = junk$xi_var,
    # etc
    sigma = sigma, regular = regular,
    AIC = NA, BIC = NA, FVE = FVE,
    # parameters for restricting the \hat{\sigma}^2_{new,2}
    rho_opt = rho, sigmanew = sigma_new_2)
  class(ans) = "FPCA"
  ans
}
