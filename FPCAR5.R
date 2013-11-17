require(caTools)
require(pracma)

# options(warn = 2)
options(error = recover)

# source local weighted least square functions
source('lwls.R')

PACE = setRefClass(Class = "PACE",
  fields = list(
      # data after binning
      data = "matrix",
      data_cutted = "matrix",    # cutted
      num_subject = "numeric",
      subject_hash = "character",
      regular = "character",
      bandwidth_mean = "numeric",
      bandwidth_cov = "numeric",
      num_bins = "numeric",
      cov_grid = "numeric",
      x_data = "numeric",
      mu_data = "numeric",
      x_cov = "numeric",
      mu_cov = "numeric",
      raw_cov = "matrix",
      y_var = "numeric",
      2d_weight = "matrix",
      kernel = "character",
      cov_surface_smoothed = "matrix", # from cutted
      cov_surface_fitted = "matrix",   # from cutted
      cor_surface_fitted = "matrix",   # from cutted
      x_data_cutted = "numeric",  # from cutted
      mu_data_cutted = "numeric", # from cutted
      x_cov_cutted = "numeric",   # from cutted
      mu_cov_cutted = "numeric",  # from cutted
      eigen_values = "numeric",   # from cutted
      eigen_func_data = "list",   # from cutted
      eigen_func_cov = "list",    # from cutted
      FVEs = "numeric",
      no_of_components = "numeric",
      bandwidth_mean_cv_mode = "character",
      bandwidth_cov_cv_mode = "character",
      FVE_threshold = "numeric",
      error = "logical",
      sigma = "numeric",
      new_ridge = "logical",
      rho_cv = "logical",
      rho = "numeric",
      out_percent = "numeric",
      datarange = "numeric",
      userrange = "numeric",
      xi_est = "list",
      y_pred = "matrix",
      xi_var = "list"             # indexed by hash
    ),
  methods = list(
    initialize = function(x , y
      bandwidth_mean = 0,
      bandwidth_cov = 0,
      num_bins = NULL,
      cov_grid = 50,
      kernel = c("Auto",
                 "Epanechnikov",
                 "Rectangular",
                 "Gaussian",
                 "Quartic",
                 "Gausvar"),
      bandwidth_mean_cv_mode = c("geo_mean", "gcv", "cv"),
      bandwidth_cov_cv_mode = c("geo_mean", "gcv", "cv"),
      FVE_threshold = 0.95,
      error = TRUE,
      out_percent = 0)
    {
      num_subject <<- length(x)
      kernel <<- match.arg(kernel)
      bandwidth_mean <<- bandwidth_mean,
      bandwidth_cov <<- bandwidth_cov,

      bandwidth_mean_cv_mode <<- match.arg(bandwidth_mean_cv_mode)
      bandwidth_mean_cv_mode = switch(bandwidth_mean_cv_mode,
                                      "cv" = 1,
                                      "gcv" = 2,
                                      "geo_mean" = 3)
      bandwidth_cov_cv_mode <<- match.arg(bandwidth_cov_cv_mode)
      bandwidth_cov_cv_mode = switch(bandwidth_cov_cv_mode,
                                     "cv" = 1,
                                     "gcv" = 2,
                                     "geo_mean" = 3)
      cov_grid <<- cov_grid
      error <<- error
      # check regular
      xx = unlist(x)
      datarange <<- range(xx)
      check_regular(length(xx) / length(unique(xx)) / num_subject)
      rm(xx)
      ## binning
      # avoid copy x and y, no separate function for binning
      # no of bins
      if (is.null(num_bins)) {
        m = median(sapply(x, length))
        mm = max(sapply(x, length))
        if (num_subject <= 5000) {
          if (mm <= 400 | m <= 20) {
            if (regular == 'sparse') {
              # irregular or sparse
              num_bins <<- 50
            }
            else {
              num_bins <<- mm
            }
          }
          else {
            num_bins <<- 400
          }
        }
        else {
          m_star = max(20, round((5000 - num_subject) * 19 / 2250 + 400))
          if (m > m_star) {
            num_bins <<- m_star
          }
          else {
            if (regular == 'sparse') {
              # irregular or sparse
              num_bins <<- 50
            }
            else {
              num_bins <<- mm
            }
          }
          rm(m_star)
        }
        rm(m, mm)
      }
      else {
        if (num_bins == 0) {
          num_bins <<- length(x[[1]])
        }
        else {
          num_bins <<- num_bins
        }

      }
      # binning grid
      d = diff(datarange) / num_bins
      x_data <<- seq(datarange[1] + d / 2,
                     datarange[2] - d / 2,
                     length.out = num_bins)
      d = diff(datarange) / cov_grid
      x_cov <<- seq(datarange[1] + d / 2,
                    datarange[2] - d / 2,
                    length.out = cov_grid)
      grids = seq(datarange[1], datarange[2], length.out = num_bins + 1)
      # binning
      for (i in 1:num_subject) {
        data[i,] <<- as.vector(by(y[[i]], cut(x[[i]], grids, include.lowest = TRUE), mean))
      }
      rm(d, grids)

      # Part I: Obtain smoothed mean curve.
      cat('Part I: Obtain smoothed mean curve.\n')
      set_mean()

      # Part II: Obtain smoothed covariance surface.
      cat('Part II: Obtain smoothed covariance surface.\n')
      # set raw cov mat
      set_raw_cov()
      # set smoothed cov mat
      set_smoothed_cov()

      # Part III: Choose number of principal components functions.
      cat('Part III: Choose number of principal components functions.\n')
      cut_off(out_percent)
      estimate_sigma()
      set_eigens()
      set_k()

      # Part IV: Perform principal components analysis.
      cat('Part IV: Perform principal components analysis.\n')
      set_ridge()
      pcs_est()
    },
    check_regular = function(f) {
      regular <<- 'sparse'
      if (f > 0.75) {
        regular <<- 'missing'
      }
      if (f == 1) {
        regular <<- 'regular'
      }
      cat('The design is ', regular, '.\n', sep = "")
    },
    intkernel = function() {
      switch(kernel,
             "Epanechnikov" = 0,
             "Rectangular" = 1,
             "Gaussian" = 2,
             "Quartic" = 3,
             "Gausvar" = 4)
    },
    set_mean = funcion() {
      if (kernel == "Auto") {
        if(regular == 'regular' & num_bins >= 20) {
          kernel = "Epanechnikov"
        }
        else {
          kernel = "Gaussian"
        }
      }
      kernel <<- kernel
      kernel = intkernel(kernel)
      # x_data = x_data
      y_data = as.double(colMeans(data, na.rm = TRUE))
      w_data = as.double(colSums(!is.na(data)))
      if (bandwidth_mean == 0) {
        bandwidth_mean = bandwidth_choice_1d(kernel = as.integer(kernel),
                                             x_in = x_data, y_in = y_data, w_in = w_data,
                                             cv_mode = bandwidth_mean_cv_mode)
      }
      bandwidth_mean <<- bandwidth_mean
      cat('Bandwidth for mean function is: ', bandwidth_mean, '.\n', sep = "")
      mu_data <<- lwls(bandwidth_mean, kernel, x_data, y_data, w_data, x_out = x_data)$output
      mu_cov <<- lwls(bandwidth_mean, kernel, x_data, y_data, w_data, x_out = x_cov)$output
    },
    set_raw_cov = function() {
      data = as.matrix(sweep(data, 2, mu_data, '-'))
      if (regular == 'regular') {
        raw_cov <<- t(data) %*% data / n
        2d_weight <<- matrix(n, ncol = ncol(data), nrow = ncol(data))
      }
      else {
        raw_cov <<- matrix(NA, p, p)
        2d_weight <<- matrix(NA, p, p)
        for (i in 1:p) {
          for (j in i:p) {
            raw_cov[i,j] <<- raw_cov[j,i] = mean(data[,i] * data[,j], na.rm = TRUE)
            2d_weight[i,j] <<- 2d_weight[j,i] = sum(!(is.na(data[,i]) | is.na(data[,j])))
          }
        }
      }
      y_var <<- diag(raw_cov)
      if (error) {
        diag(raw_cov) <<- NA
        diag(2d_weight) <<- 0
      }
    },
    set_smoothed_cov = function() {
      kernel = intkernel(kernel)
      if (bandwidth_cov == 0) {
        bandwidth_cov = bandwidth_choice_2d(kernel = as.integer(kernel),
                                            x_in = x_data, y_in = raw_cov, w_in = 2d_weight,
                                            cv_mode = bandwidth_cov_cv_mode)
      }
      bandwidth_cov <<- bandwidth_cov
      cat('Bandwidth for covariance surface are: ', bandwidth_cov[1],
          ' ', bandwidth_cov[2], '.\n', sep = "")
      cov_surface_smoothed <<- lwls_2d(bandwidth_cov, as.integer(kernel),
                                       x_in = x_data, y_in = raw_cov, w_in = 2d_weight,
                                       x_out1 = x_cov, x_out2 = x_cov)
      cov_surface_smoothed <<- (cov_surface_smoothed + t(cov_surface_smoothed)) / 2
    },
    cut_off = function(out_percent) {
      # set out_percent
      if (is.numeric(out_percent)) {
        if (out_percent < 0) {
          cat('out_percent must between 0 to 1!\n Do not leave out boundary!\n')
          out_percent = 0
        }
        if (out_percent > 0.25) {
          cat('leaving out ', out_percent, ' percent of the data in the boundary is too much, reset out_percent to 0.25 !\n', sep = "")
          out_percent = 0.25
        }
      }
      else {
        cat('out_percent input is unsupported!\n Do not leave out boundary!\n')
        out_percent = 0
      }
      out_percent <<- out_percent
      cat('Use', out_percent, 'for leave out boundary.\n')
      userrange <<- as.numeric(quantile(x_data, c(out_percent, 1 - out_percent)))
      # set
      if (out_percent != 0) {
        junk_index = x_data >= userrange[1] & x_data <= userrange[2]
        x_data_cutted <<- x_data[junk_index]
        mu_data_cutted <<- mu_data[junk_index]
        data_cutted <<- data[,junk_index]
        x_cov_cutted <<- seq(min(x_data_cutted), max(x_data_cutted), length.out = cov_grid)
        mu_cov_cutted= approx(x_cov, mu_cov, xout = x_cov_cutted)$y

        new_coordinate = expand.grid(x_cov_cutted, x_cov_cutted)
        cov_surface_smoothed <<- matrix(interp2(x_cov, x_cov, cov_surface_smoothed,
                                                new_coordinate[,2],  new_coordinate[,1],
                                                method = "linear"), nrow = length(x_cov), byrow = FALSE)
      }
      else {
        x_data_cutted <<- x_data
        mu_data_cutted <<- mu_data
        x_cov_cutted <<- x_cov
        mu_cov_cutted <<- mu_cov
      }
    },
    estimate_sigma = function() {
      if (error) {
        y_var_cut = approx(x_data, y_var, x_data_cutted)$y
        ...
      }
      else {
        sigma <<- 0
      }
    },
    set_eigens = function() {

    },
    set_ridge = function() {

    },
    pcs_est = function() {

    }
    ))
