# may not useful in future
lwls_1d = function(bandwidth, kernel, data, 
 weight_out_need = FALSE, x_out = data$newx, 
 power = 1, cv_mode = 0) {
  if (class(data) != "FDA_1d_binned_data") {
    stop("The class of data must be FDA_1d_binned_data")
  }
  lwls(bandwidth, kernel, data$newx, data$newy, data$count / sum(data$count), 
       weight_out_need, 
       x_out, power, cv_mode)
}

lwls2d = function(bandwidth, kernel, data, x_out1, x_out2, power = 1, cv_mode = 0) {
  # bandwidth: a vactor of 2 numeric
  if (class(data) != "FDA_2d_binned_data") {
    stop("The class of data must be FDA_2d_binned_data")
  }
  ans = list()
  if (missing(x_out1)) {
    x_out1 = as.numeric(colnames(data$newy))
  }
  if (missing(x_out2)) {
    x_out2 = as.numeric(rownames(data$newy))
  }
  x_in1 = as.numeric(colnames(data$newy))
  x_in2 = as.numeric(rownames(data$newy))
  non_spare_row = apply(data$count, 1, function(x) sum(x != 0)) > 3
  y = data$newy[non_spare_row,]
  count = data$count[non_spare_row,]
  # for loop for first try, may use mapply instead in future
  # native c?

  middle_stage1 = matrix(as.numeric(NA), nrow = nrow(y), ncol = length(x_out1))
  mu1 = matrix(as.numeric(NA), nrow = length(x_out2), ncol = length(x_out1))
  # # get mu1
  # for (i in 1:nrow(y)) {
  #   ind = count[i,] > 0
  #   middle_stage1[i,] = lwls(bandwidth[1], kernel, x_in1[ind], y[i,ind], count[i,ind] / sum(count[i,ind]), 
  #     x_out = x_out1, power = 1, cv_mode = cv_mode)$output
  # }
  # for (i in 1:length(x_out1)) {
  #   mu1[,i] = lwls(bandwidth[2], kernel, x_in2, middle_stage1[,i], rep(1 / nrow(mu1), nrow(mu1)), 
  #     x_out = x_out2, power = 1, cv_mode = 0)$output
  # }
  # colnames(mu1) = x_out1
  # rownames(mu1) = x_out2
  # ans$mu1 = mu1


  min_bandwidth = find_min_bandwidth(x_in1)
  max_bandwidth = find_max_bandwidth(x_in1)
  if (cv_mode > 0) {
    ad_hoc_lwls = function(bw) {
      ans = 0
      for (i in 1:nrow(y)) {
        ind = count[i,] > 0
        ans = ans + lwls(bw, kernel, x_in1[ind], y[i,ind], count[i,ind] / sum(count[i,ind]), 
          x_out = x_out1, power = 1, cv_mode = cv_mode)$cv_value
      }
      ans
    }
    # bandwidth[1] = asNumeric(optimizeR(ad_hoc_lwls, lower = min_bandwidth, upper = max_bandwidth)$minimum)
    bandwidth[1] = optimize(ad_hoc_lwls, lower = min_bandwidth, upper = max_bandwidth)$minimum
    if (cv_mode == 3) {
      bandwidth[1] = sqrt(bandwidth[1] * min_bandwidth)
    }
    cat("1st bandwidth: ", bandwidth[1], "\n")
  }

  for (i in 1:nrow(y)) {
    ind = count[i,] > 0
    # middle_stage1[i,] = lwls(bandwidth[1], kernel, x_in1[ind], y[i,ind], count[i,ind] / sum(count[i,ind]), 
    #   x_out = x_out1, power = 1, cv_mode = cv_mode)$output
    # w = sqrt(count[i,ind] / sum(count[i,ind]))
    w = sqrt(count[i,ind])
    middle_stage1[i,] = lwls(bandwidth[1], kernel, x_in1[ind], y[i,ind], w / sum(w), 
      x_out = x_out1, power = 1, cv_mode = cv_mode)$output
  }

  min_bandwidth = find_min_bandwidth(x_in2)
  max_bandwidth = find_max_bandwidth(x_in2)
  if (cv_mode > 0) {
    ad_hoc_lwls = function(bw) {
      ans = 0
      for (i in 1:length(x_out1)) {
        ans = ans + lwls(bw, kernel, x_in2, middle_stage1[,i], sqrt(count[i,ind]) / sum(sqrt(count[i,ind])), 
          x_out = x_out2, power = 1, cv_mode = 0)$cv_value
      }
      ans
    }
    # bandwidth[2] = asNumeric(optimizeR(ad_hoc_lwls, lower = min_bandwidth, upper = max_bandwidth)$minimum)
    bandwidth[2] = optimize(ad_hoc_lwls, lower = min_bandwidth, upper = max_bandwidth)$minimum
    if (cv_mode == 3) {
      bandwidth[2] = sqrt(bandwidth[2] * min_bandwidth)
    }
    cat("2nd bandwidth: ", bandwidth[2], "\n")
  }

  for (i in 1:length(x_out1)) {
    # mu1[,i] = lwls(bandwidth[2], kernel, x_in2, middle_stage1[,i], rep(1 / nrow(mu1), nrow(mu1)), 
    #   x_out = x_out2, power = 1, cv_mode = 0)$output
    # w = sqrt(count[,i] / sum(count[,i]))
    w = sqrt(count[,i])
    mu1[,i] = lwls(bandwidth[2], kernel, x_in2, middle_stage1[,i], w / sum(w), 
      x_out = x_out2, power = 1, cv_mode = 0)$output
  }
  colnames(mu1) = x_out1
  rownames(mu1) = x_out2
  
  ans$mu1 = mu1

  # smooth in another order
  #skip for now
  middle_stage2 = matrix(as.numeric(NA), nrow = length(x_out2), ncol = ncol(y))
  mu2 = matrix(as.numeric(NA), nrow = length(x_out2), ncol = length(x_out1))
  # get mu2
  for (i in 1:ncol(y)) {
    ind = count[,i] > 0
    # middle_stage2[,i] = lwls(bandwidth[2], kernel, x_in2[ind], y[ind,i], count[ind,i] / sum(count[ind,i]), 
    #   x_out = x_out2, power = 1, cv_mode = 0)$output
    # w = sqrt(count[ind,i] / sum(count[ind,i]))
    w = sqrt(count[ind,i])
    middle_stage2[,i] = lwls(bandwidth[2], kernel, x_in2[ind], y[ind,i], w / sum(w), 
      x_out = x_out2, power = 1, cv_mode = 0)$output
  }
  for (i in 1:length(x_out2)) {
    # mu2[i,] = lwls(bandwidth[1], kernel, x_in1, middle_stage2[i,], rep(1 / ncol(mu2), nrow(mu2)), 
    #   x_out = x_out1, power = 1, cv_mode = 0)$output
    # w = sqrt(count[i,] / sum(count[i,]))
    w = sqrt(count[i,])
    mu2[i,] = lwls(bandwidth[1], kernel, x_in1, middle_stage2[i,], w / sum(w), 
      x_out = x_out1, power = 1, cv_mode = 0)$output
  }
  colnames(mu2) = x_out1
  rownames(mu2) = x_out2

  ans$mu2 = mu2

  ans
}
