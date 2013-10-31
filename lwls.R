dyn.load('~/Dropbox/paCe/lwls.so')

binning = function(x, y, num_bins = 51, a0 = min(x), b0 = max(x)) {
  valid = (!is.na(x) & !is.na(y))
  x = x[valid]
  y = y[valid]
  bins = seq(a0, b0, length.out = num_bins + 1)
  ind = cut(x, bins, include.lowest = TRUE)
  count = table(ind)
  newx = ((bins[-length(bins)] + bins[-1]) / 2)[count > 0]
  newy = as.vector(tapply(y, ind, mean)[count > 0])
  count = as.vector(count)
  ans = list(newx = newx, newy = newy, count = count[count > 0])
  class(ans) = "FDA_1d_binned_data"
  ans
}

find_min_bandwidth = function(x_in, number_of_points = 3) {
  x_in = sort(unique(unlist(x_in)))
  d = x_in[number_of_points:length(x_in)] - x_in[1:(length(x_in) - number_of_points + 1)]
  min(d)
}

find_max_bandwidth = function(x_in) {
  max(x_in) - min(x_in)
}

lwls = function(bandwidth, kernel, x_in, y_in, 
 count_in = rep(1 / length(x_in), length(x_in)), 
 weight_out_need = FALSE,
 x_out, power = 1, cv_mode = 0) {
  ##kernel:
  #0 epanechnikov
  #1 rectangle
  #2 gaussian
  #3 quartic
  #4 variant of Gaussian

  ##power: 
  #1 mean
  #2 1st derivative
  #3 2nd derivative

  #cv_mode
  #0 no cv calculation
  #1 ordinary cross-validation
  #2 generalized cross-validation
  #3 geometric mean of GCV and minimum bandwidth
  
  if (missing(x_out)) {
    x_out = x_in
  }

  # user should not rely on this check
  valid = (!is.na(x_in) & !is.na(y_in) & !is.na(count_in) & !is.infinite(x_in) & !is.infinite(y_in) & !is.infinite(count_in))
  x_in = x_in[valid]
  y_in = y_in[valid]
  count_in = count_in[valid]

  # comment this part before we figure out how to calculate y error

  # out = x_out > max(x_in)
  # if (sum(out) > 0) {
  #   warning('Some x_out which is larger than valid maximum x_in will be removed!')
  # }
  # x_out = x_out[!out]
  # out = x_out < min(x_in)
  # if (sum(out) > 0) {
  #   warning('Some x_out which is smaller than valid minimum x_in will be removed!')
  # }
  # x_out = x_out[!out]

  if (any(diff(x_in) > ((max(x_in) - min(x_in)) / 4))) {
    stop('The data is too sparse!')
  }

  bandwidth = as.double(bandwidth)
  kernel = as.integer(kernel)
  x_in = as.double(x_in)
  y_in = as.double(y_in)
  count_in = as.double(count_in)
  n_in = as.integer(length(x_in))
  n_out = as.integer(length(x_out))
  # here user do not expect weight_out since it's for 2nd
  # smoothing in 2d situation
  weight_out_need = as.integer(weight_out_need)
  x_out = as.double(x_out)
  output = as.double(rep(0, length(x_out)))
  weight_out = as.double(rep(0, length(x_out)))
  power = as.integer(power)
  cv_mode = as.integer(cv_mode)
  cv_value = as.double(0)
  .C("lwls", 
    bandwidth = bandwidth, kernel = kernel, x_in = x_in, y_in = y_in, count_in = count_in, 
    n_in = n_in, n_out = n_out, weight_out_need = weight_out_need,
    x_out = x_out, output = output, weight_out = weight_out,
    power = power, cv_mode = cv_mode, cv_value = cv_value)
}

bandwidth_choice_1d = function(kernel = 0, x_in, y_in, w_in, cv_mode = 3) {
  ad_hoc_lwls = function(bw) {
    lwls(bw, kernel, x_in, y_in, w_in, 
     weight_out_need = FALSE, power = 1, cv_mode = cv_mode)$cv_value
  }
  optimize(ad_hoc_lwls, lower = find_min_bandwidth(x_in), upper = find_max_bandwidth(x_in))$minimum
}

# binning2d = function(x1, x2, y, num_bins1 = 51, num_bins2 = 51, 
#  a10 = min(x1), b10 = max(x1), a20 = min(x2), b20 = max(x2)) {
#   bins1 = seq(a10, b10, length.out = num_bins1 + 1)
#   bins2 = seq(a20, b20, length.out = num_bins2 + 1)
#   newx1 = (bins1[-length(bins1)] + bins1[-1]) / 2
#   newx2 = (bins2[-length(bins2)] + bins2[-1]) / 2
#   output = as.numeric(rep(0, num_bins1 * num_bins2))
#   count = as.integer(rep(0, num_bins1 * num_bins2))
#   ans = .C("binning2d", as.numeric(x1), as.numeric(x2), as.numeric(y), as.integer(length(y)),
#     as.numeric(a10), as.numeric(a20), 
#     as.numeric(bins1[2] - bins1[1]), as.numeric(bins2[2] - bins2[1]), 
#     as.integer(num_bins1), as.integer(num_bins2), output = output, count = count)
#   newy = matrix(ans$output, nrow = num_bins2, ncol = num_bins1, byrow = TRUE)
#   count = matrix(ans$count, nrow = num_bins2, ncol = num_bins1, byrow = TRUE)
#   colnames(newy) = newx1
#   rownames(newy) = newx2
#   colnames(count) = newx1
#   rownames(count) = newx2
# 
#   ans = list(newy = newy, count = count)
#   class(ans) = "FDA_2d_binned_data"
#   ans
# }


# lwls_seq is just a part of lwls_2d
# user should not call this function directly
# WARNING:
# lwls_seq or/and lwls_2d is wrong, DO NOT use it until fixed!
# it's worth to use lwls-seq since it's much faster than vanilla
#   user  system elapsed
#  0.107   0.000   0.107
#   user  system elapsed
#  1.637   0.000   1.639
# however vanilla matches matlab result
lwls_seq = function(bandwidth, kernel, x_in, y_in, w_in, 
 weight_out_need = TRUE,
 x_out, power = 1, cv_mode = 0) {
  if (as.integer(length(y_in) / length(x_in)) * length(x_in) != length(y_in)) {
    stop('Wrong data size!')
  }
  if (missing(x_out)) {
    x_out = x_in
  }
  bandwidth = as.double(bandwidth)
  kernel = as.integer(kernel)
  x_in = as.double(x_in)
  y_in = as.double(y_in)
  n_in = as.integer(length(y_in))
  x_out = as.double(x_out)
  n_out = as.integer(length(y_in) / length(x_in) * length(x_out))
  n_col = as.integer(length(y_in) / length(x_in))
  w_in = as.double(w_in / sum(w_in) * n_col)
  # w_in = as.double(w_in)
  weight_out_need = as.integer(weight_out_need)
  output = as.double(rep(0, n_out))
  weight_out = as.double(rep(0, n_out))
  power = as.integer(power)
  cv_mode = as.integer(cv_mode)
  cv_value = as.double(0)
  .C("lwls_seq", 
    bandwidth = bandwidth, kernel = kernel, x_in = x_in, y_in = y_in, w_in = w_in, 
    n_in = n_in, n_out = n_out, n_col = n_col, weight_out_need = weight_out_need,
    x_out = x_out, output = output, weight_out = weight_out,
    power = power, cv_mode = cv_mode, cv_value = cv_value, NAOK = TRUE)
}

lwls_2d = function(bandwidth, kernel, x_in_1, x_in_2, y_in, count_in, 
 x_out_1 = x_in_1, x_out_2 = x_in_2, power = 1, cv_mode = 0) {
  # in 1st smoothing we always need weight for 2nd step
  # x_in_1 is the index along rows
  middle_stage = lwls_seq(bandwidth[1], kernel, x_in_1, y_in, count_in, 
   weight_out_need = TRUE, x_out_1, power, cv_mode)
  # newy = matrix(middle_stage$output, ncol = length(y_in) / length(x_in_1), byrow = FALSE)
  # neww = matrix(middle_stage$weight_out, ncol = length(y_in) / length(x_in_1), byrow = FALSE)
  newy = matrix(middle_stage$output, ncol = length(y_in) / length(x_in_1), byrow = TRUE)
  neww = matrix(middle_stage$weight_out, ncol = length(y_in) / length(x_in_1), byrow = TRUE)
  # for test
  # neww = sweep(neww, 1, rowSums(neww), "/")

  ans = matrix(lwls_seq(bandwidth[2], kernel, x_in_2, newy, neww, 
                        weight_out_need = FALSE, x_out_2, power, cv_mode)$output, 
               nrow = length(x_out_1), byrow = TRUE)
  rownames(ans) = x_out_1
  colnames(ans) = x_out_2
  ans
}

# bandwidth_choice_2d_symmetrical = function()
# bandwidth_choice_2d_unsymmetrical = function()


# transfer data.frame to format lwls_2d_vanilla use
# will remove in future with lwls_2d_vanilla
# data is cov matrix
# junk is data
adhoc_regular_cov = function(cov, weight) {
  ans = list(newy = cov, count = weight)
  rownames(ans$count) = rownames(ans$newy)
  colnames(ans$count) = colnames(ans$newy)
  class(ans) = "FDA_2d_binned_data"
  ans
}

lwls_2d_vanilla = function(bandwidth, kernel, data, x_out1, x_out2, power = 1, cv_mode = c(0, 0)) {
  # bandwidth: a vactor of 2 numeric
  if (class(data) != "FDA_2d_binned_data") {
    stop("The class of data must be FDA_2d_binned_data")
  }
  ans = list()
  cv_value = 0
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

  middle_stage1 = matrix(as.numeric(NA), nrow = nrow(y), ncol = length(x_out1))
  middle_weight1 = matrix(as.numeric(NA), nrow = nrow(y), ncol = length(x_out1))
  mu1 = matrix(as.numeric(NA), nrow = length(x_out2), ncol = length(x_out1))
  for (i in 1:nrow(y)) {
    ind = count[i,] > 0
    w = count[i,ind]
    junk = lwls(bandwidth[1], kernel, x_in1[ind], y[i,ind], w / sum(w), weight_out_need = TRUE,
      x_out = x_out1, power = 1, cv_mode = cv_mode[1])
    # junk = lwls(bandwidth[1], kernel, x_in1[ind], y[i,ind], w, weight_out_need = TRUE,
    #   x_out = x_out1, power = 1, cv_mode = cv_mode)
    middle_stage1[i,] = junk$output
    middle_weight1[i,] = junk$weight_out
    # this is acceptable since it's cov matrix so there isn't too many rows
    if (cv_mode[1] > 0) {
      cv_value = junk$cv_value
    }
  }
  # should change return form in formal version
  if (cv_mode[1] > 0) {
    return(cv_value)
  }
  # middle_weight1 = sweep(middle_weight1, 1, rowSums(middle_weight1), "/")
  for (i in 1:length(x_out1)) {
    w = middle_weight1[,i]
    # mu1[,i] = lwls(bandwidth[2], kernel, x_in2, middle_stage1[,i], w / sum(w), 
    #   x_out = x_out2, power = 1, cv_mode = 0)$output
    mu1[,i] = lwls(bandwidth[2], kernel, x_in2, middle_stage1[,i], w, 
      x_out = x_out2, power = 1, cv_mode = cv_mode[2])$output
    if (cv_mode[2] > 0) {
      cv_value = junk$cv_value
    }
  }
  if (cv_mode[2] > 0) {
    return(cv_value)
  }
  colnames(mu1) = x_out1
  rownames(mu1) = x_out2

  mu1
  # ans$middle_stage1 = middle_stage1
  # ans$middle_weight1 = middle_weight1
  # ans$mu1 = mu1

  # # smooth in another order
  # #skip for now
  # middle_stage2 = matrix(as.numeric(NA), nrow = length(x_out2), ncol = ncol(y))
  # middle_weight2 = matrix(as.numeric(NA), nrow = length(x_out2), ncol = ncol(y))
  # mu2 = matrix(as.numeric(NA), nrow = length(x_out2), ncol = length(x_out1))
  # # get mu2
  # for (i in 1:ncol(y)) {
  #   ind = count[,i] > 0
  #   w = count[ind,i]
  #   junk = lwls(bandwidth[2], kernel, x_in2[ind], y[ind,i], w / sum(w), weight_out_need = TRUE,
  #     x_out = x_out2, power = 1, cv_mode = 0)
  #   # junk = lwls(bandwidth[2], kernel, x_in2[ind], y[ind,i], w, weight_out_need = TRUE,
  #   #   x_out = x_out2, power = 1, cv_mode = 0)
  #   middle_stage2[,i] = junk$output
  #   middle_weight2[,i] = junk$weight_out
  # }
  # # middle_weight2 = sweep(middle_weight2, 1, rowSums(middle_weight2), "/")
  # for (i in 1:length(x_out2)) {
  #   w = middle_weight2[i,]
  #   # mu2[i,] = lwls(bandwidth[1], kernel, x_in1, middle_stage2[i,], w / sum(w), 
  #   #   x_out = x_out1, power = 1, cv_mode = 0)$output
  #   mu2[i,] = lwls(bandwidth[1], kernel, x_in1, middle_stage2[i,], w, 
  #     x_out = x_out1, power = 1, cv_mode = 0)$output
  # }
  # colnames(mu2) = x_out1
  # rownames(mu2) = x_out2
  # ans$middle_stage2 = middle_stage2
  # ans$middle_weight2 = middle_weight2
  # ans$mu2 = mu2
  
  # ans
}

# also renew after replace vanilla
bandwidth_choice_2d = function(kernel, data, cv_mode, symmetric = TRUE) {
  # non-symmetric version
  ad_hoc_lwls = function(bw) {
    lwls_2d_vanilla(c(bw, NA), kernel, data, power = 1, cv_mode = c(cv_mode, 0))
  }
  x_in = as.numeric(rownames(data$newy))
  bw1 = optimize(ad_hoc_lwls, lower = find_min_bandwidth(x_in), upper = find_max_bandwidth(x_in))$minimum
  # ad_hoc_lwls = function(bw) {
  #   lwls_2d_vanilla(c(bw1, bw), kernel, data, power = 1, cv_mode = c(0, cv_mode))
  # }
  # bw2 = optimize(ad_hoc_lwls, lower = find_min_bandwidth(x_in), upper = find_max_bandwidth(x_in))$minimum
  c(bw1, bw1)
}
