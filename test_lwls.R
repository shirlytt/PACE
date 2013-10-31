x1 = as.numeric(readLines('x1'))
x2 = as.numeric(readLines('x2'))
y = as.numeric(readLines('y'))
# y = x2 * 3
source('lwls.R')
data = binning2d(x1, x2, y)
# ans = lwls2d(c(0.5,0.5), 2, data)
# anse = lwls2d(c(0.5,0.5), 0, data)
# ans = lwls2d(c(0.5,0.5), 2, data, cv_mode = 2)
# ans = lwls2d(c(0.5,0.5), 2, data)
bandwidth = 0.5
kernel = 0
x_in = as.numeric(rownames(data$newy))
y_in = as.numeric(data$newy)
w_in = as.numeric(data$count)
x_out = x_in
output = rep(0, length(y_in))
n_in = as.integer(length(y_in))
n_out = as.integer(length(output))
in_col_size = as.integer(nrow(data$newy))
out_col_size = as.integer(in_col_size)
power = as.integer(1)
cv_mode = as.integer(0)
cv_value = as.integer(0)

y_in_s = as.numeric(data$newy[1:length(x_in)])
w_in_s = as.numeric(data$count[1:length(x_in)])
output_s = rep(0, length(y_in_s))
n_in_s = as.integer(length(y_in_s))
n_out_s = as.integer(length(output_s))

single = .C('lwls', bandwidth = bandwidth, kernel = kernel, x_in = x_in, y_in = y_in_s, w_in = w_in_s, 
 x_out = x_out, output = output_s, n_in = n_in_s, n_out = n_out_s, 
 power = power, cv_mode = cv_mode, cv_value = cv_value)


ans = .C('lwls_seq', bandwidth = bandwidth, kernel = kernel, x_in = x_in, y_in = y_in, w_in = w_in, 
 x_out = x_out, output = output, n_in = n_in, n_out = n_out, in_col_size = in_col_size, out_col_size = out_col_size, 
 power = power, cv_mode = cv_mode, cv_value = cv_value)
