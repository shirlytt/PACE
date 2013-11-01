options(error = recover)
source('lwls.R')
data = matrix(rnorm(51*51,0,10), ncol = 51, byrow = TRUE)
colnames(data) = seq(0, 10, length.out = 51)
cov_matrix = cov(data)
bandwidth = c(0.5, 0.5)
kernel = 0
x_in_1 = rownames(cov_matrix)
x_in_2 = colnames(cov_matrix)
y_in = cov_matrix
count_in = y_in / y_in

adhoc_cov = adhoc_regular_cov(y_in, count_in)

junk1 = lwls_2d(bandwidth, kernel, x_in_1, x_in_2, y_in, count_in)
junk2 = lwls_2d_vanilla(bandwidth, kernel, adhoc_cov)
