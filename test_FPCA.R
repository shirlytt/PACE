source('simfunc.R')
source('FPCA.R')
require(caTools)
options(warn = 2)

Nsimu = 1
n = 100
# for (setting in 1:4) {
setting = 4
# for (shape in 1:2) {
shape = 2
# for (error in c(0,1,5)) {
error = 1
# for (CK in c(2, 5)) {
CK = 2
ext = paste0('n', n, 'setting', setting, 'shape', shape, 'error', error, 'CK', CK)

DataS = list()
for(i in 1:Nsimu){
  DataS[[i]] = SimuData(n,setting,shape,error,CK)
}
i = 1
y = DataS[[i]]$y
t = DataS[[i]]$x
x = FPCA(y = y, x = t, out_percent = 0.05, rho_cv = TRUE)


sigmas = list()
lambdas = list()
eigenfuncs = list()
mus = list()

for (i in 1:length(DataS)) {
  # data = matrix(unlist(DataS[[i]]$y), ncol = 51, byrow = TRUE)
  # colnames(data) = DataS[[i]]$x[[1]]
  y = DataS[[i]]$y
  t = DataS[[i]]$x
  if (error > 0) {
    junk = try(FPCA(y = y, t = t), silent = TRUE)
  }
  else {
    junk = try(FPCA(y = y, t = t, error = FALSE), silent = TRUE)
  }
  if (class(junk) != "try-error") {
    sigmas[[i]] = junk$sigma
    lambdas[[i]] = junk$lambda
    eigenfuncs[[i]] = junk$eigen_func_ins
    mus[[i]] = junk$mu_in
  }
}

if (error > 0) {
  sigma_bias = mean(unlist(sigmas)) - 1
  sigma_var = var(unlist(sigmas))
  sigma_mse = sigma_var + sigma_bias ^ 2
}

lambdas = lambdas[!sapply(lambdas, is.null)]
first_lambda_bias = mean(sapply(lambdas, function(x) x[[1]])) - 9
first_lambda_var = var(sapply(lambdas, function(x) x[[1]]))
first_lambda_mse = first_lambda_var + first_lambda_bias ^ 2

second_lambda_bias = mean(sapply(lambdas, function(x) x[[2]])) - 6
second_lambda_var = var(sapply(lambdas, function(x) x[[2]]))
second_lambda_mse = first_lambda_var + first_lambda_bias ^ 2

# no_of_comp = table(sapply(lambdas, length))

time_point = seq(0,10,length.out=21)
true_eigenfunc = xeig(time_point,10,2)
true_mu = time_point+sin(time_point)
true_mu = true_mu + dnorm(time_point,mean = 5,sd = 2/3)


mus = mus[!sapply(mus, is.null)]
eigenfuncs = eigenfuncs[!sapply(eigenfuncs, is.null)]
mu_imse = list()
eigenfunc1_imse = list()
eigenfunc2_imse = list()
for (i in 1:length(mus)) {
  mu_imse[[i]] = trapz(time_point, (mus[[i]] - true_mu)^2)
  eigenfunc1_imse[[i]] = trapz(time_point, (eigenfuncs[[i]][,1] - true_eigenfunc[1,])^2)
  eigenfunc2_imse[[i]] = trapz(time_point, (eigenfuncs[[i]][,2] - true_eigenfunc[1,])^2)
}

mu_imse = unlist(mu_imse)
eigenfunc1_imse = unlist(eigenfunc1_imse)
eigenfunc2_imse = unlist(eigenfunc2_imse)

if (error > 0) {
  save(sigma_bias,sigma_var,sigma_mse,first_lambda_bias,first_lambda_var,first_lambda_mse,second_lambda_bias,second_lambda_var,second_lambda_mse,no_of_comp,mu_imse,eigenfunc1_imse,eigenfunc2_imse, file = paste0(ext, '.rda'))
  png(paste0(ext, 'sigma.png'))
  plot(density(unlist(sigmas)))
  abline(v=1)
  dev.off()
}
else {
  save(first_lambda_bias,first_lambda_var,first_lambda_mse,second_lambda_bias,second_lambda_var,second_lambda_mse,no_of_comp,mu_imse,eigenfunc1_imse,eigenfunc2_imse, file = paste0(ext, '.rda'))
}

png(paste0(ext, 'first_lambda.png'))
plot(density(sapply(lambdas, function(x) x[[1]])))
abline(v=9)
dev.off()

png(paste0(ext, 'second_lambda.png'))
plot(density(sapply(lambdas, function(x) x[[2]])))
abline(v=6)
dev.off()

png(paste0(ext, 'noc.png'))
hist(sapply(lambdas, length))
dev.off()

# }
# }
# }
# }
