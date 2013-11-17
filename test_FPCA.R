source('simfunc.R')
source('FPCA.R')
require(caTools)
# options(warn = 2)

Nsimu = 100
n = 100
## for (setting in 1:4) {
setting = 4
## for (shape in 1:2) {
shape = 2
# for (error in c(0,1,5)) {
error = 1
# for (CK in c(2, 5)) {
CK = 2
ext = paste0('n', n, 'setting', setting, 'shape', shape, 'error', error, 'CK', CK)
if (setting == 4) {
  n = 500
}

DataS = list()
for(i in 1:Nsimu){
  DataS[[i]] = SimuData(n,setting,shape,error,CK)
}

sigmas = list()
lambdas = list()
eigenfuncs = list()
mus = list()
score_real = list()
score_est = list()
ir_x = list()

for (i in 1:length(DataS)) {
  y = DataS[[i]]$y
  t = DataS[[i]]$x
  x = try(FPCA(y = y, x = t, rho_cv = TRUE), silent = TRUE)
  if (class(x) != "try-error") {
    sigmas[[i]] = x$sigma
    lambdas[[i]] = x$lambda
    eigenfuncs[[i]] = x$eigen_func_data
    mus[[i]] = x$mu_data
    score_real[[i]] = DataS[[i]]$scores
    score_est[[i]] = x$xi_est
    ir_x[[i]] = x$x_data
  }
}

save(sigmas, lambdas, eigenfuncs, mus, score_real, score_est, ir_x, file = paste0(ext, '.rda'))
#}
#}
#}
#}

