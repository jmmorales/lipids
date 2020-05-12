library(rstan) # observe startup messages
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

C = as.matrix(read.csv('https://github.com/jmmorales/lipids/raw/master/data/C.csv'))
Ys = read.csv('https://github.com/jmmorales/lipids/raw/master/data/Ys.csv')
TT = as.matrix(read.csv('https://github.com/jmmorales/lipids/raw/master/data/TT.csv'))
X = as.matrix(read.csv('https://github.com/jmmorales/lipids/raw/master/data/X.csv'))
Dmat = as.matrix(read.csv('https://github.com/jmmorales/lipids/raw/master/data/Dmat.csv'))
DistP = as.matrix(read.csv('https://github.com/jmmorales/lipids/raw/master/data/DistP.csv'))