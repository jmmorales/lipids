library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
C = as.matrix(read.csv('https://github.com/jmmorales/lipids/raw/master/data/C.csv'))
Ys = read.csv('https://github.com/jmmorales/lipids/raw/master/data/Ys.csv')
TT = as.matrix(read.csv('https://github.com/jmmorales/lipids/raw/master/data/TT.csv'))
X = as.matrix(read.csv('https://github.com/jmmorales/lipids/raw/master/data/X.csv'))
Dmat = as.matrix(read.csv('https://github.com/jmmorales/lipids/raw/master/data/Dmat.csv'))
DistP = as.matrix(read.csv('https://github.com/jmmorales/lipids/raw/master/data/DistP.csv'))

dat <- list(N = 2309, 
            Y  = Ys$Y, 
            X = X, 
            K = 2, 
            J = 155, 
            jj = Ys$jj, 
            L=3, 
            TT = TT, 
            C=C, 
            M_1=1,
            M_2=1,
            Dmat=Dmat,
            J_1=Ys$J_1, 
            J_2=Ys$J_2,
            N_1=70,
            N_2=83,
            DistP=DistP, 
            ones=numeric(155)+1)
fit <- stan(file='stancode.stan', data = dat, iter = 10, chains = 1)
save(fit, file="fit.RData")

system(mutt -s “Master, this is an email from HAL 9000, your commands have been completed” pajarom@gmail.com < ~/R/output.txt)