library(ape)
library(geosphere)
library(viridis)

# prepare data and fit the model ----
filo <- read.nexus("data/output_bird_phylogeny_31Mar2020.nex")

# get a correlation matrix from the tree 
CC = vcv(filo, corr = TRUE)
#Cp = vcv(filop, corr = TRUE)

# sort species and re-arrange phylogenetic correlation
tmp <- dimnames(CC)
ids <- as.numeric(as.factor(tmp[[1]]))

C <- matrix(NA, ncol(CC), ncol(CC))
for(i in 1:ncol(CC)){
  for(j in 1:ncol(CC)){
    C[ids[i],ids[j]] <- CC[i,j]
  }
}

# cophenetic distance among plants
coph <- read.table("data/matriz_phylo_dist.txt", header = TRUE)

coP <- as.matrix(coph)

tmp <- dimnames(coP)
idsp <- as.numeric(as.factor(tmp[[1]]))

DistP <- matrix(NA, ncol(coP), ncol(coP))
for(i in 1:ncol(coP)){
  for(j in 1:ncol(coP)){
    DistP[idsp[i],idsp[j]] <- coP[i,j]
  }
}

datos <- read.csv("data/dados_spp_level.csv", header = TRUE)
# head(datos)

xy = unique(data.frame(site = datos$site, long=datos$lon,lat=datos$lat))
rownames(xy) = xy$site

Dmat <- distm(cbind(xy$long, xy$lat), fun = distHaversine)/1000000

visits <- (datos$visits)
lipid <- scale(datos$lipid)
site <- datos$site
spp <- datos$phylo # bird genus
sp <- sort(unique(spp))

pspp <- datos$Phylo_name # plant genus
psp <- unique(pspp)

np  <- 2    # number of parameters
nt  <- 2    # number of traits
npp <- 1
ntp <- 1

ns   <- length(sp)  # number of bird species
ns_p <- length(psp) # number of plant species
nsites <- length(unique(site))
nobs <- length(site)

# traits for birds
TT <- matrix(1, nrow = ns, ncol = nt + 1)
for(i in 1: ns){
  TT[i, 2] <- datos$Diet.Fruit[spp==sp[i]][1]
  TT[i, 3] <- log(datos$BodyMass[spp==sp[i]][1])
}
fd <- TT[,2]
bs <- TT[,3]

TT[, 2] <- scale(TT[, 2])
TT[, 3] <- scale(TT[, 3])


N <- length(visits)
Y <- visits
K <- np
X <- matrix(c(rep(1, N), lipid), N, K )

jj <- as.numeric(as.factor(datos$phylo))
J <- max(jj)

J_1 <- as.numeric(as.factor(datos$site))
J_2 <- as.numeric(as.factor(datos$Phylo_name))
N_1 <- max(J_1)
N_2 <- max(J_2)
L <- ncol(TT)
ones <- numeric(J) + 1

library(rstan)
#Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

pr_dat <- list(N = N,
               Y = Y,
               K = K,
               X = X,
               jj = jj,
               J = J,
               L = L,
               TT = TT,
               C = C,
               ones = ones,
               M_1 = 1,
               M_2 = 1,
               Dmat = Dmat,
               J_1 = J_1,
               J_2 = J_2,
               N_1 = N_1,
               N_2 = N_2,
               DistP = DistP/100)

fit <- stan(file = 'pg_nb_pp.stan', data = pr_dat,
            iter = 100000, thin = 50, chains = 4)
#, control = list(max_treedepth = 15))

# leave one out CV ----
library(loo)
# Extract pointwise log-likelihood and compute LOO
log_lik <- extract_log_lik(fit, merge_chains = FALSE)

r_eff <- relative_eff(exp(log_lik)) 

loo_1 <- loo(log_lik, r_eff = r_eff, cores = 2)
print(loo_1)
#--------------------------

# figure options
res = 300
he = 17
wi = 17

fit_summary <- summary(fit)$summary

# check R_hat
hist(fit_summary[,10], 100)

# check n_eff
hist(fit_summary[,9], 100)


# plot of regression coefficients against degree of frugivory and body size ----

mean_b <- fit_summary[grepl("b_m", rownames(fit_summary)),]
mean_z <- fit_summary[grepl("Z", rownames(fit_summary)),]
zs <- extract(fit, pars = "z")

nsim <- 1500
xp <- seq(min(TT[,2]), max(TT[,2]), length.out = 101)
xg <- seq(min(TT[,3]), max(TT[,3]), length.out = 100)
z1 <- matrix(NA, nsim, length(xp))
z2 <- matrix(NA, nsim, length(xp))
z3 <- matrix(NA, nsim, length(xg))
z4 <- matrix(NA, nsim, length(xg))


tiff(filename="sp-zetas.tif", res=res, unit="cm", height=17, width=17)
nf = layout(matrix(c(1,2,3,4,0,0),3,2,byrow=TRUE), widths=c(1,1), heights=c(1,1,0.1))
#layout.show(nf)
op <- par( mar = c(3, 3, 2, 2) + 0.1, las = 1, bty = "n", cex = 1)
plot(fd , mean_b[seq(1, (ns*2),by=2),1], pch = 19, col=rgb(0,0,0, 0.4), 
     ylab = "", xlab = "", main = "")
lines(fd, mean_z[1,1] + mean_z[3,1]*TT[,2], lwd = 3)
for(i in 1:nsim) z1[i,] <- zs$z[i,1] + zs$z[i,2]*xp
ci1 <- HPDinterval(as.mcmc(z1)) 
lines(0:100, ci1[,1])
lines(0:100, ci1[,2])

plot(fd, mean_b[seq(2, (ns*2)+1 ,by=2),1], pch = 19, col=rgb(0,0,0, 0.4),
     ylab = "", xlab = "", main = "")
lines(fd, mean_z[2,1] + mean_z[4,1]*TT[,2], lwd = 3)
for(i in 1:nsim) z2[i,] <- zs$z[i,4] + zs$z[i,5]*xp
ci2 <- HPDinterval(as.mcmc(z2)) 
lines(0:100, ci2[,1])
lines(0:100, ci2[,2])

#for(i in 1:100) lines(fd, zs$z[i,4] + zs$z[i,5]*TT[,2], col=rgb(0,0,0, 0.2) )
plot(bs, mean_b[seq(1, (ns*2),by=2),1], pch = 19, col=rgb(0,0,0, 0.4),
     ylab = "", xlab = "", main = "")
lines(bs, mean_z[1,1] + mean_z[5,1]*TT[,3], lwd = 3)
for(i in 1:nsim) z3[i,] <- zs$z[i,1] + zs$z[i,3]*xg
ci3 <- HPDinterval(as.mcmc(z3)) 
lines(seq(min(bs), max(bs), length.out = 100), ci3[,1])
lines(seq(min(bs), max(bs), length.out = 100), ci3[,2])


plot(bs, mean_b[seq(2, (ns*2)+1 ,by=2),1], pch = 19, col=rgb(0,0,0, 0.4),
     ylab = "", xlab = "",  main = "")

for(i in 1:nsim) z4[i,] <- zs$z[i,4] + zs$z[i,6]*xg
ci4 <- HPDinterval(as.mcmc(z4)) 
lines(seq(min(bs), max(bs), length.out = 100), ci4[,1])
lines(seq(min(bs), max(bs), length.out = 100), ci4[,2])

par(op)
par(mfrow=c(1,1))
dev.off()

#----
# prepare tables
gens <- sort(unique(datos$Phylo_name))
slopes <- mean_b[seq(2, (ns*2)+1 ,by=2),]
rownames(slopes) = gens
write.csv(slopes, file = "sp-slopes.csv")

intercepts <-  mean_b[seq(1, (ns*2),by=2),]
rownames(intercepts) = gens
write.csv(intercepts, file = "sp-intercepts.csv")

write.csv(mean_z, file = "sp-zetas.csv")

eta <- fit_summary[grepl("etasq", rownames(fit_summary)),]
delta <- fit_summary[grepl("delta", rownames(fit_summary)),]
rhos <- fit_summary[grepl("rhosq", rownames(fit_summary)),]

write.csv(rbind(eta, delta, rhos), file = "sp-etas.csv") 

rho <- fit_summary[grepl("rho", rownames(fit_summary)),]

sigma <- fit_summary[grepl("Sigma", rownames(fit_summary)),]
write.csv(rbind(sigma, rho[1,]), file = "sp-sigma.csv") 

#-----

# histogram of slopes  
tiff(filename="sp-hist_slopes.tif", res=res, unit="cm", height=8, width=8)
op <- par(las = 1, bty = "n", cex = 1)
hist(mean_b[seq(2, (ns*2)+1 ,by=2),1], 60, 
     main = "", xlab = "slope", col=rgb(0,0,0, 0.5))
par(op)
dev.off()


# plot response against lipids -------------------------------------------------

mz <- matrix(mean_z[,1], 3, 2, byrow = TRUE)

xl <- seq(-0.9, 2.4, by = 0.1 )
fro <- seq(0, 100, by = 20)
fr <- (fro - mean(fd))/ sd(fd)

xlt <- xl * sd(datos$lipid) + mean(datos$lipid)
# grs <- viridis(length(fr))

rwb <- colorRampPalette(colors = c("grey", "blue"))
grs <- rwb(length(fr))
#ev <- matrix(NA, length(xl), length(fr))

tiff(filename="sp-lipids.tif", res=res, unit="cm", height=10, width=11)
op <- par(las = 1, bty = "n", cex = 1)
plot(xlt, rep(-1, length(xl)), ylim = c(0,5),
     xlab = "lipids", ylab = "expected visits")
legend("bottomright", legend = paste( fro, "%"), 
       col = grs, pch = 19, bty = "n", cex = 0.6)
for(i in 1:length(fr)){
  lines(xlt, exp( mz[1,1] + mz[2,1] * fr[i]  + (mz[2,1] + mz[2,2] *fr[i]) * xl), 
        lwd = 3, col = grs[i])
  #ev[,i] = exp( mz[1,1] + mz[2,1] * fr[i]  + (mz[2,1] + mz[2,2] *fr[i]) * xl)
}
par(op)
dev.off()
#--------------