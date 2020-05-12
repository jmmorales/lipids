library(rstan) # observe startup messages
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)