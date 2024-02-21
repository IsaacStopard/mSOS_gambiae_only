rm(list = ls())

# R script to run the VC_functions.R on the MRC Centre cluster
# Author: Isaac J Stopard
# Version: 1.0.0
# Notes: this will only work with the old version of the MRC Centre cluster with access to the MRC centre network

library(tidyverse); library(rrq); library(didehpc)

setwd("/Volumes/ijs11/mSOS_gambiae_only")

##################################
##### running on the cluster #####
##################################

# some randomly sampled parameter estimates from the Stan model fits to estimate the EIP and biting rate
params_df <- readRDS(file = "params_df.rds") %>% as.data.frame()

pars <- subset(params_df, a0 == 0)

# options(didehpc.template = "GeneralNodes",
#         didehpc.cores = 1,
#         didehpc.home = "/Volumes/ijs11",
#         didehpc.username = "ijs11")

config <- didehpc::didehpc_config(credentials = "ijs11",
                                  use_rrq = TRUE,
                                  worker_timeout = 60 * 60 * 1,
                                  cluster = "wpia-hn", #fi--dideclusthn #fi--didemrchnb
                                  cores = 1,
                                  home = "/Volumes/ijs11")

packages <- list(loaded = c("zipfR", "calculus", "cubature", "rrq"))
sources <- "utils/VC_functions.R"

root <- "contexts"
ctx <- context::context_save(root, sources = sources, packages = packages)
obj <- didehpc::queue_didehpc(ctx, config = config)

obj$install_packages(packages = "mrc-ide/rrq")

rrq <- obj$rrq_controller()

workers <- obj$submit_workers(150)

# t <- obj$enqueue_bulk(pars[,c("a0", "alpha", "g1", "g2",
#                                          "shape", "rate", "mu_PL", "k")],
#                       E_j_bites_a0)
# 
# obj$stop_workers()

n_p <- ceiling(nrow(pars)/150)
n_i <- ceiling(nrow(pars) / n_p)
inds <- seq(1, nrow(pars))
s_inds <- split(inds, ceiling(seq_along(inds)/n_i))

t_out <- vector(mode = "list", length = n_p)

for(i in 1:n_p){
  print(i)
  t <- rrq$enqueue_bulk(pars[s_inds[[i]],c("a0", "alpha", "g1", "g2",
                                           "shape", "rate", "mu_PL", "k")],
                        E_j_bites_a0)
  
  t_out[[i]] <- t
  saveRDS(t, file = paste0("results/E_j_a0_0_only_",i,".rds"))
  
}

# checking the values
pars[,"E_j_a0"] <- rep(NA, nrow(pars))
pars[,"E_j_a0_conv"] <- rep(NA, nrow(pars))

for(i in 1:n_p){
  print(i)
  t_out[[i]] <- readRDS(file = paste0("results/E_j_a0_0_only_",i,".rds"))
  l_t <- length(t_out[[i]])
  values <- rep(NA, l_t)
  conv <- rep(NA, l_t)
  
  for(j in 1:l_t){
    values[j] <- t_out[[i]][[j]]$value
    conv[j] <- t_out[[i]][[j]]$cuba$returnCode
  }
  
  pars[s_inds[[i]], "E_j_a0"] <- values
  pars[s_inds[[i]], "E_j_a0_conv"] <- conv
}

unique(pars$E_j_a0_conv)

pars %>% group_by(location, temp_source, temp, feed) %>% summarise(med = median(E_j_a0),
                                                                   low = quantile(E_j_a0, probs = c(0.025))[[1]],
                                                                   up = quantile(E_j_a0, probs = c(0.975))[[1]])


pars_r <- subset(params_df, a0 != 0)
n_p_r <- ceiling(nrow(pars_r)/450)
n_i_r <- ceiling(nrow(pars_r) / n_p_r)
inds_r <- seq(1, nrow(pars_r))
s_inds_r <- split(inds_r, ceiling(seq_along(inds_r)/n_i_r))

t_out_r <- vector(mode = "list", length = n_p_r)

for(i in 12:n_p_r){
  print(i)
  t_r <- rrq$enqueue_bulk(pars_r[s_inds_r[[i]],c("a0", "alpha", "g1", "g2",
                                                 "shape", "rate", "mu_PL", "k")],
                          E_j_bites_a0)
  
  t_out_r[[i]] <- t_r
  saveRDS(t_r, file = paste0("results/E_j_a0_no_0_",i,".rds"))
  
}

rrq$worker_stop()


rrq$destroy()

workers <- obj$submit_workers(150)

res_fixed <- rrq$enqueue_bulk(pars[,c("a0", "alpha", "g1", "g2", #obj
                                      "shape", "rate", "mu_PL", "k", "EIP")], 
                              E_j_bites_a0_fixed)


rrq$message_send("TIMEOUT_SET", 0)

t$results()

t <- obj$enqueue_bulk(pars[,c("a0", "alpha", "g1", "g2", 
                              "shape", "rate", "mu_PL", "k", "EIP")], 
                      E_j_bites_a0_fixed)

t <- obj$enqueue(sessionInfo())