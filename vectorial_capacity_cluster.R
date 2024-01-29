rm(list = ls())

library(tidyverse); library(rrq); library(didehpc)

setwd("/Volumes/ijs11/mSOS_gambiae_only")

##################################
##### running on the cluster #####
##################################

# data
params_df <- readRDS(file = "params_df.rds")

pars <- params_df %>% rename(shape = shape_total_S,
                                    rate = rate_total_S,
                                    mu_PL = mu,
                             EIP = EIP_DD
                                    ) %>% as.data.frame()

config <- didehpc::didehpc_config(credentials = "ijs11",
                                  use_rrq = TRUE,
                                  worker_timeout = 60 * 60 * 10,
                                  cluster = "fi--dideclusthn",
                                  cores = 1,
                                  home = "/Volumes/ijs11")

packages <- list(loaded = c("zipfR", "calculus", "cubature", "rrq"))
sources <- "utils/VC_functions.R"

root <- "contexts"
ctx <- context::context_save(root, sources = sources, packages = packages)
obj <- didehpc::queue_didehpc(ctx, config = config)

obj$install_packages("mrc-ide/rrq")

rrq <- obj$rrq_controller()

workers <- obj$submit_workers(150)

n_p <- 10
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
  saveRDS(t, file = paste0("E_j_a0_",i,".rds"))
}

rrq$worker_stop()

obj$stop_workers()

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






