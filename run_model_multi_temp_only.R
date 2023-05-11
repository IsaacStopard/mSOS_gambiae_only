rm(list = ls())
library(tidyverse); library(rstan); library(shinystan); library(cowplot); library(zipfR); library(truncnorm);library(ggpmisc);
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source(file = "utils/functions_temp_only.R")
###########################
### reading in the data ###
###########################

s_data <- read.csv(file = "data/processed/ES_new_constant_temp_spz_processed.csv")
o_data <- read.csv(file = "data/processed/ES_new_constant_temp_oocyst_processed.csv")

s_data <- s_data[,c(2:ncol(s_data))]
s_data$index_ref <- rep(1, nrow(s_data))
s_data$gametocytemia <- round(s_data$gametocytemia, digits = 5)

# wrangling the oocyst data
o_data <- o_data[,c(2:ncol(o_data))]
o_data$gametocytemia <- round(o_data$gametocytemia, digits = 5)
o_data$index_ref <- rep(1, nrow(o_data))

s_data <- s_data[-which(s_data$temp == 17 & s_data$gametocytemia == 0.00024),]
o_data <- o_data[-which(o_data$temp == 17 & o_data$gametocytemia == 0.00024),]

all_data <- temp_g_indexing(o_data, s_data)
#saveRDS(all_data, file = "data/processed/all_data_processed.rds")
s_data_in <- generate_prevalence(all_data$sporozoite_data)
o_data_in <- oocyst_intensity_indexing(all_data$oocyst_data)

# values from previous model fitting to Anopheles stephensi data
mean_temp <- 27.9032
sd_temp <-  3.471223
m_rate_O <- 1.44
c_rate_O <- 4.19

# don't extrapolate beyond 21 degrees celsius
scaled_temp <- (all_data$unique_temp - mean_temp) / sd_temp
rate_O_prior <- scaled_temp * m_rate_O + c_rate_O
rate_O_prior[which(rate_O_prior < 1)] <- rep(rate_O_prior[5], length(which(rate_O_prior < 1)))
rate_O_prior <- rep(max(rate_O_prior), length(rate_O_prior))
# posterior predictive distribution times
PPD_times <- sort(unique(c(seq(0,49,0.5), s_data_in$DPI, o_data_in$unique_oocyst_intensity$DPI)))
length_ppd_times <- length(PPD_times)

PPD_times_O <- sort(unique(c(seq(0,34.0,0.5), o_data_in$unique_oocyst_intensity$DPI)))
length_ppd_times_O <- length(PPD_times_O)

iterations <- 5500
warmup <- 3000
chains <- 4

s_data_in_temp <- generate_prevalence_temp(subset(all_data$sporozoite_data, index_g!=1))
o_data_in_temp <- oocyst_intensity_indexing_temp(subset(all_data$oocyst_data, index_g!=1))

# pooled model
data_in_temp <- list(
  n_S = as.integer(nrow(s_data_in_temp)), 
  S_sample = as.integer(s_data_in_temp$sample),
  S_positive = as.integer(s_data_in_temp$positive), 
  S_time = as.double(s_data_in_temp$DPI),
  S_ind_temp = as.integer(s_data_in_temp$index_temp),
  
  n_O_intsy = as.integer(length(o_data_in_temp$oocyst_intensity_index)),
  O_intsy_index = as.integer(o_data_in_temp$oocyst_intensity_index),
  
  n_unq_O_intsy_index = as.integer(nrow(o_data_in_temp$unique_oocyst_intensity)),
  unq_O_intsy_time = as.double(o_data_in_temp$unique_oocyst_intensity[,"DPI"]),
  unq_O_intsy = as.double(o_data_in_temp$unique_oocyst_intensity[,"Oocyst_number"]),
  unq_O_intsy_ind_temp = as.integer(o_data_in_temp$unique_oocyst_intensity[,"index_temp"]),
  
  n_ppd_times = as.integer(length_ppd_times), 
  PPD_times = PPD_times,
  n_ppd_times_O = as.integer(length_ppd_times_O), 
  PPD_times_O = PPD_times_O,
  
  n_unq_temp = as.integer(length(all_data$unique_temp)),
  temp = as.double(all_data$unique_temp_scaled),
  rate_O_prior = rate_O_prior[1]  #rate_O_prior[all_data$unique_gt[,"index_temp"]]
)

finit_temp <- function(){
  list(shape_O = 10, #rep(10, nrow(all_data$unique_gt)), 
       rate_O = 4, # rep(4, nrow(all_data$unique_gt)),
       a_shape_S = 0,
       b_shape_S = 0,
       c_shape_S = 20,
       m_rate_S = 0, 
       c_rate_S = 5,
       a_delta = 0,
       b_delta = 0,
       c_delta = -1.1,
       a_delta_S = 0,
       b_delta_S = 0,
       c_delta_S = -1.1,
       a_mu = 1.5, 
       b_mu = 1.5,
       c_mu = 1.5,
       g_mu = 8,
       k = 3 # rep(3, length(all_data$unique_temp)) # parasite load parameters
  )
}

fit_temp <- stan(file = "model/mSOS_multi_8.stan", data = data_in_temp, iter=iterations, chains = chains, seed=12345,
                 warmup = warmup, init = finit_temp, control = list(max_treedepth = 12.5, adapt_delta = 0.99))

saveRDS(fit_temp, "fits/fit_mSOS_temp_only_f2_f3.rds")
fit_temp <- readRDS("fits/fit_mSOS_temp_only_f2_f3.rds")

p_params <- rstan::extract(fit_temp, c("shape_O", "rate_O", "k",
                                       "a_shape_S", "b_shape_S", "c_shape_S",
                                       "m_rate_S", "c_rate_S",
                           "a_delta", "b_delta", "c_delta",
                           "a_delta_S", "b_delta_S", "c_delta_S", 
                           "a_mu", "b_mu", "c_mu", "g_mu"))

bind_rows(lapply(seq(1, length(p_params)), function(i){
  p <- p_params[[i]]
  c(names(p_params[i]), round(c(median(p), quantile(p, probs = c(0.025, 0.975))), digits = 2))
}))
# independent model

data_in_all_temp <- list(
  n_S = as.integer(nrow(s_data_in_temp)), 
  S_sample = as.integer(s_data_in_temp$sample),
  S_positive = as.integer(s_data_in_temp$positive), 
  S_time = as.double(s_data_in_temp$DPI),
  S_ind_gt = as.integer(s_data_in_temp$index_temp),
  
  n_O_intsy = as.integer(length(o_data_in_temp$oocyst_intensity_index)),
  O_intsy_index = as.integer(o_data_in_temp$oocyst_intensity_index),
  
  n_unq_O_intsy_index = as.integer(nrow(o_data_in_temp$unique_oocyst_intensity)),
  unq_O_intsy_time = as.double(o_data_in_temp$unique_oocyst_intensity[,"DPI"]),
  unq_O_intsy = as.double(o_data_in_temp$unique_oocyst_intensity[,"Oocyst_number"]),
  unq_O_intsy_ind_gt = as.integer(o_data_in_temp$unique_oocyst_intensity[,"index_temp"]),
  
  n_ppd_times = as.integer(length_ppd_times), 
  PPD_times = PPD_times,
  n_ppd_times_O = as.integer(length_ppd_times_O), 
  PPD_times_O = PPD_times_O,
  
  n_unq_gt = as.integer(length(all_data$unique_temp)),
  rate_O_prior = rate_O_prior[1]  #rate_O_prior[all_data$unique_gt[,"index_temp"]]
)

finit_all <- function(){
  list(shape_O = rep(10, length(all_data$unique_temp)), 
       rate_O = rep(4, length(all_data$unique_temp)),
       shape_S = rep(20, length(all_data$unique_temp)), 
       rate_S = rep(3, length(all_data$unique_temp)), # sporozoite development parameters
       mu = rep(3, length(all_data$unique_temp)),
       k = rep(3, length(all_data$unique_temp)), # parasite load parameters
       delta = rep(0.75, length(all_data$unique_temp)),
       delta_S = rep(0.9, length(all_data$unique_temp))
  )
}

fit_all_temp <- stan(file = "model/mSOS_multi.stan", data = data_in_all_temp, iter=iterations, chains = chains, seed=12345,
                warmup = warmup, init = finit_all, control = list(max_treedepth = 12.5, adapt_delta = 0.99))

saveRDS(fit_all_temp, file = "fits/fit_mSOS_multi_temp.rds")
fit_all_temp <- readRDS(file = "fits/fit_mSOS_multi_temp.rds")

# model checking
pars <- c("shape_O", "rate_O", "a_shape_S", "b_shape_S", "c_shape_S",
          "m_rate_S", "c_rate_S",
          "a_delta", "b_delta", "c_delta", 
          "a_delta_S", "b_delta_S", "c_delta_S", 
          "a_mu", "b_mu", "c_mu",
          "k")

png("results_temp/pairs_plot.png", height = 2500, width = 2500)
pairs(fit_temp, pars = pars)
dev.off()
png("results_temp/traceplot.png", height = 750, width = 1000)
traceplot(fit_temp, pars = pars)
dev.off()
