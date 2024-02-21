# R script to read in the mosquito dissection data
# Author: Isaac J Stopard
# Version: 1.0.0

###########################
### reading in the data ###
###########################
source("utils/functions_temp_only.R")

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
s_data_in <- generate_prevalence_temp(all_data$sporozoite_data) %>% mutate(temp_label = paste0(temp, "°C")) # includes the pilot feed data
o_data_in <- oocyst_intensity_indexing(all_data$oocyst_data)
  
# values from previous model fitting to Anopheles stephensi data
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008658

mean_temp <- 27.9032
sd_temp <-  3.471223
m_rate_O <- 1.44
c_rate_O <- 4.19
  
# number of MCMC samples
iterations <- 5500
warmup <- 3000
chains <- 4
  
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
  
s_data_in_temp <- generate_prevalence_temp(subset(all_data$sporozoite_data, index_g!=1)) %>% 
    mutate(temp_label = paste0(temp, "°C"))
  
o_data_in_temp <- oocyst_intensity_indexing_temp(subset(all_data$oocyst_data, index_g!=1))


