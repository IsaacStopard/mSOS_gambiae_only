rm(list = ls())
library(tidyverse)
source(file = "utils/data_wrangling_functions.R")

s_data <- read.csv(file = "data/ES_new_constant_temp_spz_data.csv")
o_i_data <- read.csv(file = "data/ES_new_constant_temp_oocyst_intensity_infected.csv")
o_p_data <- read.csv(file = "data/ES_new_constant_temp_oocyst_prevalence.csv")

# sporozoite data
s_data <- process_prevalence_data(s_data)
write.csv(s_data, file = "data/processed/ES_new_constant_temp_spz_processed.csv")

# oocyst data
o_p_data <- process_prevalence_data(o_p_data)
o_data <- o_p_data[-which(o_p_data$presence == 1),] # removing infected
o_data$Oocyst_number <- rep(0, nrow(o_data))
colnames(o_i_data) <- c("DPI", "Cup", "Experiment", "temp", "DTR", "Oocyst_number", "presence", "gametocytemia")
o_data <- rbind(o_data, o_i_data)

write.csv(o_data, file = "data/processed/ES_new_constant_temp_oocyst_processed.csv")