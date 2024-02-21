# R script to process the parameters for the individual vectorial capacity estimates and plotting the results
# Author: Isaac J Stopard
# Version: 1.0.0
# Notes: # estimating the temp-dependent vectorial capacity with degree-day model and mSOS model comparison

rm(list = ls())
library(tidyverse);
library(cowplot); library(zipfR);
library(truncnorm); library(ggpmisc); library(patchwork);
library(doParallel); library(foreach)
library(calculus);

theme_set(theme_bw() + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  text = element_text(size = 18)))

#####################
##### functions #####
#####################

source(file = "utils/read_data.R")
source(file = "utils/VC_functions.R")

###########################
### reading in the data ###
###########################

# EIP Stan fits
fit_temp <- readRDS("fits/fit_mSOS_temp_only_f2_f3.rds")

# biting rate Stan fits
fit_br <- readRDS("fits/biting_rate_fit.rds")

params_temp <- rstan::extract(fit_temp)

params_br <- rstan::extract(fit_br)

# mean_biting_rate <- read.csv(file = "biting_rate_estimates.csv")

params_df <- rbind(data.frame(location = c("Kericho", "Kitale", "Kisumu", "Garissa"),
                     "recent" = c(17.5, 19, 23.4, 28.8),
                     "future" = c(19.5, 20.9, 25.2, 30)) %>% pivot_longer(cols = c("recent", "future"),
                                                                          names_to = "temp_source",
                                                                          values_to = "temp") %>% 
                     mutate(scaled_temp = (temp - all_data$m_temp) / all_data$sd_temp,
                            g1 = quadratic_function(t = temp, a = 0.000233, b = -0.0101, c = 0.118),
                            g2 = 0.00809 * temp - 0.111, feed = "feed_1"),
                   data.frame(location = c("Kericho", "Kitale", "Kisumu", "Garissa"),
                        "recent" = c(17.5, 19, 23.4, 28.8),
                        "future" = c(19.5, 20.9, 25.2, 30)) %>% pivot_longer(cols = c("recent", "future"),
                                                                             names_to = "temp_source",
                                                                             values_to = "temp") %>% mutate(
                                                                               scaled_temp = (temp - all_data$m_temp) / all_data$sd_temp,
                                                                               g1 = quadratic_function(t = temp, 
                                                                                                       a = 0.00005632, b = -0.00244, c = 0.028),
                                                                               g2 = 0.00851 * temp - 0.0916, feed = "feed_2"))


#params_df[,"alpha"] <- mean_biting_rate[match(params_df$temp, mean_biting_rate$temp), "med"]

extract_params <- function(i, params_df, params_temp, params_br, n_s){
  set.seed(123)
  i_pdf <- sample(seq(1, 10000), size = n_s, replace = FALSE)
  set.seed(123)
  i_br <- sample(seq(1, 5000), size = n_s, replace = FALSE)
  
  scaled_temp <- params_df[i, "scaled_temp"][[1]]
  temp <- params_df[i, "temp"][[1]]
  
  shape_S <- quadratic_function(scaled_temp, params_temp$a_shape_S, params_temp$b_shape_S, params_temp$c_shape_S) 
  rate_S <- scaled_temp * params_temp$m_rate_S + params_temp$c_rate_S
  shape_O <- params_temp$shape_O
  rate_O <- params_temp$rate_O
  mu_total_S <- (rate_O * shape_S + rate_S * shape_O) / (rate_O * rate_S)
  sigma_sq_S <- (rate_O^2 * shape_S + rate_S^2 * shape_O) / (rate_O^2 * rate_S^2)
  shape_total_S <- mu_total_S^2 / sigma_sq_S
  rate_total_S <- mu_total_S / sigma_sq_S
  mu <- logistic_function(quadratic_function(scaled_temp, params_temp$a_mu, params_temp$b_mu, params_temp$c_mu)) * params_temp$g_mu
  k <- params_temp$k
  
  Briere_out <- ifelse(temp >= params_br$T0 & temp <= params_br$Tm, params_br$a/1000 * temp * (temp - params_br$T0) * sqrt((params_br$Tm - temp)), 0)
  
  
  return(params_df[rep(i, length(i_pdf)), ] %>% 
           mutate(shape_total_S = shape_total_S[i_pdf],
                    rate_total_S = rate_total_S[i_pdf],
                    mu = mu[i_pdf],
                    k = k[i_pdf],
                    alpha = Briere_out[i_br],
                  i_pdf = i_pdf,
                  i_br = i_br)
  )
}

params_df_cluster <- lapply(seq(1, nrow(params_df)), extract_params, params_df = params_df, params_temp = params_temp, params_br = params_br, n_s = 100) %>% bind_rows()

n <- nrow(params_df_cluster)
a0 <- seq(0, 15, 1)

params_df_cluster <- params_df_cluster %>% slice(rep(1:n(), each = length(a0))) %>% mutate(a0 =rep(a0, n))

params_df_cluster <- params_df_cluster %>% mutate(EIP_DD = degree_day(temp))

params_df_cluster <- params_df_cluster %>% rename(shape = shape_total_S,
                                  rate = rate_total_S,
                                  mu_PL = mu,
                                  EIP = EIP_DD)

saveRDS(params_df_cluster, file = "/Volumes/ijs11/mSOS_gambiae_only/params_df.rds")
saveRDS(params_df_cluster, file = "fits/params_df.rds")


#################################################################################
##### calculating for vectorial capacity for the median parameter estimates #####
#################################################################################


# the vectorial capacity was calculated on the MRC Centre cluster, which requires access.
# run locally use:
# params_df <- params_df_cluster %>% mutate(E_j_a0 = v_E_j_bites_a0(a0 = a0,
#                                                           alpha = alpha,
#                                                           g1 = g1,
#                                                           g2 = g2,
#                                                           shape = shape_total_S,
#                                                           rate = rate_total_S,
#                                                           mu_PL = mu,
#                                                           k = k),
#                                   EIP_DD = degree_day(temp),
#                                   E_j_a0_fixed = v_E_j_bites_a0_fixed(a0 = a0,
#                                                                       alpha = alpha,
#                                                                       g1 = g1,
#                                                                       g2 = g2,
#                                                                       EIP = EIP_DD))
# 
# saveRDS(params_df, file = "params_df_E_j_bites_a0.rds")
# 
# 
# check_df <- params_df[1,]
# 
# saveRDS(params_df, file = "params_df_E_j_bites_a0.rds")
# 
# params_df <- readRDS(file = "params_df_E_j_bites_a0.rds")

##################################
##### processing the results #####
##################################

# reading in the results from the cluster runs

pars <- readRDS(file = "results_temp/params_df_results_EIP_mSOS.rds")

pars <- left_join(pars, subset(params_df_cluster, a0 <= 15), 
                  by = c("location", "temp_source", "temp", "scaled_temp", "g1", "g2", "feed", "shape", "rate", "mu_PL", "k", "alpha", "a0", "EIP"))

pars_sum <- pars %>% group_by(location, temp_source, temp, feed, a0) %>% summarise(med = median(E_j_a0),
                                                                   low = quantile(E_j_a0, probs = c(0.025))[[1]],
                                                                   up = quantile(E_j_a0, probs = c(0.975))[[1]],
                                                                   med_dd = median(E_j_dd),
                                                                   low_dd = quantile(E_j_dd, probs = c(0.025))[[1]],
                                                                   up_dd = quantile(E_j_dd, probs = c(0.975))[[1]]) %>% as.data.frame()



pars_sum$location <- factor(pars_sum$location, levels = c("Kericho", "Kitale", "Kisumu", "Garissa"))

pars_sum$E_j_a0 <- paste0(round(pars_sum$med, digits = 2)," (",round(pars_sum$low, digits = 2),"-",round(pars_sum$up, digits = 2),")")
pars_sum$E_j_a0_fixed <- paste0(round(pars_sum$med_dd, digits = 2)," (",round(pars_sum$low_dd, digits = 2),"-",round(pars_sum$up_dd, digits = 2),")")

pars_sum_0 <- subset(pars_sum, a0 == 0)

relative_E_j_a0_all <- left_join(
  subset(pars, temp_source == "recent")[,c("location", "feed", "a0", "i_pdf", "i_br", "E_j_a0", "E_j_dd")] %>% rename(r_E_j_a0 = E_j_a0, r_E_j_dd = E_j_dd),
  
  subset(pars, temp_source == "future")[,c("location", "feed", "a0", "i_pdf", "i_br", "E_j_a0", "E_j_dd")] %>% rename(f_E_j_a0 = E_j_a0, f_E_j_dd = E_j_dd),
                             by = c("location", "a0", "feed", "i_pdf", "i_br")) %>% 
  mutate(E_j_a0_ratio = f_E_j_a0 / r_E_j_a0,
         E_j_dd_ratio = f_E_j_dd / r_E_j_dd)
  
relative_E_j_a0 <- relative_E_j_a0_all %>% group_by(location, a0, feed) %>% 
  summarise(med = median(E_j_a0_ratio),
            med_dd = median(E_j_dd_ratio),
            low = quantile(E_j_a0_ratio, probs = c(0.025)),
            low_dd = quantile(E_j_dd_ratio, probs = c(0.025)),
            up = quantile(E_j_a0_ratio, probs = c(0.975)),
            up_dd = quantile(E_j_dd_ratio, probs = c(0.975)))

relative_E_j_a0 <- rbind(relative_E_j_a0[,c("location", "feed", "a0", "med", "low", "up")] %>% mutate(model = "Suh-Stopard"),
                         relative_E_j_a0[,c("location", "feed", "a0", "med_dd", "low_dd", "up_dd")] %>% rename(med = med_dd, low = low_dd, up = up_dd) %>% mutate(model = "Degree-day"))

# viewing the results
subset(relative_E_j_a0, a0 == 0) %>% 
  mutate(result = paste0(round(med, digits = 2), " (", round(low, digits = 2), " – ", round(up, digits = 2),")"))

# changing the order of the locations
relative_E_j_a0$location <- factor(relative_E_j_a0$location, levels = c("Kericho", "Kitale", "Kisumu", "Garissa"))

# rbind(params_df %>% dplyr::select("location", "temp_source", "feed", "a0", "E_j_a0") %>% tidyr::pivot_wider(names_from = temp_source,
#                                                                                                 values_from = E_j_a0) %>% mutate(model = "EIP_dist"),
#                          params_df %>% dplyr::select("location", "temp_source", "feed", "a0", "E_j_a0_fixed") %>% tidyr::pivot_wider(names_from = temp_source,
#                                                                                                                        values_from = E_j_a0_fixed) %>% 
#                            mutate(model = "degree_day")) %>% mutate(ratio = future / recent)

params_df <- params_df %>% 
  mutate(temp_source = ifelse(temp_source == "future", "Future mean temperature", "Recent mean temperature"))

params_df$temp_source <- factor(params_df$temp_source, levels = c("Recent mean temperature", "Future mean temperature"))


vc_plot_df <- left_join(left_join(pars_sum[,c("location", "temp_source", "temp", "feed", "a0", "med", "med_dd")] %>% pivot_longer(cols = c("med", "med_dd"), 
                                         names_to = "model", values_to = "VC_med") %>% mutate(model = ifelse(grepl("dd", model), "fixed", "dist")),
                                         
                                   pars_sum[,c("location", "temp_source", "temp", "feed", "a0", "low", "low_dd")] %>% pivot_longer(cols = c("low", "low_dd"), 
                                         names_to = "model", values_to = "VC_low") %>% mutate(model = ifelse(grepl("dd", model), "fixed", "dist")),
                                   by = c("location", "temp_source", "temp", "feed", "a0", "model")),
                                   
                                   pars_sum[,c("location", "temp_source", "temp", "feed", "a0", "up", "up_dd")] %>% 
                                     pivot_longer(cols = c("up", "up_dd"), names_to = "model", values_to = "VC_up") %>% 
                                     mutate(model = ifelse(grepl("dd", model), "fixed", "dist")),
                                   by = c("location", "temp_source", "temp", "feed", "a0", "model")) %>% 
                    mutate(temp_source = ifelse(temp_source == "future", "Future mean temperature", "Recent mean temperature"))

vc_plot_df$temp_source <- factor(vc_plot_df$temp_source, levels = c("Recent mean temperature", "Future mean temperature")) 

vc_plot <- ggplot(data = vc_plot_df,
                  aes(x = a0, y = VC_med, ymin = VC_low, ymax = VC_up, 
                      group = interaction(feed, model, temp_source, location))) + 
  geom_ribbon(aes(fill = model), alpha = 0.1) +
  #geom_point(aes(fill = model, shape = feed), size = 3, alpha = 0.5) +
  geom_line(aes(colour = model, linetype = feed), linewidth = 1) + 
  facet_grid(vars(location), vars(temp_source)) +
  ylab("Expected number of infectious bites\nper infected mosquito (z)") + 
  xlab("Mosquito age when infected (days)") +
  #scale_shape_manual(values = c(21, 22), labels = c("Feed1", "Feed2"), name = "Feed") +
  scale_linetype_manual(values = c(1, 2), name = "Feed", labels = c("Feed1", "Feed2")) +
  scale_colour_manual(values = c("#009E73", "grey40"), name = "EIP model", labels = c("Suh-Stopard", "Degree-day")) +
  scale_fill_manual(values = c("#009E73", "grey40"), name = "EIP model", labels = c("Suh-Stopard", "Degree-day"))

leg <- cowplot::get_legend(vc_plot)

r_plot <- ggplot(data = relative_E_j_a0, 
       aes(x = a0, y = med, ymin = low, ymax = up, group = interaction(model, feed, location))) +
  geom_ribbon(aes(fill = model), alpha = 0.1) +
  #geom_point(aes(fill = model, shape = feed), size = 3, alpha = 0.5) +
  geom_line(aes(colour = model, linetype = feed), linewidth = 0.75) +
  facet_wrap(vars(location), scales = "free_y") +
  scale_y_sqrt() +
  ylab("Ratio of future:recent expected infectious\nbites per infected mosquito") +
  xlab("Mosquito age when infected (days)") +
  scale_linetype_manual(values = c(1, 2), labels = c("Feed1", "Feed2"), name = "Feed") +
  scale_colour_manual(values = c("grey40", "#009E73"), name = "EIP model", labels = c("Degree-day", "Suh-Stopard")) +
  scale_fill_manual(values = c("grey40", "#009E73"), name = "EIP model", labels = c("Degree-day", "Suh-Stopard"))  

vc_plot_nl <- vc_plot + theme(legend.position = "none")
r_plot_nl <- r_plot + theme(legend.position = "none")

png(file = "results_temp/VC_plots.png", height = 1050, width = 975)
((vc_plot_nl / r_plot_nl) | leg) +
  plot_layout(widths = c(0.875, 0.175)) +
  plot_annotation(tag_levels = list(c("a", "b"), "")) &
  theme(plot.tag = element_text(face = 'bold'))
dev.off()

# checking the values
cbind(subset(params_df, a0 == 0 & temp_source == "recent" & feed == "feed_1")$E_j_a0 %>% round(digits = 2),
           subset(params_df, a0 == 0 & temp_source == "recent" & feed == "feed_1")$E_j_a0_fixed %>% round(digits = 2),
           subset(params_df, a0 == 0 & temp_source == "recent" & feed == "feed_2")$E_j_a0 %>% round(digits = 2),
           subset(params_df, a0 == 0 & temp_source == "recent" & feed == "feed_2")$E_j_a0_fixed %>% round(digits = 2)) %>% view()

cbind(subset(relative_E_j_a0, feed == "feed_1" & model == "EIP_dist" & a0 == 0)$ratio %>% round(digits = 2), 
      subset(relative_E_j_a0, feed == "feed_1" & model == "degree_day" & a0 == 0)$ratio %>% round(digits = 2),
      subset(relative_E_j_a0, feed == "feed_2" & model == "EIP_dist" & a0 == 0)$ratio %>% round(digits = 2), 
      subset(relative_E_j_a0, feed == "feed_2" & model == "degree_day" & a0 == 0)$ratio %>% round(digits = 2)) %>% view()



