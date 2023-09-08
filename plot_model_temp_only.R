rm(list = ls())

# packages
library(tidyverse); library(rstan); library(shinystan); library(cowplot); 
library(zipfR); library(truncnorm);library(ggpmisc); library(patchwork); library(DescTools);
library(calculus)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# functions and data that are used in multiple files
source(file = "utils/functions_temp_only.R")
source(file = "utils/plotting_functions_temp_only.R")

# changing the default ggplot theme
theme_set(theme_bw() + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  text = element_text(size = 18)))

# reading in the EIP model fits
fit_temp <- readRDS("fits/fit_mSOS_temp_only_f2_f3.rds")
fit_all_temp <- readRDS(file = "fits/fit_mSOS_multi_temp.rds")

# model checking
pars <- c("shape_O", "rate_O", "a_shape_S", "b_shape_S", "c_shape_S",
          "m_rate_S", "c_rate_S",
          "a_delta", "b_delta", "c_delta", 
          "a_delta_S", "b_delta_S", "c_delta_S", 
          "a_mu", "b_mu", "c_mu",
          "k")

# model diagnostic plots
png(file = "traceplot_pooled_model.png", height = 500, width = 1000)
traceplot(fit_temp, pars = pars)
dev.off()

##############################
### visualising model fits ###
##############################

### posterior predictions
### sporozoites
# extracting the posterior predicted sporozoite prevalence from the MCMC samples
S_ppd_df <- rbind(run_prop_ppd_df(fit_temp, "pooled", "S_prevalence_ppd", length_ppd_times, PPD_times),
                  run_prop_ppd_df(fit_all_temp, "independent", "S_prevalence_ppd", length_ppd_times, PPD_times)) %>% 
  mutate(temp_label = paste0(temp, "°C"))

S_plot <- ggplot(data = S_ppd_df) + 
  geom_ribbon(aes(x = DPI, ymin = lower, ymax = upper, col = model, fill = model), alpha = 0.25) +
  geom_pointrange(data = s_data_in_temp, aes(x = DPI, y = prevalence, ymin = lower, ymax = upper), col = "grey20") +
  geom_line(aes(x = DPI, y = median, col = model), linewidth = 0.75) +
  facet_wrap(~temp_label) + xlab("Days post infection") + ylab("Sporozoite prevalence") +
  scale_colour_manual(values = c("#56B4E9", "#E69F00"), name = "") +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"), name = "") +
  scale_y_continuous(labels = scales::percent)

### oocysts
# generating the oocyst prevalence estimates from the data
o_data_plot <- generate_prevalence_temp(subset(all_data$oocyst_data, index_g!=1)) %>% 
  mutate(temp_label = paste0(temp, "°C"))

# calculating the mean oocyst intensities from the data for all mosquitoes that were dissected from oocysts
oocyst_number <- generate_oocyst_intensity_summary(subset(all_data$oocyst_data, index_g!=1))
oocyst_number[,"temp"] <- all_data$unique_temp[oocyst_number$index_temp]
oocyst_number[,"temp_label"] <- paste0(oocyst_number[,"temp"],"°C")

# calculating the mean oocyst intensities from the data for infected mosquitoes (with oocysts only)
oocyst_number_inf <- generate_oocyst_intensity_summary(subset(all_data$oocyst_data, index_g!=1 & Oocyst_number>0))
oocyst_number_inf[,"temp"] <- all_data$unique_temp[oocyst_number_inf$index_temp]
oocyst_number_inf[,"temp_label"] <- paste0(oocyst_number_inf[,"temp"],"°C")

# generating the posterior predicted oocyst intensities
O_I_ppd <- rbind(run_prop_ppd_df(fit_temp, "pooled", "O_intsy_ppd", length_ppd_times_O, PPD_times_O),
                 run_prop_ppd_df(fit_all_temp, "independent", "O_intsy_ppd", length_ppd_times_O, PPD_times_O)) %>% 
  mutate(temp_label = paste0(temp, "°C"))

O_I_inf_ppd <- rbind(run_prop_ppd_df(fit_temp, "pooled", "O_intsy_inf_ppd", length_ppd_times_O, PPD_times_O),
                     run_prop_ppd_df(fit_all_temp, "independent", "O_intsy_inf_ppd", length_ppd_times_O, PPD_times_O)) %>% 
  mutate(temp_label = paste0(temp, "°C"))

O_I_plot <- ggplot(data = O_I_ppd) + geom_ribbon(aes(x = DPI, ymin = lower, ymax = upper, col = model, fill = model), alpha = 0.25) +
  geom_pointrange(data = oocyst_number, aes(x = DPI, y = mean, ymin = lower_CI, ymax = upper_CI), col = "grey20") +
  geom_line(aes(x = DPI, y = median, col = model), linewidth = 0.75) +
  facet_wrap(~temp_label) + xlab("Days post infection") + 
  ylab("Mean oocyst intensity\n(among all mosquitoes)") +
  scale_colour_manual(values = c("#56B4E9", "#E69F00"), name = "") +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"), name = "")

O_I_inf_plot <- ggplot(data = O_I_inf_ppd) + geom_ribbon(aes(x = DPI, ymin = lower, ymax = upper, col = model, fill = model), alpha = 0.25) +
  geom_pointrange(data = oocyst_number_inf, aes(x = DPI, y = mean, ymin = lower_CI, ymax = upper_CI), col = "grey20") +
  geom_line(aes(x = DPI, y = median, col = model), linewidth = 0.75) +
  facet_wrap(~temp_label) + xlab("Days post infection") + 
  ylab("Mean oocyst intensity\n(among infected mosquitoes)") +
  scale_colour_manual(values = c("#56B4E9", "#E69F00"), name = "") +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"), name = "")

legend <- get_legend(S_plot)

png("results_temp/fits_plot.png", height = 900, width = 950)
(((S_plot + theme(legend.position = "none")) / (O_I_inf_plot + theme(legend.position = "none"))) | legend) +
  plot_layout(widths = c(0.875, 0.175)) +
  plot_annotation(tag_levels = list(c("a", "b"), "")) &
  theme(plot.tag = element_text(face = 'bold'))
dev.off()

png("results_temp/fits_plot_O.png", height = 850, width = 975)
(((O_I_plot + theme(legend.position = "none")) / (O_I_inf_plot + theme(legend.position = "none"))) | legend) +
  plot_layout(widths = c(0.875, 0.175)) +
  plot_annotation(tag_levels = list(c("a", "b"), "")) &
  theme(plot.tag = element_text(face = 'bold'))
dev.off()

##############################
### actual vs fitted plots ###
##############################

# Brier score calculation for sporozoite infected mosquitoes
# individual mosquito dissection data
s_data_brier <- subset(s_data, gametocytemia != 0.00024)

# matching the individual mosquito dissection data to the predicted sporozoite prevalence values for each day post infection and temperature
s_data_brier$pooled_median <- subset(S_ppd_df, model == "pooled")[match(interaction(round(s_data_brier$DPI, digits = 2), s_data_brier$temp), 
      interaction(round(subset(S_ppd_df, model == "pooled")$DPI, digits = 2), subset(S_ppd_df, model == "pooled")$temp)), "median"]

s_data_brier$ind_median <- subset(S_ppd_df, model == "independent")[match(interaction(round(s_data_brier$DPI, digits = 2), s_data_brier$temp), 
                                                                        interaction(round(subset(S_ppd_df, model == "independent")$DPI, digits = 2), 
                                                                                    subset(S_ppd_df, model == "independent")$temp)), "median"]
# calculating the Brier scores and Brier skill scores
s_data_brier <- s_data_brier %>% rowwise() %>% mutate(res_pooled = (presence - pooled_median)^2,
                        res_ind = (presence - ind_median)^2) %>% ungroup() %>% group_by(temp) %>% 
  mutate(mean = mean(presence),
         ) %>% ungroup() %>% rowwise() %>% mutate(res_mean = (presence - mean)^2)

s_data_brier_text <- s_data_brier %>% group_by(temp) %>% summarise(brier_pooled = sum(res_pooled)/n(),
                                                                   brier_ind = sum(res_ind)/n(),
                                                                   brier_mean = sum(res_mean)/n(),
                                                                   scaled_brier_pooled = 1 - brier_pooled / brier_mean,
                                                                   scaled_brier_ind = 1 - brier_ind / brier_mean) %>% 
  mutate(temp_label = paste0(temp, "°C"))

# plotting
# function to match the actual and predicted values for both models
# data wrangling
gen_a_f_S <- function(s_data_in_temp, 
                      temp_S_ppd, 
                      model){
  inds_m <- match(interaction(round(s_data_in_temp$DPI, digits = 1), s_data_in_temp$temp), interaction(round(temp_S_ppd$DPI, digits = 1), temp_S_ppd$temp))
  a_f_df <- data.frame("median" = temp_S_ppd[inds_m, "median"],
                       "lower" = temp_S_ppd[inds_m, "lower"],
                       "upper" = temp_S_ppd[inds_m, "upper"],
                       "actual" = s_data_in_temp$prevalence,
                       "a_lower" = s_data_in_temp$lower,
                       "a_upper" = s_data_in_temp$upper,
                       "temp" = s_data_in_temp$temp,
                       "model" = rep(model, length(inds_m)))
  return(a_f_df)
}

a_f_df <- rbind(gen_a_f_S(s_data_in_temp, subset(S_ppd_df, model == "pooled"), "pooled"),
                gen_a_f_S(s_data_in_temp, subset(S_ppd_df, model == "independent"), "independent"))

s_data_brier_text <- left_join(s_data_brier_text, a_f_df %>% group_by(temp) %>% summarise(x = max(a_upper) * 0.8,
                                                                                          y = min(upper) + 0.025 * max(upper),
                                                                                          y_ = min(upper) + 0.175 * max(upper)), 
                               by = c("temp"))

actual_fitted_S <- ggplot(data = a_f_df %>% mutate(temp_label = paste0(temp, "°C"))) + 
  geom_pointrange(aes(x = actual, y = median, ymin = lower, ymax = upper, col = model), size = 0.875) + 
  geom_errorbarh(aes(xmax = a_upper, xmin = a_lower, y = median, height = 0, col = model)) + 
  # geom_smooth(formula = y ~ x, aes(x = actual, y = median, col = model), method = "lm", se = FALSE,
  #             size = 1.5) +
  geom_abline(aes(slope=1, intercept=0), linetype = 2, linewidth = 1.5, alpha = 0.75) +
  facet_wrap(~temp_label, scales = "free") +  
  # stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~"), col = model), parse = TRUE,
  #              label.x = 0.95, label.y = 0.05) +
  ylab("Fitted sporozoite prevalence") + xlab("Actual sporozoite prevalence") +
  scale_y_continuous(labels = scales::percent) + scale_x_continuous(labels = scales::percent) +
  scale_colour_manual(values = c("#56B4E9", "#E69F00"), name = "") +
  geom_label(data = s_data_brier_text,
            aes(x = x, y = y, label = paste0("BSS = ", round(scaled_brier_pooled, digits = 2))),
            fill = factor("#E69F00"), size = 5, col = "black", alpha = 0.5) +
  geom_label(data = s_data_brier_text,
             aes(x = x, y = y_, label = paste0("BSS = ", round(scaled_brier_ind, digits = 2))),
             fill = factor("#56B4E9"), size = 5, col = "black", alpha = 0.5)

# oocysts - mean intensity all
gen_a_f_O <- function(oocyst_number, temp_O_I_ppd, model, inf_){
  inds_m_pl <- match(interaction(round(oocyst_number$DPI, digits = 1), oocyst_number$temp), 
                     interaction(round(temp_O_I_ppd$DPI, digits = 1), temp_O_I_ppd$temp))
  
  a_f_df_pl <- data.frame("median" = temp_O_I_ppd[inds_m_pl, "median"],
                          "lower" = temp_O_I_ppd[inds_m_pl, "lower"],
                          "upper" = temp_O_I_ppd[inds_m_pl, "upper"],
                          "actual" = oocyst_number$mean,
                          "a_lower" = oocyst_number$lower_CI,
                          "a_upper" = oocyst_number$upper_CI,
                          "temp" = oocyst_number$temp,
                          "model" = rep(model, length(inds_m_pl)),
                          "mosq_pop" = rep(inf_, length(inds_m_pl)))
  return(a_f_df_pl)
}

a_f_df_pl_plot <- rbind(gen_a_f_O(oocyst_number, subset(O_I_ppd, model == "pooled"), "pooled", "all mosquitoes"),
                        gen_a_f_O(oocyst_number, subset(O_I_ppd, model == "independent"), "independent", "all mosquitoes"),
                        gen_a_f_O(oocyst_number_inf, subset(O_I_inf_ppd, model == "pooled"), "pooled", "infected mosquitoes only"),
                        gen_a_f_O(oocyst_number_inf, subset(O_I_inf_ppd, model == "independent"), "independent", "infected mosquitoes only"))

a_f_df_pl_plot <- a_f_df_pl_plot %>% group_by(model, mosq_pop) %>% 
  mutate(a_mean = mean(actual)) %>% ungroup() %>% rowwise() %>%
  mutate(res = (actual - median)^2,
         res_mean = (actual - a_mean)^2)

r2_o <- a_f_df_pl_plot %>% group_by(model, mosq_pop) %>% summarise(r2 = 1 - (sum(res)/sum(res_mean)),
                                                                   r_lab = paste0("R^'2' == ", round(r2, digits = 2)))

actual_fitted_O <- ggplot(data = a_f_df_pl_plot %>% mutate(temp = paste0(temp,"°C")), aes(x = actual, y = median)) + 
  geom_pointrange(aes(x = actual, y = median, ymin = lower, ymax = upper, col = temp), size = 0.875) +
  geom_errorbarh(aes(xmax = a_upper, xmin = a_lower, height = 0,
                     col = temp)) +
  facet_grid(vars(model), vars(mosq_pop), scales = "free") + 
  geom_abline(aes(slope=1, intercept=0), linetype = 2, linewidth = 1.5, alpha = 0.75) +
  xlab("Actual mean oocyst intensity") + ylab("Fitted mean oocyst intensity") +
  theme(legend.title=element_blank()) +
  geom_label(data = r2_o,
             aes(x = 50, y = 4, 
                 label = r_lab),
            size = 5, col = "black", alpha = 0.5, parse = TRUE)

png("results_temp/actual_fitted_O_S.png", height = 1100, width = 1100)
(actual_fitted_O + plot_spacer() + 
    plot_layout(widths = c(10, 1))) / actual_fitted_S + 
  plot_layout(heights = c(1, 1.5)) +
  plot_annotation(tag_levels = c("a", "b")) &
  theme(plot.tag = element_text(face = 'bold'))
dev.off()

##################
### EIP values ###
##################

# independent model fits
EIP_single <- generate_EIP(fit_all_temp, length(unique(subset(all_data$sporozoite_data, index_g!=1)$temp)))
EIP_df <- EIP_single$EIP_df
EIP_quantile <- EIP_single$EIP_quantiles

seq_temp <- seq(min(EIP_quantile$temp), max(EIP_quantile$temp), 0.01)

degree_day <- rbind(data.frame("temp" = seq_temp,
                               "EIP_" = 111 / (seq_temp - 16),
                               "EIP" = rep("EIP[10]", length(seq_temp))),
                    data.frame("temp" = seq_temp,
                               "EIP_" = 111 / (seq_temp - 16),
                               "EIP" = rep("EIP[50]", length(seq_temp))),
                    data.frame("temp" = seq_temp,
                               "EIP_" = 111 / (seq_temp - 16),
                               "EIP" = rep("EIP[90]", length(seq_temp))))

# for all temperatures
# pooled model estimates
temps <- seq(17, 30, 0.1)
scaled_temps <- (temps - all_data$m_temp) / all_data$sd_temp

params_temp <- rstan::extract(fit_temp)

EIP_index <- get_EIP(params_temp, scaled_temps, 10000)
# index <- seq(1, length(temps))
# mean_EIP <- as.data.frame(t(sapply(index, calc_mean_EIP, EIP_index = EIP_index)))
# mean_EIP$temp <- temps
# colnames(mean_EIP) <- c("mean", "lower", "median", "upper", "n_na", "temp")
# saveRDS(mean_EIP, file = "results_temp/mean_EIP.rds")
mean_EIP <- readRDS(file = "results_temp/mean_EIP.rds")

plot_df <- rbind(cbind(gen_quantiles(EIP_index$EIP_10, temps), data.frame("EIP" = rep("EIP[10]", length(temps)))),
                 cbind(gen_quantiles(EIP_index$EIP_50, temps), data.frame("EIP" = rep("EIP[50]", length(temps)))),
                 cbind(gen_quantiles(EIP_index$EIP_90, temps), data.frame("EIP" = rep("EIP[90]", length(temps)))),
                 cbind(mean_EIP[,-c(5)], data.frame("EIP" = rep("mean", length(temps)))))

write.csv(plot_df, file = "results_temp/temp_model_EIP_values.csv")

# calculating the probability the degree day model EIP is > than EIP50
p_EIP_50_df <- data.frame(temp = temps,
                          p_EIP50 = sapply(seq(1, length(temps)), function(i, EIP_index, temps){
                            sum(EIP_index$EIP_50[,i] < 111/(temps[i] - 16))/length(EIP_index$EIP_50[,i])
                          }, temps = temps, EIP_index = EIP_index))

f <- approxfun(x = subset(p_EIP_50_df, p_EIP50 > 0.02 & p_EIP50 < 0.98)$p_EIP50, 
               y = subset(p_EIP_50_df, p_EIP50 > 0.02 & p_EIP50 < 0.98)$temp, method = "linear")

f(0.95)

p_EIP_50_plot <- ggplot(data = p_EIP_50_df, aes(x = temp, y = p_EIP50)) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(breaks = seq(18, 30, 2)) +
  ylab(expression(paste(italic(p(EIP[50]<EIP[D]))))) + xlab("Temperature (°C)") +
  #geom_hline(yintercept = 0.5, linetype = 2, size = 1) +
  theme(text = element_text(size = 18)) +
  geom_line(linewidth = 1.5) +
  theme(plot.background = element_rect(colour = "grey55", fill = NA, linewidth = 1),
        plot.margin = margin(0.05, 0.25, 0.05, 0.05, "cm"))

EIP_p50_plot <- ggplot() + 
  geom_ribbon(data = subset(plot_df, EIP!="mean"), aes(x = temp, ymin = lower, ymax = upper, fill = factor(EIP)), alpha = 0.45) +
  geom_line(data = subset(plot_df, EIP!="mean"), aes(x = temp, y = median, col = factor(EIP)), linewidth = 1) +
  geom_line(data = degree_day, 
            aes(x = temp, y = EIP_), linetype = 2.5, linewidth = 1.25) +
  geom_pointrange(data = EIP_quantile, aes(x = temp, ymin = lower, y = median, ymax = upper, col = factor(EIP)), size = 0.75) + 
  #geom_line(data = degree_day, aes(x = temp, y = EIP_), linetype = 2, size = 1.25, alpha = 0.75) +
  scale_x_continuous(limits = c(17, 30), breaks = seq(17, 30, 1)) +
  scale_y_continuous(limits = c(0, 115), breaks = seq(0, 110, 10)) +
  ylab("EIP (days)") + xlab("Temperature (°C)") + 
  theme(text = element_text(size = 20)) +
  scale_colour_manual(values = c("#56B4E9", "#E69F00", "#CC79A7"), labels = c(parse(text="EIP[10]"), parse(text="EIP[50]"), parse(text="EIP[90]"))) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00", "#CC79A7"), labels = c(parse(text="EIP[10]"), parse(text="EIP[50]"), parse(text="EIP[90]"))) +
  guides(fill=guide_legend(title=""), colour=guide_legend(title="")) + 
  inset_element(p_EIP_50_plot, left = 0.55, right = 0.99, bottom = 0.5,  top = 0.99) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(face = 'bold'))

png(file = "results_temp/EIP_pEIP50_plot.png", height = 550, width = 900)
EIP_p50_plot
dev.off()

EIP_u <- calc_EIP_v(EIP_index = EIP_index, temps = temps, n_iter = 10000)
#saveRDS(EIP_u, file = "results_temp/EIP_u.rds")
EIP_u <- readRDS(file = "results_temp/EIP_u.rds")

###############################
##### variance in the EIP #####
###############################
var_df <- bind_rows(lapply(seq(1, length(temps)),
                           calc_var_mean_i_EIP,
                           EIP_index = EIP_index)) %>% mutate(temp = temps[index])
 
saveRDS(var_df, file = "results_temp/var_df.rds")
var_df <- readRDS(file = "results_temp/var_df.rds")

var_df <- var_df %>% mutate(cv = sqrt(var)/m) %>% na.omit()

var_df <- var_df %>% group_by(temp) %>% summarise(var_m = median(var),
                                        cv_m = median(cv),
                                        var_lower = quantile(var, probs = 0.025),
                                        var_upper = quantile(var, probs = 0.975),
                                        cv_lower = quantile(cv, probs = 0.025),
                                        cv_upper = quantile(cv, probs = 0.975))

range_10_90_all <- EIP_index$EIP_90 - EIP_index$EIP_10
range_10_90 <- data.frame("median" = apply(range_10_90_all, 2, median),
                          "lower" = apply(range_10_90_all, 2, quantile, probs = c(0.025)),
                          "upper" = apply(range_10_90_all, 2, quantile, probs = c(0.975)),
                          "temp" = temps)

EIP_10_EIP_90_plot <- ggplot(data = range_10_90, 
         aes(x = temp, y = median, ymin = lower, ymax = upper)) +
    geom_ribbon(alpha = 0.25) +
    geom_line(linewidth = 1.5) +
    ylab(expression(paste("Range between ", EIP[10]," and ", EIP[90], " (days)"))) + 
    xlab("Temperature (°C)") +
    scale_x_continuous(limits = c(17, 30), breaks = seq(18, 30, 2)) +
    scale_y_continuous(breaks = seq(0, 30, 5), limits = c(0, 30))

var_plot <- ggplot(data = var_df, aes(x = temp, y = var_m, ymin = var_lower, ymax = var_upper)) +
  geom_ribbon(alpha = 0.25) +
  geom_line(linewidth = 1.5) +
  ylab("Variance of the\nEIP distribution (days)") + xlab("Temperature (°C)")  +
    scale_x_continuous(limits = c(17, 30), breaks = seq(18, 30, 2)) +
    scale_y_continuous(breaks = seq(0, 120, 20), limits = c(0, 120))

CV_plot <- ggplot(data = var_df, aes(x = temp, y = cv_m, ymin = cv_lower, ymax = cv_upper)) +
  geom_ribbon(alpha = 0.25) +
  geom_line(linewidth = 1.5) +
  ylab("Coefficient of variation\nof the EIP distribution") + xlab("Temperature (°C)")  +
    scale_x_continuous(limits = c(17, 30), breaks = seq(18, 30, 2)) +
    scale_y_continuous(breaks = seq(0, 0.25, 0.05), limits = c(0, 0.25))

png(file = "results_temp/var_plot.png", height = 400, width = 1100)
wrap_plots(
  EIP_10_EIP_90_plot,
  
  var_plot,
  
  CV_plot
) + plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(face = 'bold'))
dev.off()

#########################
### vector competence ###
#########################
# pooled model values
iter <- 10000

delta_O <- as.data.frame(matrix(NA, nrow = iter, ncol = length(temps)))
delta_S <- as.data.frame(matrix(NA, nrow = iter, ncol = length(temps)))

for(i in 1:length(scaled_temps)){
  delta_O[,i] <- 1/(1 + exp(-(scaled_temps[i]^2 * params_temp$a_delta + scaled_temps[i] * params_temp$b_delta + params_temp$c_delta)))
  delta_S[,i] <- 1/(1 + exp(-(scaled_temps[i]^2 * params_temp$a_delta_S + scaled_temps[i] * params_temp$b_delta_S + params_temp$c_delta_S)))
}

delta <- gen_quantiles(delta_O * delta_S, temps)

write.csv(delta, file = "results_temp/temp_model_vector_competence_values.csv")

# independent model values

i_delta_O <- as.data.frame(rstan::extract(fit_all_temp, "delta"))
i_delta_S <- as.data.frame(rstan::extract(fit_all_temp, "delta_S"))

i_delta <- i_delta_O * i_delta_S

i_delta_q <- data.frame("temp" = all_data$unique_temp,
                        "lower" = apply(i_delta, 2, quantile, probs = c(0.025)),
                        "median" = apply(i_delta, 2, median),
                        "upper" = apply(i_delta, 2, quantile, probs = c(0.975)),
                        "mean" = apply(i_delta, 2, mean))

######################
### parameter plot ###
######################
get_param_values <- function(fit, param){
  placeholder <- as.data.frame(rstan::extract(fit, param)[[1]])[,1]
  data.frame("mean" = mean(placeholder), "median" = median(placeholder),
             "lower" = quantile(placeholder, c(0.05))[[1]], "upper" = quantile(placeholder, c(0.95))[[1]])
}

get_all_params <- function(fit, fit_all, param, all_data){
  n_t <- length(all_data$unique_temp)
  all_params <- as.data.frame(matrix(NA, nrow = n_t * 2, ncol = 5))
  
  for(i in 1:n_t){
    all_params[i,] <- cbind(get_param_values(fit, paste0(param,"[",i,"]")),
                            "temp" = all_data$unique_temp[i])
    
    all_params[i+n_t,] <- cbind(get_param_values(fit_all, paste0(param,"[",i,"]")),
                                "temp" = all_data$unique_temp[i])
  }
  
  all_params <- cbind(all_params, data.frame("model" = c(rep("pooled", n_t), rep("independent", n_t)),
                                             "param" = rep(param, nrow(all_params))))
  colnames(all_params) <- c("mean", "median", "lower", "upper", "temp", "model", "param")
  return(all_params)
}

out <- rbind(get_all_params(fit_temp, fit_all_temp, "delta", all_data),
             get_all_params(fit_temp, fit_all_temp, "delta_S", all_data),
             get_all_params(fit_temp, fit_all_temp, "shape_S", all_data),
             get_all_params(fit_temp, fit_all_temp, "rate_S", all_data),
             get_all_params(fit_temp, fit_all_temp, "mu", all_data),
             get_all_params(fit_temp, fit_all_temp, "mu_total_S", all_data))

out$param <- as.character(out$param)

out[which(out$param == "delta"), "param"] <- "Parameter: delta[O]"
out[which(out$param == "delta_S"), "param"] <- "Parameter: delta[S]"
out[which(out$param == "shape_S"), "param"] <- "Parameter: alpha[OS]"
out[which(out$param == "rate_S"), "param"] <- "Parameter: beta[OS]"
out[which(out$param == "mu"), "param"] <- "Parameter: mu"
out[which(out$param == "mu_total_S"), "param"] <- "Parameter: alpha[GS]/beta[GS]"

# smoothing the parameter values 
# for all temperatures
temps <- seq(17, 30, 0.1)
scaled_temps <- (temps - all_data$m_temp) / all_data$sd_temp

# shape parameter
smooth_S <- cbind(scaled_temps, temps, as.data.frame(matrix(NA, nrow = length(temps), ncol = 6)))
colnames(smooth_S) <- c("temp", "temp_", "mean", "median", "lower", "upper", "param", "model")

for(i in 1:nrow(smooth_S)){
  placeholder <- rstan::extract(fit_temp, "a_shape_S")[[1]] * smooth_S[i,"temp"]^2 + 
    rstan::extract(fit_temp, "b_shape_S")[[1]] * smooth_S[i,"temp"] + 
    rstan::extract(fit_temp, "c_shape_S")[[1]]
  smooth_S[i, "mean"] <- mean(placeholder)
  smooth_S[i, "median"] <- median(placeholder)
  smooth_S[i, "lower"] <- quantile(placeholder, c(0.05))[[1]]
  smooth_S[i, "upper"] <- quantile(placeholder, c(0.95))[[1]]
  smooth_S[i, "param"] <- "Parameter: alpha[OS]"
  smooth_S[i, "model"] <- "pooled"
}

# rate parameter
smooth_R <- cbind(scaled_temps, temps, as.data.frame(matrix(NA, nrow = length(temps), ncol = 6)))
colnames(smooth_R) <- c("temp", "temp_", "mean", "median", "lower", "upper", "param", "model")

for(i in 1:nrow(smooth_R)){
  placeholder <- rstan::extract(fit_temp, "m_rate_S")[[1]] * smooth_R[i,"temp"] + 
    rstan::extract(fit_temp, "c_rate_S")[[1]]
  smooth_R[i, "mean"] <- mean(placeholder)
  smooth_R[i, "median"] <- median(placeholder)
  smooth_R[i, "lower"] <- quantile(placeholder, c(0.05))[[1]]
  smooth_R[i, "upper"] <- quantile(placeholder, c(0.95))[[1]]
  smooth_R[i, "param"] <- "Parameter: beta[OS]"
  smooth_R[i, "model"] <- "pooled"
}

# delta parameter
smooth_D <- cbind(scaled_temps, temps, as.data.frame(matrix(NA, nrow = length(temps), ncol = 6)))
colnames(smooth_D) <- c("temp", "temp_", "mean", "median", "lower", "upper", "param", "model")

for(i in 1:nrow(smooth_D)){
  placeholder <- 1/(1+exp(-(rstan::extract(fit_temp, "a_delta")[[1]] * smooth_D[i,"temp"]^2 + 
                              rstan::extract(fit_temp, "b_delta")[[1]] * smooth_D[i,"temp"] + 
                              rstan::extract(fit_temp, "c_delta")[[1]])))
  smooth_D[i, "mean"] <- mean(placeholder)
  smooth_D[i, "median"] <- median(placeholder)
  smooth_D[i, "lower"] <- quantile(placeholder, c(0.05))[[1]]
  smooth_D[i, "upper"] <- quantile(placeholder, c(0.95))[[1]]
  smooth_D[i, "param"] <- "Parameter: delta[O]"
  smooth_D[i, "model"] <- "pooled"
}


smooth_M <- cbind(scaled_temps, temps, as.data.frame(matrix(NA, nrow = length(temps), ncol = 6)))
colnames(smooth_M) <- c("temp", "temp_", "mean", "median", "lower", "upper", "param", "model")

for(i in 1:nrow(smooth_M)){
  placeholder <- (rstan::extract(fit_temp, "a_shape_S")[[1]] * smooth_M[i,"temp"]^2 + 
                    rstan::extract(fit_temp, "b_shape_S")[[1]] * smooth_M[i,"temp"] + 
                    rstan::extract(fit_temp, "c_shape_S")[[1]]) /
    (rstan::extract(fit_temp, "m_rate_S")[[1]] * smooth_M[i,"temp"] + 
       rstan::extract(fit_temp, "c_rate_S")[[1]])
  smooth_M[i, "mean"] <- mean(placeholder)
  smooth_M[i, "median"] <- median(placeholder)
  smooth_M[i, "lower"] <- quantile(placeholder, c(0.05))[[1]]
  smooth_M[i, "upper"] <- quantile(placeholder, c(0.95))[[1]]
  smooth_M[i, "param"] <- "Parameter: alpha[GS]/beta[GS]"
  smooth_M[i, "model"] <- "pooled"
}

smooth_mu <- cbind(scaled_temps, temps, as.data.frame(matrix(NA, nrow = length(temps), ncol = 6)))
colnames(smooth_mu) <- c("temp", "temp_", "mean", "median", "lower", "upper", "param", "model")

for(i in 1:nrow(smooth_mu)){
  placeholder <- (1/(1+exp(-(rstan::extract(fit_temp, "a_mu")[[1]] * smooth_mu[i,"temp"]^2 + 
                               rstan::extract(fit_temp, "b_mu")[[1]] * smooth_mu[i,"temp"] + 
                               rstan::extract(fit_temp, "c_mu")[[1]])))) * rstan::extract(fit_temp, "g_mu")[[1]]
  smooth_mu[i, "mean"] <- mean(placeholder)
  smooth_mu[i, "median"] <- median(placeholder)
  smooth_mu[i, "lower"] <- quantile(placeholder, c(0.05))[[1]]
  smooth_mu[i, "upper"] <- quantile(placeholder, c(0.95))[[1]]
  smooth_mu[i, "param"] <- "Parameter: mu"
  smooth_mu[i, "model"] <- "pooled"
}

smooth_DS <- cbind(scaled_temps, temps, as.data.frame(matrix(NA, nrow = length(temps), ncol = 6)))
colnames(smooth_DS) <- c("temp", "temp_", "mean", "median", "lower", "upper", "param", "model")

for(i in 1:nrow(smooth_DS)){
  placeholder <- 1/(1+exp(-(rstan::extract(fit_temp, "a_delta_S")[[1]] * smooth_DS[i,"temp"]^2 + 
                              rstan::extract(fit_temp, "b_delta_S")[[1]] * smooth_DS[i,"temp"] + 
                              rstan::extract(fit_temp, "c_delta_S")[[1]])))
  smooth_DS[i, "mean"] <- mean(placeholder)
  smooth_DS[i, "median"] <- median(placeholder)
  smooth_DS[i, "lower"] <- quantile(placeholder, c(0.05))[[1]]
  smooth_DS[i, "upper"] <- quantile(placeholder, c(0.95))[[1]]
  smooth_DS[i, "param"] <- "Parameter: delta[S]"
  smooth_DS[i, "model"] <- "pooled"
}

smooth_params <- rbind(smooth_S, smooth_R, smooth_D, smooth_M, smooth_mu, smooth_DS)
smooth_params$label <- smooth_params$param

out$temp_ <-out$temp
out$label <- out$param

colnames(delta)[1] <- "temp_"
param_plot_df <- rbind(subset(smooth_params, param!= "Parameter: beta[OS]" & param!= "Parameter: alpha[OS]" &
                                param!= "Parameter: alpha[GS]/beta[GS]")[-1],
                       cbind(delta, data.frame("param" = rep("VC", nrow(delta)),
                                               "model" = rep("pooled", nrow(delta)),
                                               "label" = rep("delta[O]*delta[S]", nrow(delta)))))
# changing the name of mu to gamma
out[which(out$label == "Parameter: mu"), "label"] <- "Parameter: mu"
param_plot_df[which(param_plot_df$label == "Parameter: mu"), "label"] <- "Parameter: mu"

# plotting each individually
# mean oocyst load

smooth_mu <- smooth_mu %>% mutate(param = "Parameter: mu")

mu_plot <- ggplot(data = smooth_mu) +
  geom_pointrange(data = oocyst_number_inf, aes(x = temp, y = mean, ymin = lower_CI, ymax = upper_CI), alpha = 0.75, shape = 21,
                  fill = "black", col = "black", size = 0.875) +
  geom_ribbon(aes(x = temp_, ymin = lower, ymax = upper, fill = model), alpha = 0.2) +
  geom_line(aes(x = temp_, y = median, col = model), linewidth = 1.5) +
  geom_pointrange(data = subset(out, label == "Parameter: mu" & model == "independent"),  
                  aes(x = temp_, y = median, ymin = lower, ymax = upper, col = model, fill = model),
                  alpha = 0.75, shape = 21, size = 0.875) +
  scale_x_continuous(limits = c(17, 30), breaks = seq(17, 30, 1)) +
  scale_colour_manual(values = c("#56B4E9", "#E69F00"), name = "Model") + 
  scale_fill_manual(values = c("#56B4E9", "#E69F00"), name = "Model") +
  facet_wrap(~label, scales = "free", labeller =  label_parsed) + 
  scale_y_continuous(breaks = seq(0, 80, 10)) +
  ylab("Mean oocyst load\namong infected mosquitoes") + xlab("Temperature (°C)")

delta_plot <- ggplot(data = smooth_D) +
  geom_ribbon(aes(x = temp_, ymin = lower, ymax = upper, fill = model), alpha = 0.2) +
  geom_line(aes(x = temp_, y = median, col = model), linewidth = 1.5) +
  geom_pointrange(data = o_data_plot, aes(x = temp, y = prevalence, ymin = lower, ymax = upper), alpha = 0.75, shape = 21,
                  fill = "black", col = "black", size = 0.875) +
  geom_pointrange(data = subset(out, label == "Parameter: delta[O]" & model == "independent"),  
                  aes(x = temp_, y = median, ymin = lower, ymax = upper, col = model, fill = model),
                  alpha = 0.75, shape = 21, size = 0.875) +
  scale_x_continuous(limits = c(17, 30), breaks = seq(17, 30, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = scales::percent) +
  scale_colour_manual(values = c("#56B4E9", "#E69F00"), name = "Model") + 
  scale_fill_manual(values = c("#56B4E9", "#E69F00"), name = "Model") +
  facet_wrap(~label, scales = "free", labeller =  label_parsed) + 
  ylab("Human-to-mosquito transmission\nprobability (to oocyst life stage)") + 
  xlab("Temperature (°C)")

delta_S_plot <- ggplot(data = smooth_DS) +
  geom_ribbon(aes(x = temp_, ymin = lower, ymax = upper, fill = model), alpha = 0.2) +
  geom_line(aes(x = temp_, y = median, col = model), linewidth = 1.5) +
  geom_pointrange(data = subset(out, label == "Parameter: delta[S]" & model == "independent"),  
                  aes(x = temp_, y = median, ymin = lower, ymax = upper, col = model, fill = model),
                  alpha = 0.75, shape = 21, size = 0.875) +
  scale_x_continuous(limits = c(17, 30), breaks = seq(17, 30, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = scales::percent) +
  scale_colour_manual(values = c("#56B4E9", "#E69F00"), name = "Model") + 
  scale_fill_manual(values = c("#56B4E9", "#E69F00"), name = "Model") +
  facet_wrap(~label, scales = "free", labeller =  label_parsed) + 
  ylab("Conversion of oocyst to\nsalivary gland sporozoite infection") + 
  xlab("Temperature (°C)")

HMTP_plot <- ggplot(data = delta %>% mutate(model = "pooled",
                               label = "delta[O]*delta[S]")) +
  geom_ribbon(aes(x = temp_, ymin = lower, ymax = upper, fill = model), alpha = 0.2) +
  geom_line(aes(x = temp_, y = median, col = model), linewidth = 1.5) +
  geom_pointrange(data = i_delta_q %>% mutate(model = "independent"), 
                  aes(x = temp, y = median, ymin = lower, ymax = upper, fill = model, col = model), 
                  shape = 21, alpha = 0.75, size = 0.875) +
  scale_colour_manual(values = c("#56B4E9", "#E69F00"), name = "Model") + 
  scale_fill_manual(values = c("#56B4E9", "#E69F00"), name = "Model") +
  facet_wrap(~label, scales = "free", labeller =  label_parsed) +
  ylab("Human-to-mosquito transmission\nprobability (to sporozoite life stage)") + xlab("Temperature (°C)") +
  scale_x_continuous(limits = c(17, 30), breaks = seq(17, 30, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = scales::percent)

png(file = "results_temp/parameter_plot_with_data.png", width = 1250, height = 750)
mu_plot + 
    delta_plot + 
    delta_S_plot + 
    HMTP_plot + 
    plot_layout(ncol = 2) +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(face = 'bold'))
dev.off()



 