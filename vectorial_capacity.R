# estimating the temp-dependent vectorial capacity
# degree-day model and mSOS model comparison
rm(list = ls())
library(tidyverse); library(rstan); library(shinystan); 
library(cowplot); library(zipfR); library(truncnorm); library(ggpmisc); library(patchwork);
library(doParallel); library(foreach)

# Assuming mosquito survival is exponentially distributed
# survival was assumed to follow Marten's models
source(file = "utils/functions_temp_only.R")
source(file = "utils/plotting_functions_temp_only.R")

prop_surv <- function(t, mu, shape, rate, mu_PL, k){
  return(exp(-mu * t) * EIP_PDF(t, a = shape, b = rate, mu = mu_PL, k = k))
}

prop_surv_int <- function(mu, shape, rate, mu_PL, k){
  out <- tryCatch(
    {stats::integrate(prop_surv, lower = 0, upper = Inf, mu = mu, shape = shape, rate = rate, mu_PL = mu_PL, k = k)[[1]]},
    error=function(cond){return(NA)}
    )
  return(out)
}

v.prop_surv_int <- Vectorize(prop_surv_int)

# survival function models
martens_1 <- function(T){
  return(-0.0016*T^2 + 0.054*T+0.45)
}

martens_2 <- function(T){
  return(exp(-(1/(-4.4 + 1.31*T - 0.03*T^2))))
}

theme_set(theme_bw() + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  text = element_text(size = 18)))

# EIP model
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

temps <- seq(17, 30, 0.1)
scaled_temps <- (temps - all_data$m_temp) / all_data$sd_temp

fit_temp <- readRDS("fits/fit_mSOS_temp_only_f2_f3.rds")
# for all temperatures
# pooled model estimates
params_temp <- rstan::extract(fit_temp)

EIP_index <- get_EIP(params_temp, scaled_temps, 10000)

# calculating the proportion surviving the EIP
params_df <- data.frame(temp = temps) %>% mutate(M1 = (1 - martens_1(temp)),
                                            M2 = (1 - martens_2(temp)),
                                            degree_day = 111 / (temp - 16)) %>% 
  mutate(p_s_dd = exp(- M2 * degree_day))

ggplot(data = params_df) +
  geom_line(aes(x = temp, y = M1)) +
  geom_line(aes(x = temp, y = M2), col = "skyblue")

# calculating the probability of surviving the EIP
cl <- makeCluster(3)
registerDoParallel(cl)
clusterExport(cl = cl, varlist = c("params_df", "EIP_index", "prop_surv", "EIP_PDF", "v.prop_surv_int", "prop_surv_int"))
all_p_s <- foreach(i=1:nrow(params_df),
                   .packages = (.packages())) %dopar% 
  {
    v.prop_surv_int(mu = params_df[i, "M2"], 
                    shape = EIP_index$shape_total_S[,i],
                    rate = EIP_index$rate_total_S[,i],
                    mu_PL = EIP_index$mu[,i],
                    k = EIP_index$k[,i])
   
}
stopCluster(cl)

model_M2 <- bind_rows(lapply(seq(1, length(all_p_s)), function(i, all_p_s){
  if(sum(is.na(all_p_s[[i]]))<10){
    return(quantile(all_p_s[[i]], probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
  } else{
    return(NA)
  }
}, all_p_s = all_p_s)) %>% as.data.frame()

plot_grid(
  ggplot(data = cbind(params_df, model_M2)) + 
  geom_line(aes(x = temp, y = p_s_dd), linewidth = 1) +
  geom_line(aes(x = temp, y = `50%`), linewidth = 1, col = "skyblue") +
  ylab("Probability of surviving the EIP"),
  
  ggplot(data = cbind(params_df, model_M2)) + 
  geom_line(aes(x = temp, y = `50%`/p_s_dd), linewidth = 1) +
  ylab("Ratio of probability of surviving the mSOS EIP\nestimate compared to degree day model estimate") +
    scale_y_log10(),
  
  labels = c("A", "B")
)













