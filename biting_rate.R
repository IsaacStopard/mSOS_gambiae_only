library(readxl); library(tidyverse); library(rstan);

df <- readxl::read_excel("data/doi_10.5061_dryad.74839_shapiro/PLoS.Biology.DataFigures.Supplemental.xlsx",
                   sheet = "Fig4.Raw.Data")

df_long <- df %>% tidyr::pivot_longer(cols = c("clutch.1", "clutch.2"), names_to = "clutch", values_to =  "gono_length")

ggplot(data = df_long,
       aes(x = temp, y = gono_length, col = id, group = factor(id))) +
  geom_point()

ggplot(data = df %>% tidyr::pivot_longer(cols = c("clutch.1", "clutch.2"), names_to = "clutch", values_to =  "gono_length"),
       aes(x = clutch, y = gono_length, col = id, group = factor(id))) +
  geom_point() +
  geom_line() +
  facet_wrap(~temp)

ggplot(data = df %>% tidyr::pivot_longer(cols = c("clutch.1", "clutch.2"), 
                                         names_to = "clutch", values_to =  "gono_length"),
       aes(x = id, y = gono_length, group = factor(clutch))) +
  geom_point() +
  facet_wrap(~temp) +
  geom_smooth(method = "lm")

df %>% group_by(temp) %>% summarise(sum(is.na(clutch.2)))

Briere_model <- "
data{
  int<lower=0> n; // number of data points
  real temp[n]; // 
  int gono[n];
  int n_pred_temp;
  vector[n_pred_temp] pred_temp;
}

parameters{
  real a; // Briere function parameter
  real T0; // Briere function parameter
  real Tm; // Briere function parameter
}

transformed parameters{
  vector[n] Briere;
  vector[n] m; // mean biting rate
  vector<lower = 0, upper = 1>[n] p; // geometric distribution parameter
  
  for(i in 1:n){
    Briere[i] = temp[i] >= T0 && temp[i] <= Tm ? a/1000 * temp[i] * (temp[i] - T0) * sqrt((Tm - temp[i])): 0;
    m[i] = 1/Briere[i];
    p[i] = 1/(m[i] + 1);
  }
  
}

model{
// priors from Mordecai paper
a ~ normal(0.000203*1000, 1);
T0 ~ normal(11.7, 4.5);
Tm ~ normal(42.3, 4.5);

// likelihood
for(i in 1:n){
  target += neg_binomial_lpmf(gono[i] | 1, p[i]/(1-p[i]));
}
}

generated quantities{
vector[n_pred_temp] Briere_out;
  for(i in 1:n_pred_temp){
    Briere_out[i] = pred_temp[i] >= T0 && pred_temp[i] <= Tm ? a/1000 * pred_temp[i] * (pred_temp[i] - T0) * sqrt((Tm - pred_temp[i])) : 0;
  }
  
}

"
stanDso_in <- stan_model(model_code = Briere_model) 

df_s <- na.omit(df_long)

mean_temp <- mean(df_s$temp)
sd_temp <- sd(df_s$temp)

df_s$scaled_temp <- (df_s$temp - mean_temp)/sd_temp

pred_temp <- seq(0, 60, 0.1)

data_in <- list(n = nrow(df_s),
                temp = df_s$temp,
                gono = as.integer(df_s$gono_length),
                n_pred_temp = length(pred_temp),
                pred_temp = pred_temp)

finit <- function(){
  list(a = 0.00001,
       T0 = 15,
       Tm = 50
  )
}

fit <- sampling(stanDso_in, 
                data = data_in, 
                iter = 2500, 
                warmup= 1250,
                init = finit,
                cores = 4,
                control = list(max_treedepth = 14,
                               stepsize = 0.1)) 

saveRDS(fit, file = "fits/biting_rate_fit.rds")

pred_vals <- rstan::extract(fit, "Briere_out")$Briere_out

mean_biting_rate <- data.frame(temp = pred_temp, 
           med = apply(pred_vals, 2, median, na.rm = TRUE),
           lower = apply(pred_vals, 2, quantile, prob = 0.025, na.rm = TRUE),
           upper = apply(pred_vals, 2, quantile, prob = 0.975, na.rm = TRUE))

write.csv(mean_biting_rate, file = "biting_rate_estimates.csv")

ggplot(data = mean_biting_rate) +
  geom_ribbon(aes(x = temp, ymin = lower, ymax = upper), alpha = 0.25) +
  geom_line(aes(x = temp, y = med)) +
  theme_bw() +
  geom_point(data = na.omit(df_long) %>% group_by(temp) %>% summarise(mean = mean(1/gono_length),
                                                                           lower = quantile(1/gono_length, prob = 0.025),
                                                                           upper = quantile(1/gono_length, prob = 0.975)), 
             aes(x = temp, y = mean), col = "skyblue", size = 3) +
  ylab("Biting rate")
  
ggplot(data = na.omit(df_long),
       aes(x = temp, y = 1/gono_length, group = factor(temp))) +
  geom_boxplot()