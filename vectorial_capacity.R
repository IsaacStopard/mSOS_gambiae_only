# estimating the temp-dependent vectorial capacity
# degree-day model and mSOS model comparison
rm(list = ls())
library(tidyverse); library(rstan); library(shinystan); 
library(cowplot); library(zipfR); library(truncnorm); library(ggpmisc); library(patchwork);
library(doParallel); library(foreach)
library(calculus)

theme_set(theme_bw() + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  text = element_text(size = 18)))

#####################
##### functions #####
#####################

source(file = "utils/functions_temp_only.R")
source(file = "utils/plotting_functions_temp_only.R")

# mosquito mortality functions
gomp_par_g1 <- function(temp, b0, b2, b3){
  return(b0 + b2*temp + b3*temp^2)
}

gomp_par_g2 <- function(temp, y, s){
  return(y + s*temp)
}

gomp_fun <- function(a, g1, g2){
  return(g1 * exp(a * g2))
}

surv_gomp <- function(a0, a, g1, g2){
  return(exp(-g1 / g2 * exp(a0 * g2) * (exp(a * g2) - 1)))
}

surv_gomp_i <- function(a0, a, g1, g2){
  int <- stats::integrate(gomp_fun, lower = a0, upper = a0 + a, g1 = g1, g2 = g2)[[1]]
  return(exp(-int))
}

surv_gomp(0, 10, 0.01, 0.1)
surv_gomp_i(0, 10, 0.01, 0.1)

# Step 1
# probability of surviving the EIP given infection at age a0 and the Gompertz survival function
# integrand
prop_surv_gomp_int <- function(a, # time since age at infection
                           a0, # age at infection
                           g1, g2, # gompertz function parameters
                           shape, rate, mu_PL, k # EIP PDF parameters)
                           ){
  return(EIP_PDF(t = a, a = shape, b = rate, mu = mu_PL, k = k) * # EIP PDF
                        surv_gomp(a0 = a0, a = a, g1 = g1, g2 = g2))
}

# integral
prop_surv_gomp <- function(a0, # age at infection
                               g1, 
                               g2, # gompertz function parameters
                               shape, 
                               rate, 
                               mu_PL, 
                               k # EIP PDF parameters
                                ){
  out <- calculus::integral(prop_surv_gomp_int, # mosquito survival
                            bounds = list("a" = c(0, Inf)),
                            params = list("a0" = a0, 
                                          "g1" = g1, 
                                          "g2" = g2,
                                          "shape" = shape, 
                                          "rate" = rate, 
                                          "mu_PL" = mu_PL, 
                                          "k" = k),
                            method = "cuhre",
                            coordinates = "cartesian"
                            )
  return(out$value)
}

prop_surv_gomp_fixed <- function(a, # degree-day model EIP
                                 a0, # age at infection
                                 g1, 
                                 g2){
  return(surv_gomp(a0 = a0, a = a, g1 = g1, g2 = g2))
}

# Step 2
# probability a mosquito biting time
# probability of j bites given mosquito exits the EIP at age a1
# integrand
p_j_bites_EIP_a1_int <- function(a2, # age the mosquito dies
                                 j, # number of bites
                                 alpha, # Poisson distribution rate
                                 a1, # age the mosquito exits the EIP
                                 g1,
                                 g2){
  
    # out <- ((alpha ^ j * (a2 - a1) ^ j) / factorial(j)) *
    #     exp(a2 * g2 - (a2 - a1) * alpha) *
    #     g1 * exp(- g1 / g2 * (exp(a2 * g2) - exp(a1 * g2)))
    
    # out <- (alpha ^ j * (a2 - a1) ^ j) / factorial(j) * g1 * exp(-(-a1 + a2)*alpha + a2 * g2 - (g1 * (-exp(a1 * g2) + exp(a2 * g2)))/g2)
    
    # log scale for stability
    p1 <- j * log(alpha  * (a2 - a1)) - log(factorial(j))
    p2 <- a2 * g2 - (a2 - a1) * alpha
    p3 <- log(g1)

    fp1 <- a2 * g2
    fp2 <- a1 * g2

    p4 <- - g1 / g2 * (exp(fp1) - exp(fp2))

    out <- exp(p1 + p2 + p3 + p4)
    
    return(out)
}

# integral
p_j_bites_EIP_a1 <- function(j, # number of bites
                          alpha, # Poisson distribution rate
                          a1, # age the mosquito exits the EIP
                          g1,
                          g2
                      ){
  # out <- stats::integrate(p_j_bites_EIP_a1_int, 
  #                         lower = a1, 
  #                         upper = Inf, 
  #                         j = j, 
  #                         alpha = alpha, 
  #                         a1 = a1, 
  #                         g1 = g1, 
  #                         g2 = g2)[[1]]
  
  out <- calculus::integral(p_j_bites_EIP_a1_int, 
                            bounds = list("a2" = c(a1, Inf)),
                            params = list("j" = j,
                                          "alpha" = alpha,
                                          "a1" = a1,
                                          "g1" = g1,
                                          "g2" = g2),
                            method = "cuhre",
                            absTol = 1E-7,
                            coordinates = "cartesian",
                            )
  
  return(out$value)
}

v.p_j_bites_EIP_a1 <- Vectorize(p_j_bites_EIP_a1)

v.p_j_bites_EIP_a1(j = seq(0, 100), alpha = 0.5, a1 = 5, g1 = 0.01, g2 = 0.1)

# expected number of bites
E_j_bites_EIP_a1_int <- function(j, # number of bites
                             alpha, # Poisson distribution rate
                             a1, # age the mosquito exits the EIP
                             g1,
                             g2){
  
  return(j * v.p_j_bites_EIP_a1(j = j, alpha = alpha, a1 = a1, g1 = g1, g2 = g2))
  
}

E_j_bites_EIP_a1 <- function(alpha, # Poisson distribution rate
                             a1, # age the mosquito exits the EIP
                             g1,
                             g2){
  
  out <- calculus::integral(
    E_j_bites_EIP_a1_int,
    bounds = list("j" = c(0, Inf)),
    params = list("alpha" = alpha, "a1" = a1, "g1" = g1, "g2" = g2),
    method = "cuhre",
    absTol = 1E-7,
    coordinate = "cartesian"
    )$value
  
  # this code is the same as in Github
  # out <- stats::integrate(function(a2, # age the mosquito dies
  #                                  alpha, # Poisson distribution rate
  #                                  a1, # age the mosquito exits the EIP
  #                                  g1,
  #                                  g2){
  #   
  #   out <- -((a1 - a2) * alpha * g1 * exp((exp(a1 * g2) * g1 - exp(a2 * g2) * g1 + a2 * g2^2)/g2))
  #   
  #   return(out)
  # }, lower = a1, upper = Inf, alpha = alpha,a1 = a1, g1 = g1, g2 = g2)[[1]]
  
  return(out)
}

# sum approximation - won't use this for now
# E_j_bites_EIP_a1 <- function(alpha, # Poisson distribution rate
#                       a1, # age the mosquito exits the EIP
#                       g1,
#                       g2){
#   # find the maximum number of bites
#   j <- 0
#   p <- p_j_bites_EIP_a1(j = j, alpha = alpha, a1 = a1, g1 = g1, g2 = g2)
# 
#   j_vec <- j
#   p_vec <- p
# 
#   while (TRUE) {
#     j <- j + 1
# 
#     p <- p_j_bites_EIP_a1(j = j, alpha = alpha, a1 = a1, g1 = g1, g2 = g2)
# 
#     j_vec <- c(j_vec, j)
#     p_vec <- c(p_vec, p)
# 
#     if (sum(p_vec) >= 1 & p > 1E-10) break
#   }
# 
#   return(list("E" = sum(j_vec * p_vec), "df" = data.frame(j = j_vec, p = p_vec)))
# }

E_j_bites_EIP_a1(alpha = 0.5, a1 = 4, g1 = 0.01, g2 = 0.1)

# Step 3
# expected number of bites given infectious blood meal at age a0
p_j_bites_meal_a0_int <- function(a1, # age the mosquito exits the EIP
                              j,
                              a0, # age the mosquito is infected
                              alpha, # Poisson distribution rate
                              g1,
                              g2,
                              shape, 
                              rate, 
                              mu_PL, 
                              k # EIP PDF parameters
    ){
  out <- p_j_bites_EIP_a1(j = j, 
                            alpha = alpha, 
                            a1 = a1, 
                            g1 = g1, 
                            g2 = g2) *
    EIP_PDF(a1 - a0, 
            a = shape, 
            b = rate, 
            mu = mu_PL, 
            k = k) *
    surv_gomp(a0 = a0, 
              a = a1 - a0, 
              g1 = g1, 
              g2 = g2)
  return(out)
    
}

v.p_j_bites_meal_a0_int <- Vectorize(p_j_bites_meal_a0_int)

v.p_j_bites_meal_a0_int(j = 1,
                        a1 = seq(10, 100),
                        a0 = 10, 
                      alpha = 0.5, 
                      g1 = 0.01, 
                      g2 = 0.1,
                      shape = 80, 
                      rate = 4, 
                      mu_PL = 20, 
                      k = 0.5)

v.p_j_bites_meal_a0_int(j = 0,
                        a1 = seq(24, 400, 1),
                        a0 = params_df[i, "a0"],
                        alpha = params_df[i, "alpha"],
                        g1 = params_df[i, "g1"],
                        g2 = params_df[i, "g2"],
                        shape = params_df[i, "shape_total_S"],
                        rate = params_df[i, "rate_total_S"],
                        mu_PL = params_df[i, "mu"],
                        k = params_df[i, "k"])

p_j_bites_meal_a0 <- function(j,
                              a0, # age the mosquito is infected
                              alpha, # Poisson distribution rate
                              g1,
                              g2,
                              shape, 
                              rate, 
                              mu_PL, 
                              k){
  
  out <- calculus::integral(v.p_j_bites_meal_a0_int, # mosquito survival
                            bounds = list("a1" = c(a0, Inf)),
                            params = list("j" = j,
                                          "a0" = a0,
                                          "alpha" = alpha,
                                          "g1" = g1, 
                                          "g2" = g2,
                                          "shape" = shape, 
                                          "rate" = rate, 
                                          "mu_PL" = mu_PL, 
                                          "k" = k),
                            method = "cuhre",
                            absTol = 1E-7, # 1E-7
                            coordinates = "cartesian")$value
  
  if(j == 0){
    out <- out + (1 - prop_surv_gomp(a0 = a0, g1 = g1, g2 = g2, shape = shape, rate = rate, mu_PL = mu_PL, k = k))
  }
  
  return(out)
}

v.p_j_bites_meal_a0 <- Vectorize(p_j_bites_meal_a0)

v.p_j_bites_meal_a0(j = seq(0, 10),
               a0 = params_df[i, "a0"],
             alpha = params_df[i, "alpha"],
             g1 = params_df[i, "g1"],
             g2 = params_df[i, "g2"],
             shape = params_df[i, "shape_total_S"],
             rate = params_df[i, "rate_total_S"],
             mu_PL = params_df[i, "mu"],
             k = params_df[i, "k"])

# degree-day model equivalent

p_j_bites_meal_a0_fixed <- function(j,
                                    a0, # age the mosquito is infected
                                    alpha, # Poisson distribution rate
                                    g1,
                                    g2,
                                    EIP# EIP PDF parameters
){
  a1 <- a0 + EIP # age the mosquito exits the EIP
  out <- p_j_bites_EIP_a1(j = j, 
                          alpha = alpha, 
                          a1 = a1, 
                          g1 = g1, 
                          g2 = g2) *
    surv_gomp(a0 = a0, 
              a = a1 - a0, 
              g1 = g1, 
              g2 = g2)
  
  if(j == 0){
    out <- out + (1 - prop_surv_gomp_fixed(a0 = a0, a = EIP, g1 = g1, g2 = g2))
  }
  
  return(out)
  
}

v.p_j_bites_meal_a0_fixed <- Vectorize(p_j_bites_meal_a0_fixed)

v.p_j_bites_meal_a0_fixed(j = seq(1, 10),
                        a0 = 1, 
                        alpha = 0.5, 
                        g1 = 0.01, 
                        g2 = 0.1,
                        EIP = 10)

# expected number of bites
E_j_bites_a0_int <- function(j,
                             a0, # age the mosquito exits the EIP
                             alpha, # Poisson distribution rate
                             g1,
                             g2,
                             shape,
                             rate,
                             mu_PL,
                             k){
  return(j * p_j_bites_meal_a0(j = j, 
                                 a0 = a0, # age the mosquito is infected
                                 alpha = alpha, # Poisson distribution rate
                                 g1 = g1,
                                 g2 = g2,
                                 shape = shape, 
                                 rate = rate, 
                                 mu_PL = mu_PL, 
                                 k = k))
  
}

E_j_bites_a0 <- function(a0, # age the mosquito exits the EIP
                         alpha, # Poisson distribution rate
                         g1,
                         g2,
                         shape,
                         rate,
                         mu_PL,
                         k){
  
  out <- calculus::integral(E_j_bites_a0_int,
                            bounds = list("j" = c(0, Inf)),
                            params = list("a0" = a0, # age the mosquito is infected
                                          "alpha" = alpha, # Poisson distribution rate
                                          "g1" = g1,
                                          "g2" = g2,
                                          "shape" = shape, 
                                          "rate" = rate, 
                                          "mu_PL" = mu_PL, 
                                          "k" = k),
                            method = "cuhre",
                            absTol = 1E-7, # 1E-7
                            coordinates = "cartesian"
                            )$value
  return(out)
}

v.E_j_bites_a0 <- Vectorize(E_j_bites_a0)

E_j_bites_a0_int_fixed <- function(j,
                             a0, # age the mosquito exits the EIP
                             alpha, # Poisson distribution rate
                             g1,
                             g2,
                             EIP){
  
  return(j * p_j_bites_meal_a0_fixed(j = j, 
                               a0 = a0, # age the mosquito is infected
                               alpha = alpha, # Poisson distribution rate
                               g1 = g1,
                               g2 = g2,
                               EIP = EIP))
  
}

v.E_j_bites_a0_int_fixed <- Vectorize(E_j_bites_a0_int_fixed)

E_j_bites_a0_fixed <- function(a0, # age the mosquito exits the EIP
                         alpha, # Poisson distribution rate
                         g1,
                         g2,
                         EIP){
  
  out <- calculus::integral(v.E_j_bites_a0_int_fixed,
                            bounds = list("j" = c(0, Inf)),
                            params = list("a0" = a0, # age the mosquito is infected
                                          "alpha" = alpha, # Poisson distribution rate
                                          "g1" = g1,
                                          "g2" = g2,
                                          "EIP" = EIP),
                            method = "cuhre",
                            absTol = 1E-7,
                            coordinates = "cartesian"
  )$value
  return(out)
}

v.E_j_bites_a0_fixed <- Vectorize(E_j_bites_a0_fixed)

# step 4
E_j_bites_int <- function(a0, 
                          alpha, 
                          g1, g2,
                          shape, rate, mu_PL, k){
    return(E_j_bites_a0(a0 = a0, 
                        alpha = alpha, 
                        g1 = g1, 
                        g2 = g2,
                        shape = shape, 
                        rate = rate, 
                        mu_PL = mu_PL, 
                        k = k) * 
             alpha * exp(- alpha * a0))
    }

v.E_j_bites_int <- Vectorize(E_j_bites_int)

E_j_bites <- function(alpha, # Poisson distribution rate
                      g1,
                      g2,
                      shape,
                      rate,
                      mu_PL,
                      k){
  
  out <- calculus::integral(v.E_j_bites_int,
    bounds = list("a0" = c(0, Inf)),
    params = list("alpha" = alpha, # Poisson distribution rate
                  "g1" = g1, "g2" = g2,
                  "shape" = shape, "rate" = rate, "mu_PL" = mu_PL, "k" = k),
    method = "suave",
    absTol = 1E-5,
    coordinates = "cartesian"
  )$value

  return(out)
  }

v.E_j_bites <- Vectorize(E_j_bites)

# E_j_bites_a0_check <- function(a0, # age the mosquito exits the EIP
#                          alpha, # Poisson distribution rate
#                          g1,
#                          g2,
#                          shape,
#                          rate,
#                          mu_PL,
#                          k){
#   # find the maximum number of bites
#   j <- 0
#   p <- p_j_bites_meal_a0(j = j, 
#                          alpha = alpha, 
#                          a0 = a0, 
#                          g1 = g1, 
#                          g2 = g2, 
#                          shape = shape, 
#                          rate = rate, 
#                          mu_PL = mu_PL, 
#                          k = k)
#   
#   j_vec <- j
#   p_vec <- p
#   
#   while (TRUE) {
#     
#     j <- j + 1
#     p <- p_j_bites_meal_a0(j = j, 
#                            alpha = alpha, 
#                            a0 = a0, 
#                            g1 = g1, 
#                            g2 = g2, 
#                            shape = shape, 
#                            rate = rate, 
#                            mu_PL = mu_PL, 
#                            k = k)
#     
#     j_vec <- c(j_vec, j)
#     p_vec <- c(p_vec, p)
#     
#     if (j > 100) break #sum(p_vec) >= 1 & p < 10E-5
#   }
#   
#   return(list("E" = sum(j_vec * p_vec), "p_vec" = p_vec))
# }
# 
# E_j_bites_a0_check(a0 = 10, 
#              alpha = 0.5, 
#              g1 = 0.01, 
#              g2 = 0.1,
#              shape = 40, 
#              rate = 4, 
#              mu_PL = 10, 
#              k = 0.5)

quadratic_function <- function(t, a, b, c){
  out <- (t^2 * a) + (t * b) + c
  return(out)
}

logistic_function <- function(c){
  out <- 1 / (1 + exp(-c))
  return(out)
}

degree_day <- function(temp){
  return(111 / (temp - 16))
}

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

fit_temp <- readRDS("fits/fit_mSOS_temp_only_f2_f3.rds")

params_temp <- rstan::extract(fit_temp)

mean_biting_rate <- read.csv(file = "biting_rate_estimates.csv")

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


params_df[,"alpha"] <- mean_biting_rate[match(params_df$temp, mean_biting_rate$temp), "med"]

# params_df %>% rowwise() %>% mutate(max_a0 = round(uniroot(
#   function(x){surv_gomp(a0 = 0, a = x, g1 = g1, g2 = g2) - 0.25},  lower = 0, upper = 250)[[1]], digits = 0))

extract_params <- function(scaled_temp, params_temp){
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
  
  return(data.frame(shape_total_S = median(shape_total_S),
                    rate_total_S = median(rate_total_S),
                    mu = median(mu),
                    k = median(k)))
}

params_df <- cbind(params_df, lapply(params_df$scaled_temp, extract_params, params_temp = params_temp) %>% bind_rows %>% as.data.frame())

n <- nrow(params_df)
a0 <- seq(0, 20, 1)

params_df <- params_df %>% slice(rep(1:n(), each = length(a0))) %>% mutate(a0 =rep(a0, n))

params_df <- params_df %>% mutate(E_j_a0 = v.E_j_bites_a0(a0 = a0,
                                                          alpha = alpha,
                                                          g1 = g1,
                                                          g2 = g2,
                                                          shape = shape_total_S,
                                                          rate = rate_total_S,
                                                          mu_PL = mu,
                                                          k = k),
                                  EIP_DD = degree_day(temp),
                                  E_j_a0_fixed = v.E_j_bites_a0_fixed(a0 = a0,
                                                                      alpha = alpha,
                                                                      g1 = g1,
                                                                      g2 = g2,
                                                                      EIP = EIP_DD))

saveRDS(params_df, file = "params_df_E_j_bites_a0.rds")

params_df <- readRDS(file = "params_df_E_j_bites_a0.rds")

params_df$location <- factor(params_df$location, levels = c("Kericho", "Kitale", "Kisumu", "Garissa"))

relative_E_j_a0 <- rbind(params_df %>% dplyr::select("location", "temp_source", "feed", "a0", "E_j_a0") %>% tidyr::pivot_wider(names_from = temp_source,
                                                                                                values_from = E_j_a0) %>% mutate(model = "EIP_dist"),
                         params_df %>% dplyr::select("location", "temp_source", "feed", "a0", "E_j_a0_fixed") %>% tidyr::pivot_wider(names_from = temp_source,
                                                                                                                       values_from = E_j_a0_fixed) %>% 
                           mutate(model = "degree_day")) %>% mutate(ratio = future / recent)

params_df <- params_df %>% 
  mutate(temp_source = ifelse(temp_source == "future", "Future mean temperature", "Recent mean temperature"))

params_df$temp_source <- factor(params_df$temp_source, levels = c("Recent mean temperature", "Future mean temperature"))

vc_plot <- ggplot(data = params_df %>% pivot_longer(cols = c("E_j_a0", "E_j_a0_fixed"), 
                                         names_to = "model", values_to = "VC") ,
       aes(x = a0, y = VC, group = interaction(feed, model, temp_source, location))) + 
  geom_point(aes(fill = model, shape = feed), size = 3, alpha = 0.5) +
  geom_line(aes(colour = model), linewidth = 0.75) + 
  facet_grid(vars(location), vars(temp_source)) +
  ylab("Expected number of infectious bites\nfrom per infected mosquito (z)") + 
  xlab("Mosquito age when infected (days)") +
  scale_shape_manual(values = c(21, 22), labels = c("Feed1", "Feed2"), name = "Feed") +
  scale_colour_manual(values = c("#009E73", "grey40"), name = "EIP model", labels = c("Suh-Stopard", "Degree-day")) +
  scale_fill_manual(values = c("#009E73", "grey40"), name = "EIP model", labels = c("Suh-Stopard", "Degree-day"))

leg <- cowplot::get_legend(vc_plot)

r_plot <- ggplot(data = relative_E_j_a0, 
       aes(x = a0, y = ratio, group = interaction(model, feed, location))) +
  geom_point(aes(fill = model, shape = feed), size = 3, alpha = 0.5) +
  geom_line(aes(colour = model), linewidth = 0.75) +
  facet_wrap(vars(location), scales = "free_y") +
  scale_y_sqrt() +
  ylab("Ratio of future:recent expected infectious\nbites per infected mosquito (rVC)") +
  xlab("Mosquito age when infected (days)") +
  scale_shape_manual(values = c(21, 22), labels = c("Feed1", "Feed2"), name = "Feed") +
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

cbind(subset(params_df, a0 == 0 & temp_source == "recent" & feed == "feed_1")$E_j_a0 %>% round(digits = 2),
           subset(params_df, a0 == 0 & temp_source == "recent" & feed == "feed_1")$E_j_a0_fixed %>% round(digits = 2),
           subset(params_df, a0 == 0 & temp_source == "recent" & feed == "feed_2")$E_j_a0 %>% round(digits = 2),
           subset(params_df, a0 == 0 & temp_source == "recent" & feed == "feed_2")$E_j_a0_fixed %>% round(digits = 2)) %>% view()

cbind(subset(relative_E_j_a0, feed == "feed_1" & model == "EIP_dist" & a0 == 0)$ratio %>% round(digits = 2), 
      subset(relative_E_j_a0, feed == "feed_1" & model == "degree_day" & a0 == 0)$ratio %>% round(digits = 2),
      subset(relative_E_j_a0, feed == "feed_2" & model == "EIP_dist" & a0 == 0)$ratio %>% round(digits = 2), 
      subset(relative_E_j_a0, feed == "feed_2" & model == "degree_day" & a0 == 0)$ratio %>% round(digits = 2)) %>% view()


# Assuming mosquito survival is exponentially distributed
# survival was assumed to follow Marten's models

prop_surv <- function(t, mu, shape, rate, mu_PL, k){
  return(exp(-mu * t) * EIP_PDF(t, a = shape, b = rate, mu = mu_PL, k = k))
}

prop_surv_int <- function(mu, shape, rate, mu_PL, k){
  out <- tryCatch(
    {stats::integrate(prop_surv, lower = 0, upper = Inf, mu = mu, shape = shape, rate = rate, mu_PL = mu_PL, k = k,
                      rel.tol = 1e-15)[[1]]},
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













