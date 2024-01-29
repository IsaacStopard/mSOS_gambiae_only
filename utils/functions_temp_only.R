#################
### functions ###
#################

# function that stratifies by temperature and feed and calculates the prevalence
generate_prevalence <- function(data){
  totals <- unique(data[,c("DPI","index_temp", "temp", "index_g", "index_gt", "gametocytemia")])
  
  for(i in 1:nrow(totals)){
    totals[i,"sample"] <- length(which(data[,"DPI"] == totals[i, "DPI"]
                                       & data[,"index_temp"] == totals[i, "index_temp"] &
                                         data[,"index_g"] == totals[i, "index_g"] &
                                         data[,"index_gt"] == totals[i, "index_gt"] &
                                         data[,"gametocytemia"] == totals[i, "gametocytemia"]))
    
    totals[i,"positive"] <- length(which(data[,"DPI"] == totals[i, "DPI"] 
                                         & data[,"presence"] > 0 & data[,"index_temp"] == totals[i, "index_temp"] &
                                           data[,"index_g"] == totals[i, "index_g"] &                                             
                                           data[,"index_gt"] == totals[i, "index_gt"] &
                                           data[,"gametocytemia"] == totals[i, "gametocytemia"]))
  }
  
  rm(i)
  
  totals <- mutate(totals, prevalence = positive / sample) # prevalence
  totals <- mutate(totals, lower = prevalence - (1.96 * sqrt(prevalence * (1 - prevalence) / sample))) # 5% CI
  totals <- mutate(totals, upper = prevalence + (1.96 * sqrt(prevalence * (1 - prevalence) / sample))) # 95% CI
  totals[which(totals$lower < 0), "lower"] <- 0 # preventing the lower confidence interval being below 0
  
  return(totals)
}

# function that stratifies by temperature and calculates the prevalence
generate_prevalence_temp <- function(data){
  totals <- unique(data[,c("DPI","index_temp", "temp")])
  
  for(i in 1:nrow(totals)){
    totals[i,"sample"] <- length(which(data[,"DPI"] == totals[i, "DPI"]
                                       & data[,"index_temp"] == totals[i, "index_temp"]))
    
    totals[i,"positive"] <- length(which(data[,"DPI"] == totals[i, "DPI"] 
                                         & data[,"presence"] > 0 & data[,"index_temp"] == totals[i, "index_temp"]))
  }
  
  rm(i)
  
  conf <- DescTools::BinomCI(totals$positive, totals$sample, conf.level = 0.95, sides = "two.sided") %>% as.data.frame()
  
  totals <- mutate(totals, prevalence = positive / sample,
                   lower = conf[,"lwr.ci"],
                   upper = conf[,"upr.ci"]
                   ) # prevalence
  
  return(totals)
}

# function that tables the number of unique oocyst counts
oocyst_intensity_indexing_temp <- function(data){
  U_data <- unique(data[,c("DPI", "Oocyst_number", "index_temp")])
  index <- rep(NA, nrow(data))
  for(i in 1:nrow(U_data)){
    index[which(data$DPI == U_data[i, "DPI"] &
                  data$Oocyst_number == U_data[i, "Oocyst_number"] &
                  data$index_temp == U_data[i, "index_temp"])] <- i
  }
  out <- list("unique_oocyst_intensity" = U_data, "oocyst_intensity_index" = index)
  rm(list = c("U_data", "index"))
  return(out)
}

temp_g_indexing <- function(M_data, S_data){
  
  # temperature indexing
  temps <- append(M_data$temp, S_data$temp)
  U_temp <- sort(unique(temps))
  U_temp_scaled <- sort(unique(scale(temps, center = TRUE, scale = TRUE)[,1])) # converting temperature onto the same scale as used in the other model
  
  sd_temp <- sd(temps)
  m_temp <- mean(temps)
  
  for(i in 1:length(U_temp)){
    M_data[which(M_data$temp == U_temp[i]),"temp_scaled"] <- U_temp_scaled[i]
    S_data[which(S_data$temp == U_temp[i]),"temp_scaled"] <- U_temp_scaled[i]
    M_data[which(M_data$temp == U_temp[i]),"index_temp"] <- i
    S_data[which(S_data$temp == U_temp[i]),"index_temp"] <- i
  }
  
  # g indexing - not for survival
  # specifying the g index
  g_ind <- unique(rbind(M_data[,c("gametocytemia", "index_ref")], S_data[,c("gametocytemia", "index_ref")]))
  g_ind$index <- seq(1, nrow(g_ind), 1) 
  
  for(i in 1:nrow(g_ind)){
    M_data[which(M_data$index_ref == g_ind[i, "index_ref"] & M_data$gametocytemia == g_ind[i,"gametocytemia"]), "index_g"] <- g_ind[i, "index"]
    S_data[which(S_data$index_ref == g_ind[i, "index_ref"] & S_data$gametocytemia == g_ind[i,"gametocytemia"]), "index_g"] <- g_ind[i, "index"]
  }
  
  unique_gt <- unique(rbind(M_data[,c("index_temp", "index_g")], S_data[,c("index_temp", "index_g")]))
  unique_gt$index <- seq(1, nrow(unique_gt), 1)
  
  for(i in 1:nrow(unique_gt)){
    M_data[which(M_data$index_temp == unique_gt[i,"index_temp"] & 
                   M_data$index_g == unique_gt[i, "index_g"]), "index_gt"] <- unique_gt[i, "index"]
    
    S_data[which(S_data$index_temp == unique_gt[i,"index_temp"] & 
                   S_data$index_g == unique_gt[i, "index_g"]), "index_gt"] <- unique_gt[i, "index"]
  }
  
  
  out <- list("unique_temp" = U_temp, "unique_temp_scaled" = U_temp_scaled, "unique_g" = g_ind, "unique_gt" = unique_gt,
              "oocyst_data" = M_data, "sporozoite_data" = S_data, "m_temp" = m_temp, "sd_temp" = sd_temp)
  
  rm(list = c("temps", "U_temp", "U_temp_scaled", "M_data", "S_data", "m_temp", "sd_temp"))
  
  return(out)
}

oocyst_intensity_indexing <- function(data){
  U_data <- unique(data[,c("DPI", "Oocyst_number", "index_gt")])
  index <- rep(NA, nrow(data))
  for(i in 1:nrow(U_data)){
    index[which(data$DPI == U_data[i, "DPI"] &
                  data$Oocyst_number == U_data[i, "Oocyst_number"] &
                  data$index_gt == U_data[i, "index_gt"])] <- i
  }
  out <- list("unique_oocyst_intensity" = U_data, "oocyst_intensity_index" = index)
  rm(list = c("U_data", "index"))
  return(out)
}


# function to calculate the mean EIP from the get_EIP function
calc_mean_EIP <- function(EIP_index, i){
  p <- data.frame("a" = EIP_index$shape_total_S[,i],
                  "b" = EIP_index$rate_total_S[,i],
                  "mu" = EIP_index$mu[,i],
                  "k" = EIP_index$k[,i]) %>% rowwise() %>% mutate(mean = tryCatch({integrate(mean_func, a = a, b = b,
                                                                                             mu = mu, k = k, lower = 0, upper = Inf)[[1]]},
                                                                                  error = function(e){
                                                                                    NA
                                                                                  }))
  return(c("mean" = mean(na.omit(p$mean)), quantile(na.omit(p$mean), c(0.025, 0.5, 0.975)),
           "n_na" = sum(is.na(p$mean))))
}

calc_var_mean_i_EIP <- function(EIP_index, i){
  m <- data.frame("a" = EIP_index$shape_total_S[,i],
                  "b" = EIP_index$rate_total_S[,i],
                  "mu" = EIP_index$mu[,i],
                  "k" = EIP_index$k[,i]) %>% rowwise() %>% mutate(mean = tryCatch({calculus::integral(mean_func, 
                                                                                                      bounds = list("t" = c(0, Inf)),
                                                                                                      params = list("a" = a,
                                                                                                                    "b" = b,
                                                                                                                    "mu" = mu,
                                                                                                                    "k" = k),
                                                                                                      method = "cuhre",
                                                                                                      coordinates = "cartesian")[[1]]},
                                                                                  error = function(e){
                                                                                    NA
                                                                                  }),
                                                                  var = tryCatch({calculus::integral(int_var_func, 
                                                                                           bounds = list("t" = c(0, Inf)),
                                                                                           params = list("a" = a,
                                                                                                         "b" = b,
                                                                                                         "mu" = mu,
                                                                                                         "k" = k),
                                                                                           method = "cuhre",
                                                                                           coordinates = "cartesian")[[1]]},
                                                                                 error = function(e){
                                                                                   NA
                                                                                 }))
  out <- data.frame("m" = m$mean, "v" = m$var) %>% mutate("var" = v - m^2,
                                                 "index" = i)
  return(out)
}

# new method for extracting the EIP
generate_EIP <- function(fit_all, l_n = 11){
  fit_ <- rstan::extract(fit_all)
  EIP_df <- data.frame("a" = fit_$shape_total_S[,1], "b" = fit_$rate_total_S[,1], "mu" = fit_$mu[,1], "k" = fit_$k[,1], "index_temp" = rep(1, (iterations - warmup)*chains))
  EIP_df[,"EIP[10]"] <- Inv_EIP_CDF(EIP_df$a, EIP_df$b, EIP_df$mu, EIP_df$k, 0.1)
  EIP_df[,"EIP[50]"] <- Inv_EIP_CDF(EIP_df$a, EIP_df$b, EIP_df$mu, EIP_df$k, 0.5)
  EIP_df[,"EIP[90]"] <- Inv_EIP_CDF(EIP_df$a, EIP_df$b, EIP_df$mu, EIP_df$k, 0.9)
  
  EIP_quantiles <-cbind(data.frame("index_temp" = rep(1, 3), "EIP" = c("EIP[10]", "EIP[50]", "EIP[90]"), "mean" = c(mean(EIP_df[,"EIP[10]"]),
                                                                                                                    mean(EIP_df[,"EIP[50]"]),
                                                                                                                    mean(EIP_df[,"EIP[90]"]))),
                        as.data.frame(rbind(quantile(EIP_df[,"EIP[10]"], c(0.025, 0.5, 0.975)),
                                            quantile(EIP_df[,"EIP[50]"], c(0.025, 0.5, 0.975)),
                                            quantile(EIP_df[,"EIP[90]"], c(0.025, 0.5, 0.975)))))
  colnames(EIP_quantiles) <- c("index_temp", "EIP", "mean", "lower", "median", "upper")
  
  for(i in 2:l_n){
    placeholder_df <- data.frame("a" = fit_$shape_total_S[,i], "b" = fit_$rate_total_S[,i], "mu" = fit_$mu[,i], "k" = fit_$k[,i], "index_temp" = rep(i, (iterations - warmup)*chains))
    placeholder_df[,"EIP[10]"] <- Inv_EIP_CDF(placeholder_df$a, placeholder_df$b, placeholder_df$mu, placeholder_df$k, 0.1)
    placeholder_df[,"EIP[50]"] <- Inv_EIP_CDF(placeholder_df$a, placeholder_df$b, placeholder_df$mu, placeholder_df$k, 0.5)
    placeholder_df[,"EIP[90]"] <- Inv_EIP_CDF(placeholder_df$a, placeholder_df$b, placeholder_df$mu, placeholder_df$k, 0.9)
    EIP_df <- rbind(EIP_df, placeholder_df)
    
    placeholder_quantiles <-cbind(data.frame("index_temp" = rep(i, 3), "EIP" = c("EIP[10]", "EIP[50]", "EIP[90]"), "mean" = c(mean(placeholder_df[,"EIP[10]"]),
                                                                                                                              mean(placeholder_df[,"EIP[50]"]),
                                                                                                                              mean(placeholder_df[,"EIP[90]"]))),
                                  as.data.frame(rbind(quantile(placeholder_df[,"EIP[10]"], c(0.025, 0.5, 0.975)),
                                                      quantile(placeholder_df[,"EIP[50]"], c(0.025, 0.5, 0.975)),
                                                      quantile(placeholder_df[,"EIP[90]"], c(0.025, 0.5, 0.975)))))
    colnames(placeholder_quantiles) <- c("index_temp", "EIP", "mean", "lower", "median", "upper")
    EIP_quantiles <- rbind(EIP_quantiles, placeholder_quantiles)
    rm(list = c("placeholder_df", "placeholder_quantiles"))
  }
  
  EIP_quantiles$temp <- all_data$unique_temp[EIP_quantiles$index_temp]
  return(list("EIP_df" = EIP_df, "EIP_quantiles" = EIP_quantiles))
}

EIP_PDF <- function(t, a, b, mu, k){
  return(-(((1/b)^-a) * exp(-b * t) * (t^(-1+a)) * mu * ((k / (k + mu))^(1 + k)) * (((k + mu * zipfR::Rgamma(a, b*t)) / (k + mu)) ^ (-1 - k))
           / ((-1 + (k / (k + mu))^k) * gamma(a))))
}

EIP_PDF_pca <- function(t, a, b, mu, k){
  return(-(((1/b)^-a) * exp(-b * t) * (t^(-1+a)) * mu * ((k / (k + mu))^(1 + k)) * (((k + mu * (1 -incgam(b*t, a)/gamma(a))) / (k + mu)) ^ (-1 - k))
           / ((-1 + (k / (k + mu))^k) * gamma(a))))
}

EIP_PDF_2 <- function(t, a, b, mu, k){
  return((b*exp(-b*t)*((b*t)^(-1+a))*k^(1+k)*mu*((k+mu*Rgamma(a, b*t))^(-1-k)))/((1 - (k / (k + mu))^k) * gamma(a)))
}

EIP_CDF <- function(a, b, mu, k, t){
  return(
    (1 - k^k * (k + mu * Rgamma(a, b * t, lower = TRUE))^-k) / 
      (1 - (k / (k + mu))^k)
  )
}

# function to calculate the mean of the EIP PDF 
# to get the mean this function must be integrated
mean_func <- function(t, a, b, mu, k){
  return(t * EIP_PDF(t, a, b, mu, k))
}

# function to calculate the variance of the EIP PDF 
# to get the variance this function must be integrated
int_var_func <- function(t, a, b, mu, k){
  return(t^2 * EIP_PDF(t, a, b, mu, k))
}

Inv_EIP_CDF <- function(a, b, mu, k, p){
  one <- ((-(p * (1 - (k / (k + mu))^k) - 1)/(k^k))^(-1/k) - k) / mu
  return((zipfR::Rgamma.inv(a, one, lower = TRUE)/b))
}

prop_ppd_function <- function(fit_df, n_unq_gt, length_ppd_times, PPD_times, iterations, warmup, chains, Stan_data_name){
  prop_ppd <- array(NaN, c(length_ppd_times, ((iterations - warmup) * chains), n_unq_gt))
  for(i in 1:length_ppd_times){
    for(j in 1:n_unq_gt){
      prop_ppd[i,,j] <- fit_df[,paste0(Stan_data_name,"[",i,",",j,"]")]
    }
  }
  prop_ppd_df <- data.frame()
  labs_gt_ind <- c()
  for(i in 1:n_unq_gt){
    place <- as.data.frame(prop_ppd[,,i])
    prop_ppd_df <- rbind(prop_ppd_df, place)
    labs_gt_ind <-  append(labs_gt_ind, rep(i, length_ppd_times))
  }
  
  day_post_inf <- rep(PPD_times, n_unq_gt)
  prop_ppd_df[,"DPI"] <- day_post_inf
  prop_ppd_df[,"index_gt"] <- labs_gt_ind
  prop_ppd_df$median <- apply(prop_ppd_df[,1:((iterations - warmup)*chains)], 1, median)
  prop_ppd_df$lower <- apply(prop_ppd_df[,1:((iterations - warmup)*chains)], 1, quantile, probs = c(0.025))
  prop_ppd_df$upper <- apply(prop_ppd_df[,1:((iterations - warmup)*chains)], 1, quantile, probs = c(0.975))
  prop_ppd_df$mean <- apply(prop_ppd_df[,1:((iterations - warmup)*chains)], 1, mean)
  
  x <- ((iterations - warmup) * chains) 
  prop_ppd_df <- prop_ppd_df%>%gather("iteration", "value", 1:x)
  prop_ppd_df$percentile <- rep("iteration", nrow(prop_ppd_df))
  
  prop_quantile_ppd_df <- subset(prop_ppd_df, iteration == "V1")
  prop_quantile_ppd_df <- prop_quantile_ppd_df[,c("DPI", "index_gt", "median", "lower", "upper", "mean")]
  
  return(prop_quantile_ppd_df)
}

generate_oocyst_intensity_summary <- function(data){
  out <- unique(data[,c("DPI", "index_temp")])
  out$mean <- rep(NA, nrow(out))
  out$median <- rep(NA, nrow(out))
  out$lower <- rep(NA, nrow(out))
  out$upper <- rep(NA, nrow(out))
  
  for(i in 1:nrow(out)){
    placeholder <- subset(data, DPI == out[i,"DPI"] & index_temp == out[i, "index_temp"])
    out[i, "mean"] <- mean(placeholder[,"Oocyst_number"])
    out[i, "median"] <- quantile(placeholder[,"Oocyst_number"], 0.5)
    out[i, "lower"] <- quantile(placeholder[,"Oocyst_number"], 0.025)
    out[i, "upper"] <- quantile(placeholder[,"Oocyst_number"], 0.975)
    l.model <- lm(Oocyst_number ~ 1, data = placeholder)
    out[i, "lower_CI"] <- confint(l.model, level=0.95)[[1]]
    out[i, "upper_CI"] <- confint(l.model, level=0.95)[[2]]
  }
  return(out)
}

run_prop_ppd_df <- function(fit, model, Stan_data_name, length_ppd_times, PPD_times){
  df <- prop_ppd_function(as.data.frame(fit),length(all_data$unique_temp), length_ppd_times, PPD_times, iterations, warmup, chains, Stan_data_name)
  df$model <- rep(model, nrow(df))
  df$temp <- all_data$unique_temp[df$index_gt]
  return(df)
}

get_EIP <- function(params, temps, iter){
  shape_S <- as.data.frame(matrix(NA, nrow = iter, ncol = length(temps)))
  rate_S <- as.data.frame(matrix(NA, nrow = iter, ncol = length(temps)))
  mu <- as.data.frame(matrix(NA, nrow = iter, ncol = length(temps)))
  
  shape_O <- as.data.frame(replicate(length(temps), params$shape_O))
  rate_O <- as.data.frame(replicate(length(temps), params$rate_O))
  
  k <- as.data.frame(replicate(length(temps), params$k))
  
  for(i in 1:length(temps)){
    shape_S[,i] <- temps[i]^2 * params$a_shape_S + temps[i] * params$b_shape_S + params$c_shape_S
    rate_S[,i] <- temps[i] * params$m_rate_S + params$c_rate_S
    mu[,i] <- (1 / (1 + exp(-(temps[i]^2 * params$a_mu + temps[i] * params$b_mu + params$c_mu)))) * params$g_mu
    #mu[,i] <- params$a_mu * exp(-((temps[i] - params$b_mu)^2 / (2 * params$c_mu)^2))
  }
  
  mu_total_S <- (rate_O * shape_S + rate_S * shape_O) / (rate_O * rate_S)
  sigma_sq_S  <- (rate_O^2 * shape_S + rate_S^2 * shape_O) / (rate_O^2 * rate_S^2)
  shape_total_S <- mu_total_S^2 / sigma_sq_S;
  rate_total_S <- mu_total_S / sigma_sq_S;
  mean_total_S <- shape_total_S / rate_total_S
  
  EIP_10 <- as.data.frame(matrix(NA, nrow = iter, ncol = length(temps)))
  EIP_50 <- as.data.frame(matrix(NA, nrow = iter, ncol = length(temps)))
  EIP_90 <- as.data.frame(matrix(NA, nrow = iter, ncol = length(temps)))
  
  for(i in 1:length(temps)){
    EIP_10[,i] <- Inv_EIP_CDF(shape_total_S[,i], rate_total_S[,i], mu[,i], k[,i], 0.1)
    EIP_50[,i] <- Inv_EIP_CDF(shape_total_S[,i], rate_total_S[,i], mu[,i], k[,i], 0.5)
    EIP_90[,i] <- Inv_EIP_CDF(shape_total_S[,i], rate_total_S[,i], mu[,i], k[,i], 0.9)
  }
  
  return(list("shape_O" = shape_O, "rate_O" = rate_O, "shape_S" = shape_S, "rate_S" = rate_S, "shape_total_S" = shape_total_S, 
              "rate_total_S" = rate_total_S, "mean_total_S" = mean_total_S, "EIP_10" = EIP_10, "EIP_50" = EIP_50, "EIP_90" = EIP_90,
              "mu" = mu, "k" = k))
}

gen_quantiles <- function(df, temps){
  out <- as.data.frame(matrix(NA, nrow = length(temps), ncol = 5))
  
  for(i in 1:length(temps)){
    out[i,] <- c(temps[i], quantile(df[,i], c(0.025, 0.5, 0.975)), mean(df[,i]))
  }
  colnames(out) <- c("temp", "lower", "median", "upper", "mean")
  return(out)
}
