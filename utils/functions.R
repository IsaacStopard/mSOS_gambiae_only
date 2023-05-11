#################
### functions ###
#################
EIP_PDF <- function(t, a, b, mu, k){
  return(-(((1/b)^-a) * exp(-b * t) * (t^(-1+a)) * mu * ((k / (k + mu))^(1 + k)) * (((k + mu * Rgamma(a, b*t)) / (k + mu)) ^ (-1 - k))
           / ((-1 + (k / (k + mu))^k) * gamma(a))))
}

EIP_CDF <- function(a, b, mu, k, t){
  return(
    (1 - k^k * (k + mu * Rgamma(a, b * t, lower = TRUE))^-k) / 
      (1 - (k / (k + mu))^k)
  )
}

Inv_EIP_CDF <- function(a, b, mu, k, p){
  one <- ((-(p * (1 - (k / (k + mu))^k) - 1)/(k^k))^(-1/k) - k) / mu
  return((Rgamma.inv(a, one, lower = TRUE)/b))
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

# functions for running the multiple model
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


### generate oocyst intensity data quantiles
generate_oocyst_intensity_summary <- function(data){
  out <- unique(data[,c("DPI", "index_gt", "index_temp", "index_g")])
  out$mean <- rep(NA, nrow(out))
  out$median <- rep(NA, nrow(out))
  out$lower <- rep(NA, nrow(out))
  out$upper <- rep(NA, nrow(out))
  
  for(i in 1:nrow(out)){
    placeholder <- subset(data, DPI == out[i,"DPI"] & index_gt == out[i, "index_gt"])
    out[i, "mean"] <- mean(placeholder[,"Oocyst_number"])
    out[i, "median"] <- quantile(placeholder[,"Oocyst_number"], 0.5)
    out[i, "lower"] <- quantile(placeholder[,"Oocyst_number"], 0.025)
    out[i, "upper"] <- quantile(placeholder[,"Oocyst_number"], 0.975)
  }
  return(out)
}

### generate the mean oocyst intensity among infected mosquitoes only - mean over all DPI
gen_O_intsy_sum_inf <- function(data){
  data <- subset(data, Oocyst_number!=0)
  out <- unique(data[,c("index_gt", "index_temp", "index_g")])
  out$mean <- rep(NA, nrow(out))
  out$median <- rep(NA, nrow(out))
  out$lower <- rep(NA, nrow(out))
  out$upper <- rep(NA, nrow(out))
  
  for(i in 1:nrow(out)){
    placeholder <- subset(data, index_gt == out[i, "index_gt"])
    out[i, "mean"] <- mean(placeholder[,"Oocyst_number"])
    out[i, "median"] <- quantile(placeholder[,"Oocyst_number"], 0.5)
    out[i, "lower"] <- quantile(placeholder[,"Oocyst_number"], 0.025)
    out[i, "upper"] <- quantile(placeholder[,"Oocyst_number"], 0.975)
  }
  return(out)
}
