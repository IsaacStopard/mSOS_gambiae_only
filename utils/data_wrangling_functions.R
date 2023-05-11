########################
### process raw data ###
########################
process_prevalence_data <- function(data){
  out <- data.frame("Experiment" = integer(),
                    "Cup" = integer(),
                    "DPI" = double(),
                    "presence" = integer(),
                    "temp" = double(),
                    "DTR" = double(),
                    "gametocytemia" = double())
  
  for(i in 1:nrow(data)){
    n_inf <- data[i, "Infected"]
    n_un <- data[i, "Dissected"] - data[i, "Infected"]
    if(n_inf > 0){
      place_inf <- data.frame("Experiment" = rep(data[i, "Experiment"], n_inf),
                              "Cup" = rep(data[i, "Cup.ID"], n_inf),
                              "DPI" = rep(data[i, "DPI"], n_inf),
                              "presence" = rep(1, n_inf),
                              "temp" = rep(data[i, "temp"], n_inf),
                              "DTR" = rep(data[i, "DTR"], n_inf),
                              "gametocytemia" = rep(data[i, "gametocytemia"], n_inf))
      out <- rbind(out, place_inf)
    }
    if(n_un > 0){
      place_un <- data.frame("Experiment" = rep(data[i, "Experiment"], n_un),
                             "Cup" = rep(data[i, "Cup.ID"], n_un),
                             "DPI" = rep(data[i, "DPI"], n_un),
                             "presence" = rep(0, n_un),
                             "temp" = rep(data[i, "temp"], n_un),
                             "DTR" = rep(data[i, "DTR"], n_un),
                             "gametocytemia" = rep(data[i, "gametocytemia"], n_un))
      out <- rbind(out, place_un)
    }
  }
  
  # checking the numbers match up
  if(nrow(out) == sum(data$Dissected) &  sum(data$Infected) == sum(out$presence)){
    return(out)
  } else{
    return(print("numbers do not match up"))
  }
}

### generate oocyst intensity data quantiles
generate_oocyst_intensity_summary <- function(data){
  out <- unique(data[,c("DPI", "temp", "gametocytemia")])
  out$mean <- rep(NA, nrow(out))
  out$median <- rep(NA, nrow(out))
  out$lower <- rep(NA, nrow(out))
  out$upper <- rep(NA, nrow(out))
  
  for(i in 1:nrow(out)){
    placeholder <- subset(data, DPI == out[i,"DPI"] & temp == out[i, "temp"] & gametocytemia == out[i, "gametocytemia"])
    out[i, "mean"] <- mean(placeholder[,"Oocyst_number"])
    out[i, "median"] <- quantile(placeholder[,"Oocyst_number"], 0.5)
    out[i, "lower"] <- quantile(placeholder[,"Oocyst_number"], 0.025)
    out[i, "upper"] <- quantile(placeholder[,"Oocyst_number"], 0.975)
  }
  return(out)
}

##############################################
### get data into correct format for model ###
##############################################

# parasite data only
# scales the temperature so centred and the standard deviation is 1
# generates index for the temperature and gametocytemia
temp_g_indexing <- function(M_data, S_data){
  
  # temperature indexing
  temps <- append(M_data$temp, S_data$temp)
  U_temp <- sort(unique(temps))
  U_temp_scaled <- sort(unique(scale(temps, center = TRUE, scale = TRUE)[,1])) # converting temperature onto the same scale as used in the other model
  
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
             "oocyst_data" = M_data, "sporozoite_data" = S_data)
  
  rm(list = c("temps", "U_temp", "U_temp_scaled", "M_data", "S_data"))
  
  return(out)
}

#
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

oocyst_intensity_indexing <- function(data){
  U_data <- unique(data[,c("DPI", "Oocyst_number", "index_temp", "index_g", "index_gt")])
    
  index <- rep(NA, nrow(data))
  for(i in 1:nrow(U_data)){
    index[which(data$DPI == U_data[i, "DPI"] &
                  data$Oocyst_number == U_data[i, "Oocyst_number"] &
                  data$index_temp == U_data[i, "index_temp"] &                    
                  data$index_g == U_data[i, "index_g"] &
                  data$index_gt == U_data[i, "index_gt"])] <- i
  }
  out <- list("unique_oocyst_intensity" = U_data, "oocyst_intensity_index" = index)
  rm(list = c("U_data", "index"))
  return(out)
}

#########################################
### posterior predictive distribution ###
#########################################

# function that extracts the posterior predictive means from the MCMC output
prop_ppd_function <- function(fit_df, unique_temp, length_ppd_times, iterations, warmup, chains, length_unique_temp, n_gt_index, unique_gt, Stan_data_name){
  prop_ppd <- array(NaN, c(length_ppd_times, ((iterations - warmup) * chains), n_gt_index))
  for(i in 1:length_ppd_times){
    for(j in 1:n_gt_index){
      prop_ppd[i,,j] <- fit_df[,paste0(Stan_data_name,"[",i,",",j,"]")]
    }
  }
  prop_ppd_df <- data.frame()
  labs_g_ind <- c()
  labs_temp <- c()
  labs_gt_ind <- c()
  for(i in 1:n_gt_index){
    place <- as.data.frame(prop_ppd[,,i])
    prop_ppd_df <- rbind(prop_ppd_df, place)
    
    labs_g_ind <- append(labs_g_ind, rep(unique_gt[i, "index_g"], length_ppd_times))
    labs_temp <-  append(labs_temp, rep(unique_temp[unique_gt[i, "index_temp"]], length_ppd_times))
    labs_gt_ind <- append(labs_gt_ind, rep(unique_gt[i, "index"], length_ppd_times))
  }
  
  day_post_inf <- rep(PPD_times, nrow(unique_gt))
  prop_ppd_df[,"temp"] <- labs_temp
  prop_ppd_df[,"index_g"] <- labs_g_ind
  prop_ppd_df[,"index_gt"] <- labs_gt_ind
  prop_ppd_df[,"DPI"] <- day_post_inf
  prop_ppd_df$median <- apply(prop_ppd_df[,1:((iterations - warmup)*chains)], 1, median)
  #prop_ppd_df$lower <- apply(prop_ppd_df[,1:((iterations - warmup)*chains)], 1, quantile, probs = c(0.025), na.rm)
  #prop_ppd_df$upper <- apply(prop_ppd_df[,1:((iterations - warmup)*chains)], 1, quantile, probs = c(0.975))
  #prop_ppd_df$mean <- apply(prop_ppd_df[,1:((iterations - warmup)*chains)], 1, mean)
  
  x <- ((iterations - warmup) * chains) 
  prop_ppd_df <- prop_ppd_df%>%gather("iteration", "value", 1:x)
  prop_ppd_df$percentile <- rep("iteration", nrow(prop_ppd_df))
  
  prop_quantile_ppd_df <- subset(prop_ppd_df, iteration == "V1")
  prop_quantile_ppd_df <- prop_quantile_ppd_df[,c("DPI", "temp", "index_g", "index_gt", "median")]
  
  return(prop_quantile_ppd_df)
}

