# Helper functions to process the raw data
# Author: Isaac J Stopard
# Version: 1.0.0 

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