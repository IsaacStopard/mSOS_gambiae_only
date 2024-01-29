# EIP PDF function
EIP_PDF <- function(t, a, b, mu, k){
  return(-(((1/b)^-a) * exp(-b * t) * (t^(-1+a)) * mu * ((k / (k + mu))^(1 + k)) * (((k + mu * zipfR::Rgamma(a, b*t)) / (k + mu)) ^ (-1 - k))
           / ((-1 + (k / (k + mu))^k) * gamma(a))))
}

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

# Step 1
# probability of surviving the EIP given infection at age a0 and the Gompertz survival function
# integrand
# check if this is vectorised
prop_surv_gomp_int <- function(a, # time since age at infection
                               a0, # age at infection
                               g1, g2, # gompertz function parameters
                               shape, rate, mu_PL, k # EIP PDF parameters)
){
  return(EIP_PDF(t = a, a = shape, b = rate, mu = mu_PL, k = k) * # EIP PDF
           surv_gomp(a0 = a0, a = a, g1 = g1, g2 = g2))
}

v_prop_surv_gomp_int <- Vectorize(prop_surv_gomp_int)

# integral
prop_surv_gomp <- function(a0, # age at infection
                           g1, 
                           g2, # gompertz function parameters
                           shape, 
                           rate, 
                           mu_PL, 
                           k # EIP PDF parameters
){
  out <- calculus::integral(v_prop_surv_gomp_int, # mosquito survival
                            bounds = list("a" = c(0, Inf)),
                            params = list("a0" = a0, 
                                          "g1" = g1, 
                                          "g2" = g2,
                                          "shape" = shape, 
                                          "rate" = rate, 
                                          "mu_PL" = mu_PL, 
                                          "k" = k),
                            method = "cuhre",
                            coordinates = "cartesian",
                            vectorize = FALSE,
                            drop = TRUE
  )
  return(out$value)
}

v_prop_surv_gomp <- Vectorize(prop_surv_gomp)

prop_surv_gomp_fixed <- function(a, # degree-day model EIP
                                 a0, # age at infection
                                 g1, 
                                 g2){
  return(surv_gomp(a0 = a0, a = a, g1 = g1, g2 = g2))
}

v_prop_surv_gomp_fixed <- Vectorize(prop_surv_gomp_fixed)

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

v_p_j_bites_EIP_a1_int <- Vectorize(p_j_bites_EIP_a1_int)

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
  
  out <- calculus::integral(v_p_j_bites_EIP_a1_int, 
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

v_p_j_bites_EIP_a1 <- Vectorize(p_j_bites_EIP_a1)

# expected number of bites
E_j_bites_EIP_a1_int <- function(j, # number of bites
                                 alpha, # Poisson distribution rate
                                 a1, # age the mosquito exits the EIP
                                 g1,
                                 g2){
  
  return(j * v_p_j_bites_EIP_a1(j = j, alpha = alpha, a1 = a1, g1 = g1, g2 = g2))
  
}

v_E_j_bites_EIP_a1_int <- Vectorize(E_j_bites_EIP_a1_int)

E_j_bites_EIP_a1 <- function(alpha, # Poisson distribution rate
                             a1, # age the mosquito exits the EIP
                             g1,
                             g2){
  
  out <- calculus::integral(
    v_E_j_bites_EIP_a1_int,
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

v_E_j_bites_EIP_a1 <- Vectorize(E_j_bites_EIP_a1)

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
  out <- v_p_j_bites_EIP_a1(j = j, 
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

v_p_j_bites_meal_a0_int <- Vectorize(p_j_bites_meal_a0_int)

p_j_bites_meal_a0 <- function(j,
                              a0, # age the mosquito is infected
                              alpha, # Poisson distribution rate
                              g1,
                              g2,
                              shape, 
                              rate, 
                              mu_PL, 
                              k){
  
  out <- calculus::integral(v_p_j_bites_meal_a0_int, # mosquito survival
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

v_p_j_bites_meal_a0 <- Vectorize(p_j_bites_meal_a0)

# degree-day model equivalent

p_j_bites_meal_a0_fixed <- function(j,
                                    a0, # age the mosquito is infected
                                    alpha, # Poisson distribution rate
                                    g1,
                                    g2,
                                    EIP# EIP PDF parameters
){
  a1 <- a0 + EIP # age the mosquito exits the EIP
  out <- v_p_j_bites_EIP_a1(j = j, 
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

v_p_j_bites_meal_a0_fixed <- Vectorize(p_j_bites_meal_a0_fixed)

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
  return(j * v_p_j_bites_meal_a0(j = j, 
                               a0 = a0, # age the mosquito is infected
                               alpha = alpha, # Poisson distribution rate
                               g1 = g1,
                               g2 = g2,
                               shape = shape, 
                               rate = rate, 
                               mu_PL = mu_PL, 
                               k = k))
  
}

v_E_j_bites_a0_int <- Vectorize(E_j_bites_a0_int)

E_j_bites_a0 <- function(a0, # age the mosquito exits the EIP
                         alpha, # Poisson distribution rate
                         g1,
                         g2,
                         shape,
                         rate,
                         mu_PL,
                         k){
  
  out <- calculus::integral(v_E_j_bites_a0_int,
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

v_E_j_bites_a0 <- Vectorize(E_j_bites_a0)



E_j_bites_a0_int_fixed <- function(j,
                                   a0, # age the mosquito exits the EIP
                                   alpha, # Poisson distribution rate
                                   g1,
                                   g2,
                                   EIP){
  
  return(j * v_p_j_bites_meal_a0_fixed(j = j, 
                                     a0 = a0, # age the mosquito is infected
                                     alpha = alpha, # Poisson distribution rate
                                     g1 = g1,
                                     g2 = g2,
                                     EIP = EIP))
  
}

v_E_j_bites_a0_int_fixed <- Vectorize(E_j_bites_a0_int_fixed)

E_j_bites_a0_fixed <- function(a0, # age the mosquito exits the EIP
                               alpha, # Poisson distribution rate
                               g1,
                               g2,
                               EIP){
  
  out <- calculus::integral(v_E_j_bites_a0_int_fixed,
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

v_E_j_bites_a0_fixed <- Vectorize(E_j_bites_a0_fixed)

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
