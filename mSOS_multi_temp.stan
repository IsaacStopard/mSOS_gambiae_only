// mSOS: a multiscale stochastic model of sporogony - all temperature model
// Author: Isaac J Stopard
// Version: 0.01 
// Last updated: 28/10/2020
// Notes: oocyst load varies with temperature and gametocytemia treatment
//        no mosquito survival data is included
//        relationship between parasite development rate and temperature is linear
//        parasite development rate allowed to vary with oocyst load and temperature
//        relationship between the human to mosquito transmission probability and temperature is logistic-quadratic
//        human to mosquito transmission probability allowed to vary with oocyst load and temperature
//        no random effects
//        fits to oocyst intensity and sporozoite prevalence data
//        fits separate parameters for mu and delta for each temperature / gametocyte treatment

functions{
  // incomplete regularised gamma function
  real regularized_gamma(real a, real z_0, real z_1){
    real result = gamma_p(a, z_1) - gamma_p(a, z_0);
    return result;
  }
  
  // Parasite prevalence estimate (oocyst or sporozoite) 
  // inputs: t (time), mu_NB (negative binomial distribution mean parameter), k (negative binomial distribution overdispersion parameter), 
  //             a (gamma distribution shape parameter) and b (gamma distribution rate parameter)
  real parasite_prevalence(real t, real mu_NB, real k, real a, real b){
    real result;
    result = (1 - (k^k) * ((k + mu_NB * regularized_gamma(a, 0, t * b))^-k)) / (1 - (k / (k + mu_NB))^k);
    return result;
  }
  
  // log PDF for the oocyst intensity data
  // inputs: intensity (observed oocyst load within individual mosquitoes), t (time), mu_NB (negative binomial distribution mean parameter),
  //         k (negative binomial distribution overdispersion parameter), a (gamma distribution shape parameter), b (gamma distribution rate parameter),
  //         delta (probability of viable infection)

  real parasite_intensity_PDF(real intensity, real t, real mu, real k, real a, real b, real delta){
    real result;
    real theta = regularized_gamma(a, 0, t * b);
    if(intensity == 0){
      result = log((1 - delta) + delta * ((k / (k + mu))^k - (k / (k + theta * mu))^k) / (-1 + (k / (k + mu))^k));
    } else{
      result = (log(delta) + k * log(k) + intensity * log(theta * mu) + (-intensity - k) * log(k + theta * mu) + lchoose(-1 + intensity + k, -1 + k)) - log(1 - (k/(k + mu))^k);
    }
    return result;
  }
  
  // functions to describe the relationship between temperature and different model parameters
  real logistic_function(real x){
    real result = 1 / (1 + exp(-x));
    return result;
  }
  
  real quadratic_function(real a, real b, real c, real temp){
    real result = (a * temp^2) + (b * temp) + c;
    return result;
  }
  
  real linear_function(real m, real c, real temp){
    real result = (m * temp) + c;
    return result;
  }
  
  real modified_logistic_function(real a, real m, real c, real b, real temp){
    real result = a * (exp(b * (m - temp))/(1 + exp(c * b * (m - temp))));
    return result;
  }
}

data{

  // sporozoite prevalence data
  int n_S;
  int S_sample[n_S]; // n
  int S_positive[n_S]; // X
  int S_ind_temp[n_S];
  real S_time[n_S];
  
  // oocyst intensity data
  int n_O_intsy;
  int O_intsy_index[n_O_intsy];
  
  // uniqe oocyst intensity combinations
  int n_unq_O_intsy_index;
  real unq_O_intsy_time[n_unq_O_intsy_index];
  real unq_O_intsy[n_unq_O_intsy_index];
  int unq_O_intsy_ind_temp[n_unq_O_intsy_index];

  // variable values for indexing
  // temperature
  int n_unq_temp;
  real unq_temp[n_unq_temp];
  // gametocytemia treatment

  // posterior predicitive distribution times
  int n_ppd_times;
  vector[n_ppd_times] PPD_times;

}

parameters{
  // development time parameters
  real<lower=0> shape_O;
  real m_rate_O;
  real c_rate_O;
  real a_rate_O;
 
  real<lower=0>  shape_S;
  real<lower=0> rate_S;
  
  // parasite intensity population parameters
  real a_mu;
  real m_mu;
  real c_mu;
  real<lower=0> k;

  // proportion of mosquitoes infected parameters (probability of viable infection)
  real a_delta;
  real m_delta;
  real c_delta;
  //real b_delta;
}

transformed parameters{
  // infection parameters
  real<lower=0> mu[n_unq_temp]; // varies with g and temp
  real<lower=0, upper=1> delta[n_unq_temp]; // varies with g and temp
  
  // parasite development parameters
  real<lower=0> rate_O[n_unq_temp]; // varies with g and temp
  real<lower=0> mu_total_S[n_unq_temp];
  real<lower=0> sigma_sq_S[n_unq_temp];
  real<lower=0> shape_total_S[n_unq_temp];
  real<lower=0> rate_total_S[n_unq_temp];
  
  // oocyst intensity data
  vector[n_unq_O_intsy_index] O_intsy_probability_lookup;

  // to collect the results
  vector[n_S] S_likelihoods;
  vector[n_S] S_prevalence;
  vector[n_O_intsy] O_intsy_likelihoods;
  
  // oocyst load
  for(i in 1:n_unq_temp){
    mu[i] = quadratic_function(a_mu, m_mu, c_mu, unq_temp[i]);
    delta[i] = logistic_function(quadratic_function(a_delta, m_delta, c_delta, unq_temp[i]));
    //delta[i] = modified_logistic_function(a_delta, m_delta, c_delta, b_delta, unq_temp[i]);
    // parasite development parameters
    rate_O[i] = quadratic_function(a_rate_O, m_rate_O, c_rate_O, unq_temp[i]);
    mu_total_S[i] = (rate_O[i] * shape_S + rate_S * shape_O) / (rate_O[i] * rate_S);
    sigma_sq_S[i] = (rate_O[i]^2 * shape_S + rate_S^2 * shape_O) / (rate_O[i]^2 * rate_S^2);
    shape_total_S[i] = mu_total_S[i]^2 / sigma_sq_S[i];
    rate_total_S[i] = mu_total_S[i] / sigma_sq_S[i];
  }
 
  // oocyst intensity data
  for(i in 1:n_unq_O_intsy_index){
    O_intsy_probability_lookup[i] = parasite_intensity_PDF(unq_O_intsy[i], unq_O_intsy_time[i], 
                                            mu[unq_O_intsy_ind_temp[i]], k, shape_O, rate_O[unq_O_intsy_ind_temp[i]], 
                                            delta[unq_O_intsy_ind_temp[i]]);
  }
  
  for(i in 1:n_S){
    S_prevalence[i] = parasite_prevalence(S_time[i], mu[S_ind_temp[i]], k, shape_total_S[S_ind_temp[i]], rate_total_S[S_ind_temp[i]]);
    S_likelihoods[i] = binomial_lpmf(S_positive[i] | S_sample[i], (delta[S_ind_temp[i]] * S_prevalence[i]));
  }
  
  for(i in 1:n_O_intsy){
    O_intsy_likelihoods[i] = O_intsy_probability_lookup[O_intsy_index[i]];
  }
}

model{
  // priors
  // development time parameters
  shape_O ~ normal(13.7, 1.0);
  m_rate_O ~ normal(0, 2.5);
  c_rate_O ~ normal(4.0, 2.5);
  a_rate_O ~ normal(0, 2.5);
  
  shape_S ~ normal(17.9, 0.25);
  rate_S ~ normal(2.2, 0.25);
  
  a_mu ~ normal(0, 2.5);
  m_mu ~ normal(0, 2.5);
  c_mu ~ normal(5, 2.5);
  k ~ normal(3, 2.5);
  
  a_delta ~ normal(0, 2.5);
  m_delta ~ normal(0, 2.5);
  c_delta ~ normal(0.75, 2.5);
  //b_delta ~ normal(0, 3.5);

  // likelihood
  target += (sum(S_likelihoods) + sum(O_intsy_likelihoods));
}

generated quantities{
  real O_prevalence_ppd[n_ppd_times, n_unq_temp];
  real S_prevalence_ppd[n_ppd_times, n_unq_temp];
  real O_intsy_ppd[n_ppd_times, n_unq_temp];
  
  real mean_O[n_unq_temp];
  vector[n_S + n_O_intsy] log_lik;
  
  log_lik[1:n_S] = S_likelihoods;
  log_lik[(n_S+1):(n_S + n_O_intsy)] = O_intsy_likelihoods;
  
  for(i in 1:n_unq_temp){
    mean_O[i] = mu[i] / (1 - (k / (k + mu[i]))^k); 
  }
  
  for(i in 1:n_ppd_times){
    for(j in 1:n_unq_temp){
      O_prevalence_ppd[i, j] = parasite_prevalence(PPD_times[i], mu[j], k, shape_O, rate_O[j]) * delta[j];
      
      S_prevalence_ppd[i, j] = parasite_prevalence(PPD_times[i], mu[j], k, shape_total_S[j], 
                                                    rate_total_S[j]) * delta[j];
      
      O_intsy_ppd[i, j] = delta[j] * mean_O[j] * regularized_gamma(shape_O, 0, (PPD_times[i] * rate_O[j]));
    }
  }
}
