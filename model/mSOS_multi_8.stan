// mSOS: a multiscale stochastic model of sporogony - pooled model
// Author: Isaac J Stopard
// Version: 0.01 
// Last updated: 03/12/2020
// Notes: single parameters for oocyst load, parasite developments rate and the transmission probability
//        no mosquito survival data is included
//        fits to oocyst intensity and sporozoite prevalence data
//        single shape_O and rate_O for all temperature and g combinations
//        relationship between shape_S / rate_S and temperature / g
//        relationship between delta and temperature / g
//        relationship between k and g
//        relationship between mu and temperature / g


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
  real parasite_intensity_PDF_unstable(real intensity, real t, real mu_NB, real k, real a, real b, real delta){
    real result;
    real theta = regularized_gamma(a, 0, t * b);
    if(intensity == 0){
      result = log((1 - delta) + delta * ((k / (k + mu_NB))^k - (k / (k + theta * mu_NB))^k) / (-1 + (k / (k + mu_NB))^k));
    } else{
      result = log(delta * ((k^k * (theta * mu_NB)^intensity * (k + theta * mu_NB)^(-intensity - k) * exp(lchoose(-1 + intensity + k, -1 + k))) / (1 - (k / (k + mu_NB))^k)));
    }
    return result;
  }
  
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
  
  // rewrite the power bit - this
  // built a duplicate and check they line up
  
  real sporozoite_prevalence_PDF(real t, real mu, real k, real a, real b, real delta, int S_positive, int S_sample){
    real S_prevalence = parasite_prevalence(t, mu, k, a, b);
    real S_likelihoods = binomial_lpmf(S_positive | S_sample, (delta * S_prevalence)); 
    return S_likelihoods;
  }
  
  real quadratic_function(real t, real a, real b, real c){
    real out = (t^2 * a) + (t * b) + c;
    return out;
  }
  
  real logistic_function(real c){
    real out = 1 / (1 + exp(-c));
    return out;
  }
  
  real power_ricker_function(real t, real a, real b, real c){
    real out = b * ((t^2)/a) * exp(1 - ((t^2)/a))^c;
    return out;
  }
}

data{

  // sporozoite prevalence data
  int n_S;
  int S_sample[n_S]; // n
  int S_positive[n_S]; // X
  real S_time[n_S];
  int S_ind_temp[n_S];
  
  // oocyst intensity data
  int n_O_intsy;
  int O_intsy_index[n_O_intsy];
  
  // uniqe oocyst intensity combinations
  int n_unq_O_intsy_index;
  real unq_O_intsy_time[n_unq_O_intsy_index];
  real unq_O_intsy[n_unq_O_intsy_index];
  int unq_O_intsy_ind_temp[n_unq_O_intsy_index];

  // posterior predicitive distribution times
  int n_ppd_times;
  vector[n_ppd_times] PPD_times;
  
  int n_ppd_times_O;
  vector[n_ppd_times_O] PPD_times_O;
  
  //
  int n_unq_temp;
  real temp[n_unq_temp];
 
  // priors
  real rate_O_prior;
  
}

parameters{
  // development time parameters
  real<lower = 0> shape_O;
  real<lower = 0> rate_O;
  
  real a_shape_S;
  real b_shape_S;
  real c_shape_S;
  
  real m_rate_S;
  real c_rate_S;
  
  real a_delta;
  real b_delta;
  real c_delta;
  
  real a_delta_S;
  real b_delta_S;
  real c_delta_S;
  
  // parasite intensity population parameters
  real a_mu;
  real b_mu;
  real c_mu;
  real g_mu;
  //real<lower = 0> mu_[n_unq_gt];
  
  //real<lower = 0> k_[n_unq_g];
  real<lower = 0> k;
  
}

transformed parameters{
  // infection parameters
  // parasite development parameters
  real<lower = 0> shape_S[n_unq_temp];
  real<lower = 0> rate_S[n_unq_temp];
  
  real<lower = 0> mu_total_S[n_unq_temp];
  real<lower = 0> sigma_sq_S[n_unq_temp];
  real<lower = 0> shape_total_S[n_unq_temp];
  real<lower = 0> rate_total_S[n_unq_temp];
  
  real<lower = 0, upper = 1> delta[n_unq_temp];
  // proportion of mosquitoes infected parameters (probability of viable infection)
  real<lower = 0, upper = 1> delta_S[n_unq_temp];
  
  //real<lower = 0> mu[n_unq_gt];
  real<lower = 0> mu[n_unq_temp];
  //real<lower = 0> k[n_unq_gt];
  
  // oocyst intensity data
  vector[n_unq_O_intsy_index] O_intsy_probability_lookup;

  // to collect the results
  vector[n_S] S_likelihoods;
  vector[n_S] S_prevalence;
  vector[n_O_intsy] O_intsy_likelihoods;
  
  for(i in 1:n_unq_temp){
    shape_S[i] = quadratic_function(temp[i], a_shape_S, b_shape_S, c_shape_S);
    rate_S[i] = temp[i] * m_rate_S + c_rate_S;
    mu_total_S[i] = (rate_O * shape_S[i] + rate_S[i] * shape_O) / (rate_O * rate_S[i]);
    sigma_sq_S[i]  = (rate_O^2 * shape_S[i] + rate_S[i]^2 * shape_O) / (rate_O^2 * rate_S[i]^2);
    shape_total_S[i] = mu_total_S[i]^2 / sigma_sq_S[i];
    rate_total_S[i] = mu_total_S[i] / sigma_sq_S[i];
    delta[i] = logistic_function(quadratic_function(temp[i], a_delta, b_delta, c_delta));
    delta_S[i] = logistic_function(quadratic_function(temp[i], a_delta_S, b_delta_S, c_delta_S));
    mu[i] = logistic_function(quadratic_function(temp[i], a_mu, b_mu, c_mu)) * g_mu;
    //mu[i] = a_mu * exp(-((temp[i] - b_mu)^2 / (2 * c_mu)^2));
    //k[i] = k_[g_unq_gt[i]];
  }
  
  // oocyst intensity data
  for(i in 1:n_unq_O_intsy_index){
    O_intsy_probability_lookup[i] = parasite_intensity_PDF(unq_O_intsy[i], unq_O_intsy_time[i], 
                                            mu[unq_O_intsy_ind_temp[i]], k, shape_O, 
                                            rate_O, delta[unq_O_intsy_ind_temp[i]]);
  }
  
  for(i in 1:n_S){
    S_prevalence[i] = parasite_prevalence(S_time[i], mu[S_ind_temp[i]], k, shape_total_S[S_ind_temp[i]], rate_total_S[S_ind_temp[i]]);
    S_likelihoods[i] = binomial_lpmf(S_positive[i] | S_sample[i], (delta[S_ind_temp[i]] * delta_S[S_ind_temp[i]] * S_prevalence[i]));
  }
  
  for(i in 1:n_O_intsy){
    O_intsy_likelihoods[i] = O_intsy_probability_lookup[O_intsy_index[i]];
  }
}

model{
  // priors
  // development time parameters
  
  shape_O ~ normal(13.7, 0.1);
  rate_O ~ normal(rate_O_prior, 0.1);
  a_shape_S ~ normal(0, 5.0);
  b_shape_S ~ normal(0, 5.0);
  c_shape_S ~ normal(20, 4.0); // (20, 4) 
  
  m_rate_S ~ normal(0, 5.0);
  c_rate_S ~ normal(2.2, 4.0); // (2.2, 4)
  
  a_delta ~ normal(0, 5.0);
  b_delta ~ normal(0, 5.0);
  c_delta ~ normal(0, 4.0); // (10, 4)
  
  a_delta_S ~ normal(0, 5.0);
  b_delta_S ~ normal(0, 5.0);
  c_delta_S ~ normal(0, 5.0); // assumes no difference between O and S prevalence
  
  a_mu ~ normal(0, 5.0);
  b_mu ~ normal(0, 5.0);
  c_mu ~ normal(10, 5.0); // (10, 5)
  g_mu ~ normal(8, 4.0); // mean number of oocysts 8
  
 
  k ~ lognormal(1, 0.75);
  
  // likelihood
  target += (sum(S_likelihoods) + sum(O_intsy_likelihoods));
}

generated quantities{
  //real O_prevalence_ppd[n_ppd_times, n_unq_temp];
  real S_prevalence_ppd[n_ppd_times, n_unq_temp];
  real O_intsy_ppd[n_ppd_times_O, n_unq_temp];
  real O_intsy_inf_ppd[n_ppd_times_O, n_unq_temp];
  
  real<lower = 0> mean_O[n_unq_temp];
  vector[n_S + n_O_intsy] log_lik;
  
  for(i in 1:n_unq_temp){
    mean_O[i] = mu[i] / (1 - (k / (k + mu[i]))^k); 
  }
  
  
  log_lik[1:n_S] = S_likelihoods;
  log_lik[(n_S+1):(n_S + n_O_intsy)] = O_intsy_likelihoods;
  
  for(i in 1:n_ppd_times_O){
    for(j in 1:n_unq_temp){
        //O_prevalence_ppd[i, j] = parasite_prevalence(PPD_times[i], mu[j], k[j], shape_O, rate_O) * delta[j];
        O_intsy_ppd[i, j] = delta[j] * mean_O[j] * regularized_gamma(shape_O, 0, (PPD_times_O[i] * rate_O));
        O_intsy_inf_ppd[i, j] = mean_O[j] * regularized_gamma(shape_O, 0, (PPD_times_O[i] * rate_O));
    }
  }
  
  for(i in 1:n_ppd_times){
    for(j in 1:n_unq_temp){
        S_prevalence_ppd[i, j] = parasite_prevalence(PPD_times[i], mu[j], k, shape_total_S[j], rate_total_S[j]) * delta[j] * delta_S[j];
    }
  }
}
