// mSOS: a multiscale stochastic model of sporogony - single temperature model
// Author: Isaac J Stopard
// Version: 0.01 
// Last updated: 28/10/2020
// Notes: single parameters for oocyst load, parasite developments rate and the transmission probability
//        no mosquito survival data is included
//        fits to oocyst intensity and sporozoite prevalence data

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
}

data{

  // sporozoite prevalence data
  int n_S;
  int S_sample[n_S]; // n
  int S_positive[n_S]; // X
  real S_time[n_S];
  
  // oocyst intensity data
  int n_O_intsy;
  int O_intsy_index[n_O_intsy];
  
  // uniqe oocyst intensity combinations
  int n_unq_O_intsy_index;
  real unq_O_intsy_time[n_unq_O_intsy_index];
  real unq_O_intsy[n_unq_O_intsy_index];

  // posterior predicitive distribution times
  int n_ppd_times;
  vector[n_ppd_times] PPD_times;
  
  // priors
  real rate_O_prior;
}

parameters{
  // development time parameters
  real<lower = 0> shape_O;
  real<lower = 0> rate_O;

  real<lower = 0>  shape_S;
  real<lower = 0> rate_S;
  
  // parasite intensity population parameters
  real<lower = 0> mu;
  real<lower = 0> k;
  
  // proportion of mosquitoes infected parameters (probability of viable infection)
  real<lower = 0, upper = 1> delta;
}

transformed parameters{
  // infection parameters
  // parasite development parameters
  real<lower = 0> mu_total_S = (rate_O * shape_S + rate_S * shape_O) / (rate_O * rate_S);
  real<lower = 0> sigma_sq_S  = (rate_O^2 * shape_S + rate_S^2 * shape_O) / (rate_O^2 * rate_S^2);
  real<lower = 0> shape_total_S = mu_total_S^2 / sigma_sq_S;
  real<lower = 0> rate_total_S = mu_total_S / sigma_sq_S;
  
  // oocyst intensity data
  vector[n_unq_O_intsy_index] O_intsy_probability_lookup;

  // to collect the results
  vector[n_S] S_likelihoods;
  vector[n_S] S_prevalence;
  vector[n_O_intsy] O_intsy_likelihoods;

  // oocyst intensity data
  for(i in 1:n_unq_O_intsy_index){
    O_intsy_probability_lookup[i] = parasite_intensity_PDF(unq_O_intsy[i], unq_O_intsy_time[i], 
                                            mu, k, shape_O, rate_O, delta);
  }
  
  for(i in 1:n_S){
    S_prevalence[i] = parasite_prevalence(S_time[i], mu, k, shape_total_S, rate_total_S);
    S_likelihoods[i] = binomial_lpmf(S_positive[i] | S_sample[i], (delta * S_prevalence[i]));
  }
  
  for(i in 1:n_O_intsy){
    O_intsy_likelihoods[i] = O_intsy_probability_lookup[O_intsy_index[i]];
  }
}

model{
  // priors
  // development time parameters
  shape_O ~ normal(13.7, 0.05);
  rate_O ~ normal(rate_O_prior, 0.1);

  shape_S ~ normal(17.9, 2.0);
  rate_S ~ normal(2.2, 2.0); // fix the rate or shape - set to the mean from the stephensi analysis
  
  mu ~ lognormal(3, 1);  // increasing the SD of these priors causes divergent iterations (3.0, 1.5) // normal(5, 2.0)
  k ~ normal(5.0, 2.0); //normal(5, 2.0)// try a lognormal prior for mu 
  // need summary that compares the priors and the posteriors
  // can summarise - calculate the Kullback-Leibler divergence between the prior and the posterior
  // change the confidence intervals shown in the raw oocyst number data 
  // so it is more comparable to the mean
  // Kullback-Leibler divergence: need to write a function to calculate posterior density, 
  // then integrate the function - look up what to integrate - log(posterior) * log(prior)
  // helpful summary to gage how informative the priors are
  // maybe set a different prior for mu
  // do prior predictive distributions - for one value of temp / g combination
  
  // fix parameters to reasonable values and increase mu with increasing values
  // plot the predicted prevalence and intensity with the different values that are near to failing
  // 

  
  delta ~ beta(2.5, 2.5); // use a beta prior // normal(0.5, 2.0)
  
  // likelihood
  target += (sum(S_likelihoods) + sum(O_intsy_likelihoods));
}

generated quantities{
  real O_prevalence_ppd[n_ppd_times];
  real S_prevalence_ppd[n_ppd_times];
  real O_intsy_ppd[n_ppd_times];
  
  real mean_O = mu / (1 - (k / (k + mu))^k); 
  vector[n_S + n_O_intsy] log_lik;
  
  log_lik[1:n_S] = S_likelihoods;
  log_lik[(n_S+1):(n_S + n_O_intsy)] = O_intsy_likelihoods;
  
  for(i in 1:n_ppd_times){
    O_prevalence_ppd[i] = parasite_prevalence(PPD_times[i], mu, k, shape_O, rate_O) * delta;
    S_prevalence_ppd[i] = parasite_prevalence(PPD_times[i], mu, k, shape_total_S, rate_total_S) * delta;
    O_intsy_ppd[i] = delta * mean_O * regularized_gamma(shape_O, 0, (PPD_times[i] * rate_O));
  }
}
