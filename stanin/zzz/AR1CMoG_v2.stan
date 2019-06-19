/*
  Serial order of bigrams autoregression model
Mixture model of first order and second order autoregression model
Random intercepts for subject 
Conditional autocorrelation mixture
*/
  
data {
  int<lower=0> N;             // total number of observations
  int<lower=1> nS;              //number of subjects
  int nB[nS];                 // total number of bigrams produced by ppt
  int<lower=1> maxB;          // Max number of bigrams produced by ppt for matrix
  matrix[nS,maxB] y;            //outcome: for each subject one IKI per column (in order as produced)
  int<lower=1, upper=nS> subj[N];   //subject id
}

transformed data{
  matrix[nS, maxB] logy;
  
  for(s in 1:nS){
    int nBS = nB[s];
    for(b in 1:nBS){ 
      logy[s,b] = log(y[s,b]);
    }
  }
}

parameters {
  real phi;
  
  simplex[2] theta; // mixing proportion
  
  // Parameters for non-centering
  real<lower = 0> alpha_mu;
  real<lower = 0> alpha_sigma;	
  real alpha_raw;			// distributions
  
  real<lower=0> sigma;		// residual sd
  real<lower=0> sigma_diff;
  
  // For random effects
  vector[nS] u; //subject intercepts
  real<lower=0> sigma_u;	// subj sd
  
}


transformed parameters{
  real alpha = alpha_mu + alpha_sigma * alpha_raw;
  vector[2] log_theta;
  vector[N] RE; // Random effects
  real sigmap_e;
  real sigma_e;
  
  log_theta[1] = log(theta[1]);
  log_theta[2] = log1m(theta[1]);
  
  sigmap_e = sigma + sigma_diff;
  sigma_e = sigma - sigma_diff;
  for(n in 1:N){
    RE[n] = u[subj[n]];
  }
}

model {
  int n = 0;
  vector[2] lp_parts;
  
  // Priors
  alpha_mu ~ cauchy(6, 2.5);
  alpha_sigma ~ cauchy(0, 2.5);
  alpha_raw ~ normal(0, 2);
  phi ~ normal(0, 1);
  theta ~ beta(2, 2);
  sigma ~ cauchy(0, 2.5);
  sigma_diff ~ normal(0, 1);
  
  // REs priors
  sigma_u ~ normal(0, 2.5);
  u ~ normal(0, sigma_u); //subj random effects
  
  
  // Likelihood	
  for(s in 1:nS){
    int nBS = nB[s];
    for(b in 1:nBS){ 
      n += 1;
      lp_parts[1] = log_theta[1] + normal_lpdf(logy[s,b] | alpha + RE[n], sigmap_e);
      if(b == 1){
        lp_parts[2] = log_theta[2] + normal_lpdf(logy[s,b] | alpha + RE[n], sigma_e);
      }
      if(b > 1){
        lp_parts[2] = log_theta[2] + normal_lpdf(logy[s,b] | alpha + phi*logy[s, b-1] + RE[n], sigma_e);
      }
      target += log_sum_exp(lp_parts); 
    }
  }
}

generated quantities{
  vector[N] log_lik;
  vector[N] y_tilde;
  vector[2] lp_parts;
  real<lower=0,upper=1> theta_tilde; 
  int n = 0;
  
  // likelihood: 
    for(s in 1:nS){
      int nBS = nB[s];
      for(b in 1:nBS){
        n += 1;
        lp_parts[1] = log_theta[1] + normal_lpdf(logy[s,b] | alpha + RE[n], sigmap_e);
        if(b == 1){
          lp_parts[2] = log_theta[2] + normal_lpdf(logy[s,b] | alpha + RE[n], sigma_e);
        }
        if(b > 1){      
          lp_parts[2] = log_theta[2] + normal_lpdf(logy[s,b] | alpha + phi*logy[s, b-1] + RE[n], sigma_e);
        }
        log_lik[n] = log_sum_exp(lp_parts);
        theta_tilde = bernoulli_rng(theta[1]); 
        if(theta_tilde) {
            y_tilde[n] = normal_rng(alpha + RE[n], sigmap_e);
          }
        if(!theta_tilde){
          if(b == 1){
            y_tilde[n] = normal_rng(alpha + RE[n], sigma_e);    
          }
          if(b > 1){
            y_tilde[n] = normal_rng(alpha + phi*logy[s, b-1] + RE[n], sigma_e);
          }
        }
      }
    }
}

