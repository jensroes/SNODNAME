/*
  Serial order of bigrams auto regression model (first order)
  Random intercepts for subject 
*/
  
data {
  int<lower=0> N;

  int<lower=1> nS;            //number of subjects
  int<lower=1> nI;            // number of images

	int<lower=1, upper=nS> subj[N];   //subject id
	int<lower=1, upper=nI> image[N];   //subject id
  
  int<lower=1> nK;            // length of next vector
  int<lower=1> nKeys[nK];     // numbers of keystroke for each word within ppt
  int<lower=1> maxB;          // Max number of bigrams for matrix
  matrix[nK, maxB] y;            //outcome: for each ppt/image one IKI per column (in order as produced)
  
  vector[N] spell_H; // naming diversity

}

transformed data{
  matrix[nK, maxB] logy;
  
  for(s in 1:nK){
    int nKey = nKeys[s];
    for(b in 1:nKey){ 
      logy[s,b] = log(y[s,b]);
    }
  }
}


parameters {
  real phi; // for autoregression
  real<lower=0> sigma;		// residual sd
  real<lower = 0> delta; 
	real<lower = .001, upper = .999> theta; // mixing proportion
	real<lower=0> sigma_diff;
  
  
  vector[2] b_spell; // for name_H effect before and after onset

  // hyper prior for effects
  real b_spell_mu;
  real<lower =0> b_spell_sigma;

	// Parameters for non-centering
	real<lower = 0> alpha_mu;
  real<lower = 0> alpha_sigma;	
	vector[2] alpha_raw;			// distributions

  
  // For random effects
  vector[nS] u; //subject intercepts
  real<lower=0> sigma_u;	// subj sd
  vector[nI] w; //image intercepts
  real<lower=0> sigma_w;	// image sd
}

transformed parameters{
	vector[2] log_theta;
  vector[N] RE; // Random effects
	real sigmap_e;
	real sigma_e;
	vector[2] alpha = alpha_mu + alpha_sigma * alpha_raw;

	log_theta[1] = log(theta);
  log_theta[2] = log1m(theta);

  sigmap_e = sigma + sigma_diff;
  sigma_e = sigma - sigma_diff;

  for(n in 1:N){
    RE[n] = u[subj[n]] + w[image[n]]; 
  } 
}



model {
  vector[2] lp_parts;
  int n = 0;
  
  // Priors
  alpha_mu ~ cauchy(6, 2.5);
  alpha_sigma ~ cauchy(0, 2.5);
  alpha_raw ~ normal(0, 2);

  b_spell_mu ~ cauchy(0, 2.5);
  b_spell_sigma ~ cauchy(0, 2.5);
  b_spell ~ normal(b_spell_mu, b_spell_sigma);

  delta ~ cauchy(0, 5);
  theta ~ beta(2, 2);
  sigma_diff ~ normal(0, 1);

  phi ~ normal(0, 1);
  sigma ~ cauchy(0, 2.5);
  
  // REs priors
  sigma_u ~ normal(0,2.5);
  u ~ normal(0, sigma_u); //subj random effects

  sigma_w ~ normal(0,2.5);
  w ~ normal(0, sigma_w); //image random effects


  // Likelihood	
  for(s in 1:nK){
  int nKey = nKeys[s];
    for(b in 1:nKey){
      n += 1;
      if(b == 1){
        lp_parts[1] = log_theta[1] + normal_lpdf( logy[s, b] | alpha[1] + (b_spell[1] + delta)*spell_H[n] + RE[n], sigmap_e); 
        lp_parts[2] = log_theta[2] + normal_lpdf( logy[s, b] | alpha[1] + b_spell[1]*spell_H[n] + RE[n], sigma_e); 
      }
      if(b > 1){
         logy[s, b] ~ normal(alpha[2] + b_spell[2]*spell_H[n] + phi * logy[s, b-1] + RE[n], sigma_e); 
     }
     target += log_sum_exp(lp_parts); 
    }
  }
}

generated quantities{
  vector[N] log_lik;
  vector[N] y_tilde;
  real<lower=0,upper=1> theta_tilde; 
  int n = 0;

  for(s in 1:nK){
  int nKey = nKeys[s];
    for(b in 1:nKey){
      n += 1;
      if(b == 1){
        log_lik[n] = log_sum_exp(
          log_theta[1] + normal_lpdf(logy[s, b] | alpha[1] + (b_spell[1] + delta)*spell_H[n] + RE[n], sigmap_e),
          log_theta[2] + normal_lpdf(logy[s, b] | alpha[1] + b_spell[1]*spell_H[n] + RE[n], sigma_e)); 
        theta_tilde = bernoulli_rng(theta); 
        if(theta_tilde) { 
          y_tilde[n] = normal_rng(alpha[1] + (b_spell[1] + delta)*spell_H[n] + RE[n], sigmap_e);
        }
        else{
          y_tilde[n] = normal_rng(alpha[1] + b_spell[1]*spell_H[n] + RE[n], sigma_e);
        }
      }
      if(b > 1){
        log_lik[n] = normal_lpdf(logy[s, b] | alpha[2] + b_spell[2]*spell_H[n] + phi * logy[s, b-1] + RE[n], sigma_e); 
        y_tilde[n] = normal_rng(alpha[2] + b_spell[2]*spell_H[n] + phi * logy[s, b-1] + RE[n], sigma_e);
      }
    }
  }
}
