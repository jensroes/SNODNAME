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
  
  vector[N] name_H; // naming diversity

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
  
  vector[2] b_name; // for name_H effect before and after onset
  vector[2] d_name; // effect when naming age is 0

  // hyper prior for effects
  real b_name_mu;
  real<lower =0> b_name_sigma;

  real d_name_mu;
  real<lower =0> d_name_sigma;

	// Parameters for non-centering
	real<lower = 0> alpha_mu;
  real<lower = 0> alpha_sigma;	
	real alpha_raw;			// distributions

  
  // For random effects
  vector[nS] u; //subject intercepts
  real<lower=0> sigma_u;	// subj sd
  vector[nI] w; //image intercepts
  real<lower=0> sigma_w;	// image sd
}

transformed parameters{
	real alpha = alpha_mu + alpha_sigma * alpha_raw;
  vector[N] RE; // Random effects

  for(n in 1:N){
    RE[n] = u[subj[n]] + w[image[n]]; 
  } 
}



model {
  int n = 0;
  
  // Priors
  alpha_mu ~ cauchy(6, 2.5);
  alpha_sigma ~ cauchy(0, 2.5);
  alpha_raw ~ normal(0, 2);

  b_name_mu ~ cauchy(0, 2.5);
  b_name_sigma ~ cauchy(0, 2.5);
  b_name ~ normal(b_name_mu, b_name_sigma);

  d_name_mu ~ cauchy(0, 2.5);
  d_name_sigma ~ cauchy(0, 2.5);
  d_name ~ normal(d_name_mu, d_name_sigma);


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
      if(name_H[n] == 0){
        if(b == 1){
          logy[s, b] ~ normal(alpha + d_name[1] + RE[n], sigma); 
        }
        if(b > 1){
          logy[s, b] ~ normal(alpha + d_name[2] + phi * logy[s, b-1] + RE[n], sigma); 
        }
      }
      if(name_H[n] != 0){
        if(b == 1){
          logy[s, b] ~ normal(alpha + b_name[1]*name_H[n] + RE[n], sigma); 
        }
        if(b > 1){
          logy[s, b] ~ normal(alpha + b_name[2]*name_H[n] + phi * logy[s, b-1] + RE[n], sigma); 
        }
      }
    }
  }
}

generated quantities{
  vector[N] log_lik;
  vector[N] y_tilde;
  int n = 0;

  for(s in 1:nK){
  int nKey = nKeys[s];
    for(b in 1:nKey){
      n += 1;
      if(name_H[n] == 0){
        if(b == 1){
          log_lik[n] = normal_lpdf(logy[s, b] | alpha + d_name[1] + RE[n], sigma); 
          y_tilde[n] = normal_rng(alpha + d_name[1] + RE[n], sigma);
        }
        if(b > 1){
          log_lik[n] = normal_lpdf(logy[s, b] | alpha + d_name[2] + phi * logy[s, b-1] + RE[n], sigma); 
          y_tilde[n] = normal_rng(alpha + d_name[2] + phi * logy[s, b-1] + RE[n], sigma);
        }
      }
      if(name_H[n] != 0){
        if(b == 1){
          log_lik[n] = normal_lpdf(logy[s, b] | alpha + b_name[1]*name_H[n] + RE[n], sigma); 
          y_tilde[n] = normal_rng(alpha + b_name[1]*name_H[n] + RE[n], sigma);
        }
        if(b > 1){
          log_lik[n] = normal_lpdf(logy[s, b] | alpha + b_name[2]*name_H[n] + phi * logy[s, b-1] + RE[n], sigma); 
          y_tilde[n] = normal_rng(alpha + b_name[2]*name_H[n] + phi * logy[s, b-1] + RE[n], sigma);
        }
      }
    }
  }
}
