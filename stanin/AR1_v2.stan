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
  vector[N] spRT;
  vector[N] spell_H;
  vector[N] word_length;
  vector[N] word_frequency;
  vector[N] letter_frq;
  vector[N] dig_frq_pre;
  vector[N] dig_frq;
  vector[N] dig_frq_cond;
  
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
  vector[2] b_spoken;
  vector[2] b_spell;
  vector[2] b_wlen;
  vector[2] b_freq;
  vector[2] b_dig_freq_init;
  vector[2] b_dig_freq_mid;
  vector[2] b_dig_cond_prob;
  vector[2] b_letter_freq;

  // hyper prior for effects
  real name_mu;
  real<lower =0> name_sigma;

  real spoken_mu;
  real<lower =0> spoken_sigma;

  real spell_mu;
  real<lower =0> spell_sigma;

  real wlen_mu;
  real<lower =0> wlen_sigma;

  real freq_mu;
  real<lower =0> freq_sigma;

  real dfi_mu;
  real<lower =0> dfi_sigma;

  real dfm_mu;
  real<lower =0> dfm_sigma;

  real dcp_mu;
  real<lower =0> dcp_sigma;
  
  real letter_mu;
  real<lower =0> letter_sigma;


  
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
	vector[2] alpha = alpha_mu + alpha_sigma * alpha_raw;
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

  name_mu ~ cauchy(0, 2.5);
  name_sigma ~ cauchy(0, 2.5);
  b_name ~ normal(name_mu, name_sigma);

  spoken_mu ~ cauchy(0, 2.5);
  spoken_sigma ~ cauchy(0, 2.5);
  b_spoken ~ normal(spoken_mu, spoken_sigma);

  spell_mu ~ cauchy(0, 2.5);
  spell_sigma ~ cauchy(0, 2.5);
  b_spell ~ normal(spell_mu, spell_sigma);
  
  wlen_mu ~ cauchy(0, 2.5);
  wlen_sigma ~ cauchy(0, 2.5);
  b_wlen ~ normal(wlen_mu, wlen_sigma);
  
  freq_mu ~ cauchy(0, 2.5);
  freq_sigma ~ cauchy(0, 2.5);
  b_freq ~ normal(freq_mu, freq_sigma);
  
  dfi_mu ~ cauchy(0, 2.5);
  dfi_sigma ~ cauchy(0, 2.5);
  b_dig_freq_init ~ normal(dfi_mu, dfi_sigma);

  dfm_mu ~ cauchy(0, 2.5);
  dfm_sigma ~ cauchy(0, 2.5);
  b_dig_freq_mid ~ normal(dfm_mu, dfm_sigma);

  dcp_mu ~ cauchy(0, 2.5);
  dcp_sigma ~ cauchy(0, 2.5);
  b_dig_cond_prob ~ normal(dcp_mu, dcp_sigma);

  letter_mu ~ cauchy(0, 2.5);
  letter_sigma ~ cauchy(0, 2.5);
  b_letter_freq ~ normal(letter_mu, letter_sigma);

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
        logy[s, b] ~ normal(alpha[1] + b_name[1]*name_H[n] + 
                                    b_spoken[1]*spRT[n] +
                                    b_spell[1]*spell_H[n] +
                                    b_wlen[1]*word_length[n] +
                                    b_freq[1]*word_frequency[n] +
                                    b_dig_freq_init[1]*dig_frq_pre[n] +
                                    b_dig_freq_mid[1]*dig_frq[n] +
                                    b_dig_cond_prob[1]*dig_frq_cond[n] +
                                    b_letter_freq[1]*letter_frq[n] +
        RE[n], sigma); 
      }
      if(b > 1){
        logy[s, b] ~ normal(alpha[2] + b_name[2]*name_H[n] + 
                                    b_spoken[2]*spRT[n] +
                                    b_spell[2]*spell_H[n] +
                                    b_wlen[2]*word_length[n] +
                                    b_freq[2]*word_frequency[n] +
                                    b_dig_freq_init[2]*dig_frq_pre[n] +
                                    b_dig_freq_mid[2]*dig_frq[n] +
                                    b_dig_cond_prob[2]*dig_frq_cond[n] +
                                    b_letter_freq[2]*letter_frq[n] +
        phi * logy[s, b-1] + RE[n], sigma); 
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
      if(b == 1){
        log_lik[n] = normal_lpdf(logy[s, b] | alpha[1] + 
                                    b_name[1]*name_H[n] + 
                                    b_spoken[1]*spRT[n] +
                                    b_spell[1]*spell_H[n] +
                                    b_wlen[1]*word_length[n] +
                                    b_freq[1]*word_frequency[n] +
                                    b_dig_freq_init[1]*dig_frq_pre[n] +
                                    b_dig_freq_mid[1]*dig_frq[n] +
                                    b_dig_cond_prob[1]*dig_frq_cond[n] + 
                                    b_letter_freq[1]*letter_frq[n] +
                                    RE[n], sigma); 
        y_tilde[n] = normal_rng(alpha[1] + 
                                    b_name[1]*name_H[n] + 
                                    b_spoken[1]*spRT[n] +
                                    b_spell[1]*spell_H[n] +
                                    b_wlen[1]*word_length[n] +
                                    b_freq[1]*word_frequency[n] +
                                    b_dig_freq_init[1]*dig_frq_pre[n] +
                                    b_dig_freq_mid[1]*dig_frq[n] +
                                    b_dig_cond_prob[1]*dig_frq_cond[n] +
                                    b_letter_freq[1]*letter_frq[n] +
                                    RE[n], sigma);
      }
      if(b > 1){
        log_lik[n] = normal_lpdf(logy[s, b] | alpha[2] +
                                    b_name[2]*name_H[n] + 
                                    b_spoken[2]*spRT[n] +
                                    b_spell[2]*spell_H[n] +
                                    b_wlen[2]*word_length[n] +
                                    b_freq[2]*word_frequency[n] +
                                    b_dig_freq_init[2]*dig_frq_pre[n] +
                                    b_dig_freq_mid[2]*dig_frq[n] +
                                    b_dig_cond_prob[2]*dig_frq_cond[n] +
                                    b_letter_freq[2]*letter_frq[n] +
                                    phi * logy[s, b-1] + RE[n], sigma); 
        y_tilde[n] = normal_rng(alpha[2] + 
                                    b_name[2]*name_H[n] + 
                                    b_spoken[2]*spRT[n] +
                                    b_spell[2]*spell_H[n] +
                                    b_wlen[2]*word_length[n] +
                                    b_freq[2]*word_frequency[n] +
                                    b_dig_freq_init[2]*dig_frq_pre[n] +
                                    b_dig_freq_mid[2]*dig_frq[n] +
                                    b_dig_cond_prob[2]*dig_frq_cond[n] +
                                    b_letter_freq[2]*letter_frq[n] +
                                    phi * logy[s, b-1] + RE[n], sigma);
      }
    }
  }
}
