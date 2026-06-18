// This model does not run well for Stenella spp. 
// Stan returned the following warning:
// Warning: 8000 of 8000 (100.0%) transitions hit the maximum treedepth limit of 10.
//
// See https://mc-stan.org/misc/warnings for details.Warning: 3 of 4 chains had an E-BFMI less than 0.3.
//
// See https://mc-stan.org/misc/warnings for details.
// Also ESS values were very low.
//
// Here are problems per Gemini:
// 1. The Log-Scale Explosion (Why Treedepth hit 100%)In JAGS, a dnorm(0, 0.001) prior meant a mean of 0 and a variance of 1,000. In Stan, I translated that to normal(0, 10) (mean of 0, standard deviation of 10).However, a1, alpha1, and a2 are drawn from a lognormal distribution. If the hyper-mean (mu_a1) explores a value of just 8, the actual growth rate it passes to the equation is $e^8 = 2,980$. A value of -8 becomes $0.0003$.Because a standard deviation of 10 on the log scale is astronomically huge, Stan was trying to calculate gradients across a likelihood surface that spanned from microscopic fractions to billions. It had to take the maximum number of tiny steps (1,024 steps, or treedepth 10) for every single iteration just to cross that massive, flat mathematical plain.The Fix: We must tighten the log-scale priors. Growth rates in the Laird model are usually between 0.01 and 1.0. On the log scale, a prior of normal(-1, 2) perfectly covers this biological reality without sending the sampler to infinity.

// 2. The Hierarchical Funnel (Why E-BFMI was low)You used a "Centered Parameterization" for your hierarchical draws:a1[s, x] ~ lognormal(mu_a1[x], sd_a1);In Stan, when group-level data is noisy, this creates a literal funnel shape in the posterior probability space. The sampler slides down into the narrow neck of the funnel, gets trapped, and throws a "Low E-BFMI" (Energy Bayesian Fraction of Missing Information) warning.The Fix: We use the Non-Centered Parameterization (NCP). Instead of drawing the parameter directly from the distribution, we draw a raw, unit-normal variable (mean 0, SD 1) and mathematically transform it. This is the single most powerful trick in all of Stan coding.

// The upgraded Stan model is called model_growth_model_ncp.stan


data {
  int<lower=1> N_growth;
  int<lower=1> N_species;
  int<lower=1> N_sexes;
  
  vector[N_growth] length_growth;
  vector[N_growth] age_growth;
  array[N_growth] int<lower=1, upper=N_species> sp_growth;
  array[N_growth] int<lower=1, upper=N_sexes> sex_growth;
  
  // Note: Pass STANDARD DEVIATION in your R data list, not precision!
  vector[N_species] mu_L0_data;
  vector[N_species] sd_L0_data; 
  
  // Hard switch: 1 for two-phase, 0 for single-phase
  array[N_species] int<lower=0, upper=1> has_phase2; 
}

parameters {
  // Bounding parameters here handles the "truncation" natively
  vector<lower=0>[N_species] L0;
  
  matrix<lower=0>[N_species, N_sexes] a1;
  matrix<lower=0>[N_species, N_sexes] alpha1;
  matrix<lower=0>[N_species, N_sexes] tc;
  
  matrix<lower=0>[N_species, N_sexes] a2;
  matrix<lower=150, upper=260>[N_species, N_sexes] L_inf;

  // Global Hyper-means
  vector[N_sexes] mu_a1;
  vector[N_sexes] mu_alpha1;
  vector<lower=0>[N_sexes] mu_tc;
  vector<lower=150, upper=260>[N_sexes] mu_L_inf;
  matrix[N_species, N_sexes] mu_a2;

  // Global Hyper-Standard Deviations
  real<lower=0, upper=10> sd_a1;
  real<lower=0, upper=10> sd_alpha1;
  real<lower=0, upper=10> sd_tc;
  real<lower=0, upper=10> sd_a2;
  real<lower=0, upper=50> sd_L_inf;

  // Residual Standard Deviations
  real<lower=0, upper=30> sigma_juv;
  real<lower=0, upper=50> sigma_adult;
}

transformed parameters {
  matrix[N_species, N_sexes] Lc;
  matrix[N_species, N_sexes] alpha2;

  // Calculate Lc and deterministic alpha2
  for (s in 1:N_species) {
    for (x in 1:N_sexes) {
      Lc[s, x] = L0[s] * exp((a1[s, x] / alpha1[s, x]) * (1.0 - exp(-alpha1[s, x] * tc[s, x])));
      
      // fmax is Stan's version of the max() safety net
      real denom = fmax(log(L_inf[s, x]) - log(Lc[s, x]), 0.001);
      alpha2[s, x] = a2[s, x] / denom;
    }
  }
}

model {
  // ===========================================================
  // PRIORS
  // ===========================================================
  // L0 prior using input data
  L0 ~ normal(mu_L0_data, sd_L0_data);

  // Hyper-means (using SD of 10-31 instead of JAGS precision 0.01-0.001)
  mu_a1 ~ normal(0, 10);
  mu_alpha1 ~ normal(0, 10);
  mu_tc ~ normal(6, 10);
  mu_L_inf ~ normal(200, 30);
  
  for (s in 1:N_species) {
    for (x in 1:N_sexes) {
      mu_a2[s, x] ~ normal(0, 10);
    }
  }

  // Hierarchical Parameter Draws
  for (s in 1:N_species) {
    for (x in 1:N_sexes) {
      // Stan lognormal takes the mean and SD on the log scale
      a1[s, x] ~ lognormal(mu_a1[x], sd_a1);
      alpha1[s, x] ~ lognormal(mu_alpha1[x], sd_alpha1);
      
      // If the species has Phase 2, fit biological priors
      if (has_phase2[s] == 1) {
        tc[s, x] ~ normal(mu_tc[x], sd_tc);
        a2[s, x] ~ lognormal(mu_a2[s, x], sd_a2);
        L_inf[s, x] ~ normal(mu_L_inf[x], sd_L_inf);
      } else {
        // Dummy priors to keep the sampler happy for Species 2
        tc[s, x] ~ normal(10, 1);
        a2[s, x] ~ lognormal(0, 1);
        L_inf[s, x] ~ normal(200, 10);
      }
    }
  }

  // ===========================================================
  // LIKELIHOOD (OBSERVATION MODEL)
  // ===========================================================
  for (i in 1:N_growth) {
    int s = sp_growth[i];
    int x = sex_growth[i];
    real mu;
    real sigma_i;

    // Base Phase 1 calculation
    real mu1 = L0[s] * exp((a1[s, x] / alpha1[s, x]) * (1.0 - exp(-alpha1[s, x] * age_growth[i])));

    if (has_phase2[s] == 1) {
      if (age_growth[i] <= tc[s, x]) {
        mu = mu1;
        sigma_i = sigma_juv;
      } else {
        // Phase 2 calculation
        mu = Lc[s, x] * exp((a2[s, x] / alpha2[s, x]) * (1.0 - exp(-alpha2[s, x] * (age_growth[i] - tc[s, x]))));
        sigma_i = sigma_adult;
      }
    } else {
      // Single-phase overwrite
      mu = mu1;
      sigma_i = sigma_juv; 
    }

    // Likelihood calculation
    length_growth[i] ~ normal(mu, sigma_i);
  }
}

generated quantities {
  // Calculate pointwise log-likelihood for LOO-CV natively
  vector[N_growth] log_lik;
  
  for (i in 1:N_growth) {
    int s = sp_growth[i];
    int x = sex_growth[i];
    real mu;
    real sigma_i;

    real mu1 = L0[s] * exp((a1[s, x] / alpha1[s, x]) * (1.0 - exp(-alpha1[s, x] * age_growth[i])));

    if (has_phase2[s] == 1) {
      if (age_growth[i] <= tc[s, x]) {
        mu = mu1;
        sigma_i = sigma_juv;
      } else {
        mu = Lc[s, x] * exp((a2[s, x] / alpha2[s, x]) * (1.0 - exp(-alpha2[s, x] * (age_growth[i] - tc[s, x]))));
        sigma_i = sigma_adult;
      }
    } else {
      mu = mu1;
      sigma_i = sigma_juv;
    }

    // Stores the log likelihood for the 'loo' package
    log_lik[i] = normal_lpdf(length_growth[i] | mu, sigma_i);
  }
}
