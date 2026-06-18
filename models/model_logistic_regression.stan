// Stan model to fit logistic curves to sexual maturity-age relationships
// of Stenella with n.sp species (attenuata [1], coeruleoalba [2], and
// longirostris [3]). 

data {
  int<lower=1> N;                               // Total number of fetuses + neonates
  int<lower=1> N_species;                       // Number of species
  
  vector[N] X;                             // Length or age of the individual
  array[N] int<lower=0, upper=1> y;        // 0 = inmature, 1 = mature
  array[N] int<lower=1, upper=N_species> species; // Species index
}

parameters {
  vector[N_species] beta0;          // Intercept for each species
  
  // Constrain beta1 to be positive: probability of being born MUST increase with length
  vector<lower=0>[N_species] beta1; 
}

model {
  // Priors
  beta0 ~ normal(0, 50);
  beta1 ~ normal(0, 10);
  
  // Likelihood
  for (i in 1:N) {
    y[i] ~ bernoulli_logit(beta0[species[i]] + beta1[species[i]] * X[i]);
  }
}

generated quantities {
  // Calculate L0 directly from the logistic regression coefficients
  vector[N_species] L0;
  
  for (s in 1:N_species) {
    L0[s] = -beta0[s] / beta1[s];
  }
}
