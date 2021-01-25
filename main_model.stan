data {
  int<lower=0> N;                                       // Number of replicats
  int<lower=0> Survived[N];                             // Number of survived individuals
  int<lower=0> Total[N];                                // Total number of individuals
  int<lower=0> N_populations;                           // Number of populations
  int<lower=0,upper=N_populations> Population[N];       // Population ID
  int<lower=0> N_clones;                                // Number of clones
  int<lower=0,upper=N_clones> Clone[N];                 // Clone ID
  int<lower=0,upper=N_populations> Clonepop[N_clones];  // Population ID for each
  int<lower=0> N_X;                                     // Number of regression terms
  matrix[N,N_X] X;                                      // Design matrix
  int<lower=0> N_treatments;                            // Number of unique treatments
  matrix[N_treatments, N_X] X_treatments;               // Matrix of unique treatments
  real prior_beta[3];                                   // Parameters of the prior on the regression parameters
  real prior_sd_population[3];                          // Parameters of the prior on the population-level s.d.
  real prior_sd_clone[3];                               // Parameters of the prior on the clone-level s.d.
  real prior_nu[2];                                     // Parameters of the prior on the Student's t d.f.
}

parameters {
  vector[N_X] beta;                   // Regression coefficients
  vector[N_populations] z_population; // Standardised population random effect
  vector[N_clones] z_clone;           // Standardised clone random effect
  real<lower=0> sd_population;        // Population-level standard deviation
  real<lower=0> sd_clone;             // Clone-level standard deviation
  real<lower=1> nu;                   // Student's t degrees of freedom
}

transformed parameters {
  // Fixed linear predictor part
  vector[N] mu_fixed = X * beta;
  // Random linear predictor part
  vector[N] mu_random = z_population[Population]*sd_population + z_clone[Clone]*sd_clone;
  // Linear predictor
  vector[N] mu = mu_fixed + to_vector(X[,1]) .* mu_random;
}

model {
  // Priors
  target += student_t_lpdf(beta | prior_beta[1], prior_beta[2], prior_beta[3]);
  target += std_normal_lpdf(z_population);
  target += student_t_lpdf(z_clone | nu, 0, 1);
  target += gamma_lpdf(nu | prior_nu[1], prior_nu[2]);
  target += student_t_lpdf(sd_population | prior_sd_population[1], prior_sd_population[2], prior_sd_population[3]);
  target += student_t_lpdf(sd_clone | prior_sd_clone[1], prior_sd_clone[2], prior_sd_clone[3]);
  // Likelihood
  target += binomial_logit_lpmf(Survived | Total, mu);
}

generated quantities {
  // Posterior predictive distribution for the survival
  int<lower=0> Survived_rep[N] = binomial_rng(Total, inv_logit(mu));
  // Population effects
  vector[N_populations] Population_effects = z_population*sd_population;
  // Clone effects
  vector[N_clones] Clone_effects = Population_effects[Clonepop[1:N_clones]] + z_clone[1:N_clones]*sd_clone;
  // Treatment effects
  vector[N_treatments] Treatment_effects = inv_logit(X_treatments * beta);
}
