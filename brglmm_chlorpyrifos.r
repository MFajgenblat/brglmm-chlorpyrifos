################################################################################
# Packages
################################################################################

library(tidyverse)
library(rstan)
library(tidybayes)

################################################################################
# Setting seed and reading data
################################################################################

set.seed(1234)
survival_data <- read.csv("chlorpyrifos_data.csv", sep=";")

################################################################################
# Data exploration
################################################################################

# Visualization of the raw data per treatment
ggplot(survival_data) +
  geom_jitter(aes(x = Origin, y = Temperature, fill = Survived/Total), shape = 21, color = "grey75", size = 3) +
  facet_wrap(~ CPF) +
  scale_fill_distiller("Survival", palette = "YlOrRd", direction = 1, breaks = seq(0, 1, by = 0.2)) +
  scale_x_discrete(expand = c(0.25,0.25)) + scale_y_discrete(expand = c(0.25,0.25)) +
  coord_equal() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 9),
        axis.text = element_text(size = 9),
        axis.title = element_text(face = "bold", size = 9),
        legend.title = element_text(face = "bold", size = 9, vjust = 4),
        legend.text = element_text(size = 9))

# Visualizing of the average survival per treatment
survival_data %>%
  group_by(Origin, Temperature, CPF) %>%
  summarise(Survival = sum(Survived)/sum(Total)) %>%
  mutate(Origin = factor(Origin, levels = c("Rural", "Urban"))) %>%
  ggplot() +
  geom_bar(aes(x = CPF, y = Survival, fill = Origin), stat = "identity", position = "dodge") +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2), expand = c(0,0)) +
  scale_fill_manual("Temperature", values = c("#44a7bd", "#bd4460")) +
  facet_wrap(~ Temperature, strip.position = "bottom") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 9),
        strip.placement = "outside",
        axis.text = element_text(size = 9),
        axis.title.y = element_text(face = "bold", size = 9),
        axis.title.x = element_blank(),
        axis.line = element_line(color = "black"),
        legend.title = element_text(face = "bold", size = 9, vjust = 4),
        legend.text = element_text(size = 9),
        legend.key = element_blank())

# Visualizing clone-level survival
ggplot(survival_data) +
  geom_jitter(aes(x = Clone, y = Survived/Total, color = Temperature), height = 0.05) +
  scale_x_discrete(guide = guide_axis(n.dodge = 3)) +
  scale_y_continuous("Survival") +
  scale_color_manual(values = c("#17b6ff", "#ff5c17")) +
  facet_grid(CPF ~ Origin, scales = "free_x", switch = "x") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 7),
        legend.text = element_text(size = 7),
        legend.key = element_blank(),
        axis.title = element_text(face = "bold", size = 7),
        axis.text = element_text(size = 7),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 7),
        strip.placement = "outside")

################################################################################
# Model fitting
################################################################################

# Preparing the input data
data <- list(N = nrow(survival_data),
             Total = survival_data$Total,
             Survived = survival_data$Survived,
             N_populations = length(unique(survival_data$Population)),
             Population = as.numeric(as.factor(survival_data$Population)),
             N_clones = length(unique(survival_data$Clone)),
             Clone = as.numeric(as.factor(survival_data$Clone)),
             Clonepop = as.numeric(factor(substr(levels(as.factor(survival_data$Clone)), 1, 2),
                                          levels = levels(as.factor(survival_data$Population)))),
             N_X = ncol(cbind(1, model.matrix(~ CPF*Origin*Temperature, data = survival_data)[,-1] - 0.5)),
             X = cbind(1, model.matrix(~ CPF*Origin*Temperature, data = survival_data)[,-1] - 0.5),
             N_treatments = 8,
             X_treatments = unique(cbind(1, model.matrix(~ CPF*Origin*Temperature, data = survival_data)[,-1] - 0.5))[c(5,6,1,2,7,8,3,4),],
             prior_beta = c(3, 0, 5),
             prior_sd_population = c(3, 0, 5),
             prior_sd_clone = c(3, 0, 5),
             prior_nu = c(2, 0.1),
             Origin = as.numeric(factor(substr(levels(factor(survival_data$Clone)), 1, 1))))

# Using all available CPU cores and saving a compiled version
options(mc.cores = parallel::detectCores(), auto_write = TRUE)

# Choosing the model file name
# Uncomment the second line to use an extended model version that addresses overdispersion
modelname <- "main_model.stan"
# modelname <- "oVerdispersion_model.stan"

# Fitting the Stan model
fit <- stan(file = modelname, data = data, iter = 10000)

################################################################################
# Assessing convergence, prior sensitivity analysis & model criticism
################################################################################

# For clarity, these parts are moved to the bottom part of this script

################################################################################
# Model interpretation
################################################################################

# Posterior summaries for the main and interaction effects
fit %>%
  spread_draws(beta[effect]) %>%
  group_by(effect) %>%
  summarise(sign = sign(mean(beta)),
            effect = mean(beta),
            CI_lower = quantile(beta, 0.025),
            CI_upper = quantile(beta, 0.975),
            absolute_effect = abs(mean(beta)),
            probability = mean(sign*beta > 0))

# Posteriors densities for the main and interaction effects
fit %>%
  spread_draws(beta[effect]) %>%
  mutate(Highlight = case_when(effect %in% c(1,2,5) ~ "Yes",
                               effect %in% c(3,4,6,7,8) ~ "No")) %>%
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  stat_halfeye(aes(x = beta, y = effect, fill = Highlight),
               .width = c(0.8, 0.95), interval_size_range = c(0.3, 0.6),
               point_size = 0.6, orientation = "horizontal") +
  scale_fill_manual(values = c("grey", "#f0b181"), guide = F) +
  scale_x_continuous("Effect (logit scale)", breaks = seq(-20, 20, by = 2)) +
  scale_y_reverse(breaks = 1:8,
                  labels = c("Intercept", "Pesticide", "Urbanization",
                             "Temperature", "Pesticide \u00D7\nUrbanization",
                             "Pesticide \u00D7\nTemperature",
                             "Urbanization \u00D7\nTemperature",
                             "Pesticide \u00D7\nUrbanization \u00D7\nTemperature")) +
  theme(axis.line.x = element_line(color = "black"),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(face = "bold", size = 7, color = "black"),
        axis.title.x = element_text(face = "bold", size = 7),
        axis.title.y = element_blank(),
        panel.grid.major = element_line(color = "grey96"),
        plot.background = element_rect(fill = NA, color = NA), 
        panel.background = element_rect(fill = NA, color = NA),
        panel.border = element_blank())
ggsave("Figure_1.png", width = 11, height = 9, units = "cm", dpi = 600)

# Estimates for the 8 experimental conditions
fit %>%
  spread_draws(Treatment_effects[treatment]) %>%
  mutate(Origin = case_when(treatment %in% c(3,4,7,8) ~ "Rural",
                            treatment %in% c(1,2,5,6) ~ "Urban"),
         Origin = factor(Origin, levels = c("Rural", "Urban")),
         Temperature  = case_when(treatment %in% c(1,3,5,7) ~ "20 °C",
                                  treatment %in% c(2,4,6,8) ~ "24 °C"),
         CPF = case_when(treatment %in% c(1,2,3,4) ~ "Control",
                         treatment %in% c(5,6,7,8) ~ "Pesticide")) %>%
  ggplot() +
  stat_summary(aes(y = Treatment_effects, x = CPF, fill = Origin), fun = "median", geom = "bar", position = "dodge") +
  stat_pointinterval(aes(x = CPF, y = Treatment_effects, color = Origin), position = "dodge", .width = c(0.8, 0.95),
                     interval_size_range = c(0.3, 0.6), point_size = 0.3) +
  scale_y_continuous("Estimated survival probability", limits = c(0,1.02), breaks = seq(0,1,0.2), expand = c(0,0)) +
  scale_color_manual(values = c("black", "black"), guide = F) +
  scale_fill_manual(values = c("#44a7bd", "#bd4460")) +
  facet_wrap(~ Temperature, strip.position = "bottom") +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(color = "grey90"),
        strip.background = element_blank(),
        strip.text = element_text(size = 7),
        strip.placement = "outside",
        axis.text = element_text(size = 7, color = "black"),
        axis.title.y = element_text(face = "bold", size = 7),
        axis.title.x = element_blank(),
        axis.line = element_line(color = "black", size = 0.3),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.key = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.position = "bottom",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-8,0,0,0))
ggsave("Figure_2.png", width = 8, height = 6.5, units = "cm", dpi = 600)

# Clone-specific effects
fit %>%
  spread_draws(Clone_effects[Clone]) %>%
  merge(data.frame(Clone = 1:49,
                   Clonepop = factor(substr(levels(as.factor(survival_data$Clone)), 1, 2),
                                     levels = levels(as.factor(survival_data$Population))),
                   Clonename = factor(substr(levels(as.factor(survival_data$Clone)), 4, 4), levels = 5:1)))  %>%
  mutate(Clone = factor(Clone, levels = 49:1)) %>%
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  stat_interval(aes(x = Clone_effects, y = Clonename), .width = c(0.5, 0.8, 0.95, 0.99), size = 1.4) +
  scale_color_brewer("Credible\nintervals") +
  scale_x_continuous("Estimated effect (logit scale)", expand = c(0,0), seq(-10, 10, by = 2)) +
  scale_y_discrete("Population and clone") +
  facet_wrap(~ Clonepop, ncol = 1, scales = "free_y", strip.position = "left") +
  coord_cartesian(xlim = c(-7,7)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        panel.grid = element_line(colour = "grey95"),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = 7, face = "bold"),
        axis.text = element_text(size = 7),
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 7),
        legend.key = element_blank(),
        legend.justification = "top",
        legend.key.width = unit(0.4, "cm"),
        strip.background = element_blank(),
        strip.text.y.left = element_text(size = 7, face = "bold", angle = 0),
        strip.placement = "outside")
ggsave("Figure_clone_effects.png", width = 8, height = 18, units = "cm", dpi = 600)

# Computing specific posterior probabilities for the results section
expit <- function(x) {return(exp(x)/(1+exp(x)))}
## Effect of chlorpyrifos exposure
median(as.matrix(fit, pars = "beta[2]"))
quantile(as.matrix(fit, pars = "beta[2]"), c(0.025, 0.975))
## Marginal survival probability without exposure
median(expit(as.matrix(fit, pars = "beta[1]") - 0.5*as.matrix(fit, pars = "beta[2]")))
quantile(expit(as.matrix(fit, pars = "beta[1]") - 0.5*as.matrix(fit, pars = "beta[2]")), c(0.025, 0.975))
## Marginal survival probability with exposure
median(expit(as.matrix(fit, pars = "beta[1]") + 0.5*as.matrix(fit, pars = "beta[2]")))
quantile(expit(as.matrix(fit, pars = "beta[1]") + 0.5*as.matrix(fit, pars = "beta[2]")), c(0.025, 0.975))
## Reduction in marginal survival probability upon exposure
median(expit(as.matrix(fit, pars = "beta[1]") - 0.5*as.matrix(fit, pars = "beta[2]"))-expit(as.matrix(fit, pars = "beta[1]") + 0.5*as.matrix(fit, pars = "beta[2]")))
quantile(expit(as.matrix(fit, pars = "beta[1]") - 0.5*as.matrix(fit, pars = "beta[2]"))-expit(as.matrix(fit, pars = "beta[1]") + 0.5*as.matrix(fit, pars = "beta[2]")), c(0.025, 0.975))
## Effect of population origin
median(as.matrix(fit, pars = "beta[3]"))
quantile(as.matrix(fit, pars = "beta[3]"), c(0.025, 0.975))
## Effect of temperature
median(as.matrix(fit, pars = "beta[4]"))
quantile(as.matrix(fit, pars = "beta[4]"), c(0.025, 0.975))
## Reduction in marginal survival probability for urban population origin
median(expit(as.matrix(fit, pars = "beta[1]") - 0.5*as.matrix(fit, pars = "beta[3]"))-expit(as.matrix(fit, pars = "beta[1]") + 0.5*as.matrix(fit, pars = "beta[3]")))
quantile(expit(as.matrix(fit, pars = "beta[1]") - 0.5*as.matrix(fit, pars = "beta[3]"))-expit(as.matrix(fit, pars = "beta[1]") + 0.5*as.matrix(fit, pars = "beta[3]")), c(0.025, 0.975))
## Reduction in marginal survival probability for 24C temperature
median(expit(as.matrix(fit, pars = "beta[1]") - 0.5*as.matrix(fit, pars = "beta[4]"))-expit(as.matrix(fit, pars = "beta[1]") + 0.5*as.matrix(fit, pars = "beta[4]")))
quantile(expit(as.matrix(fit, pars = "beta[1]") - 0.5*as.matrix(fit, pars = "beta[4]"))-expit(as.matrix(fit, pars = "beta[1]") + 0.5*as.matrix(fit, pars = "beta[4]")), c(0.025, 0.975))
## Probability of a positive interaction effect between chlorpyrifos exposure and population origin
mean(as.matrix(fit, pars = "beta[5]")>0)
## Interaction effect of chlorpyrifos exposure and population origin
median(as.matrix(fit, pars = "beta[5]"))
quantile(as.matrix(fit, pars = "beta[5]"), c(0.025, 0.975))
## Marginal survival probability with exposure and urban origin
median(expit(as.matrix(fit, pars = "beta[1]") + 0.5*as.matrix(fit, pars = "beta[2]") + 0.5*as.matrix(fit, pars = "beta[3]") + 0.5*as.matrix(fit, pars = "beta[5]")))
quantile(expit(as.matrix(fit, pars = "beta[1]") + 0.5*as.matrix(fit, pars = "beta[2]") + 0.5*as.matrix(fit, pars = "beta[3]") + 0.5*as.matrix(fit, pars = "beta[5]")), c(0.025, 0.975))
## Marginal survival probability with exposure and rural origin
median(expit(as.matrix(fit, pars = "beta[1]") + 0.5*as.matrix(fit, pars = "beta[2]") - 0.5*as.matrix(fit, pars = "beta[3]") - 0.5*as.matrix(fit, pars = "beta[5]")))
quantile(expit(as.matrix(fit, pars = "beta[1]") + 0.5*as.matrix(fit, pars = "beta[2]") - 0.5*as.matrix(fit, pars = "beta[3]") - 0.5*as.matrix(fit, pars = "beta[5]")), c(0.025, 0.975))
## Interaction effect of chlorpyrifos exposure and temperature
median(as.matrix(fit, pars = "beta[6]"))
quantile(as.matrix(fit, pars = "beta[6]"), c(0.025, 0.975))
## Interaction effect of population origin and temperature
median(as.matrix(fit, pars = "beta[7]"))
quantile(as.matrix(fit, pars = "beta[7]"), c(0.025, 0.975))
## Probability of a negative interaction effect of chlorpyrifos exposure, population origin and temperature
mean(as.matrix(fit, pars = "beta[8]")<0)
## Interaction effect of chlorpyrifos exposure, population origin and temperature
median(as.matrix(fit, pars = "beta[8]"))
quantile(as.matrix(fit, pars = "beta[8]"), c(0.025, 0.975))
## Population-level standard deviation
median(as.matrix(fit, pars = "sd_population"))
quantile(as.matrix(fit, pars = "sd_population"), c(0.025, 0.975))
## Clone-level standard deviation
median(as.matrix(fit, pars = "sd_clone"))
quantile(as.matrix(fit, pars = "sd_clone"), c(0.025, 0.975))
## Number of degrees of freedom for the Student's t-distributed clone-level random effects
median(as.matrix(fit, pars = "nu"))
quantile(as.matrix(fit, pars = "nu"), c(0.025, 0.975))

################################################################################
# Assessing model convergence
################################################################################

# Distribution of Rhat values (Potential Scale Reduction Factors)
hist(summary(fit)$summary[,"Rhat"])

# Highest Rhat value (Potential Scale Reduction Factors)
max(summary(fit)$summary[,"Rhat"])

# Distribution of effective sample sizes
hist(summary(fit)$summary[,"n_eff"])

# Lowest number of effective sample sizes
min(summary(fit)$summary[,"n_eff"], na.rm = T)

# Traceplots

## Regression parameters
fit %>%
  spread_draws(beta[treatment]) %>%
  mutate(treatment = c("Intercept", "Pesticide", "Urbanization",
                       "Temperature", "Pesticide \u00D7 Urbanization",
                       "Pesticide \u00D7 Temperature",
                       "Urbanization \u00D7 Temperature",
                       "Pesticide \u00D7 Urbanization \u00D7 Temperature")[treatment]) %>%
  ggplot() +
  geom_line(aes(x = .iteration, y = beta, color = factor(.chain)), alpha = 0.5, size = 0.3) +
  scale_color_brewer(palette = "Dark2", guide = F) +
  scale_x_continuous("Iteration", expand = c(0,0), limits = c(-99,5100)) +
  scale_y_continuous("Estimate") +
  facet_wrap(~ treatment, scales = "free_y", ncol = 2) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "grey25", fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 7),
        axis.title = element_text(face = "bold", size = 7),
        axis.text = element_text(size = 7))
ggsave("Traceplots_regression_parameters.png", width = 16, height = 9, units = "cm", dpi = 600)

## Random effect s.d. & number of degrees of freedom
fit %>%
  spread_draws(sd_population, sd_clone, nu) %>%
  pivot_longer(c(sd_population, sd_clone, nu), names_to = "variable", values_to = "value") %>%
  mutate(variable = case_when(variable == "sd_population" ~ "Population-level s.d.",
                              variable == "sd_clone" ~ "Clone-level s.d.",
                              variable == "nu" ~ "Number of degrees of freedom"),
         variable = factor(variable, levels = c("Population-level s.d.", "Clone-level s.d.", "Number of degrees of freedom"))) %>%
  ggplot() +
  geom_line(aes(x = .iteration, y = value, color = factor(.chain)), alpha = 0.5, size = 0.3) +
  scale_color_brewer(palette = "Dark2", guide = F) +
  scale_x_continuous("Iteration", expand = c(0,0), limits = c(-99,5100)) +
  scale_y_continuous("Estimate") +
  facet_wrap(~ variable, scales = "free_y", ncol = 2) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "grey25", fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 7),
        axis.title = element_text(face = "bold", size = 7),
        axis.text = element_text(size = 7))
ggsave("Traceplots_sd_nu.png", width = 16, height = 5, units = "cm", dpi = 600)

## Random subset of standardized population effects
randompopulations <- sort(sample(1:10, 8, replace = F))
fit %>%
  spread_draws(z_population[population]) %>%
  filter(population %in% randompopulations) %>%
  mutate(population = paste("Population", levels(factor(survival_data$Population))[population])) %>%
  ggplot() +
  geom_line(aes(x = .iteration, y = z_population, color = factor(.chain)), alpha = 0.5, size = 0.3) +
  scale_color_brewer(palette = "Dark2", guide = F) +
  scale_x_continuous("Iteration", expand = c(0,0), limits = c(-99,5100)) +
  scale_y_continuous("Estimate") +
  facet_wrap(~ population, scales = "free_y", ncol = 2) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "grey25", fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 7),
        axis.title = element_text(face = "bold", size = 7),
        axis.text = element_text(size = 7))
ggsave("Traceplots_population_effects.png", width = 16, height = 9, units = "cm", dpi = 600)

## Random subset of standardized clone effects
randomclones <- sort(sample(1:49, 8, replace = F))
fit %>%
  spread_draws(z_clone[clone]) %>%
  filter(clone %in% randomclones) %>%
  mutate(clone = paste("Clone", levels(factor(survival_data$Clone))[clone])) %>%
  ggplot() +
  geom_line(aes(x = .iteration, y = z_clone, color = factor(.chain)), alpha = 0.5, size = 0.3) +
  scale_color_brewer(palette = "Dark2", guide = F) +
  scale_x_continuous("Iteration", expand = c(0,0), limits = c(-99,5100)) +
  scale_y_continuous("Estimate") +
  facet_wrap(~ clone, scales = "free_y", ncol = 2) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "grey25", fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 7),
        axis.title = element_text(face = "bold", size = 7),
        axis.text = element_text(size = 7))
ggsave("Traceplots_clone_effects.png", width = 16, height = 9, units = "cm", dpi = 600)

##  Log probability (up to an additive constant)
fit %>%
  spread_draws(lp__) %>%
  ggplot() +
  geom_line(aes(x = .iteration, y = lp__, color = factor(.chain)), alpha = 0.5, size = 0.3) +
  scale_color_brewer(palette = "Dark2", guide = F) +
  scale_x_continuous("Iteration", expand = c(0,0), limits = c(-99,5100)) +
  scale_y_continuous("Estimate") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "grey25", fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 7),
        axis.title = element_text(face = "bold", size = 7),
        axis.text = element_text(size = 7))
ggsave("Traceplots_lp.png", width = 16, height = 4, units = "cm", dpi = 600)

################################################################################
# Prior sensitivity analysis
################################################################################

# Scenario A1: alternative prior choice for the regression parameters only
data$prior_beta = c(3, 0, 1); data$prior_sd_population = c(3, 0, 5); data$prior_sd_clone = c(3, 0, 5); data$prior_nu = c(2, 0.1)
fit_a <- stan(file = "CPF_binomial_robust3.stan", data = data, iter = 2000, pars = c("beta", "sd_population", "sd_clone", "nu"))

data$prior_beta = c(3, 0, 2.5); data$prior_sd_population = c(3, 0, 5); data$prior_sd_clone = c(3, 0, 5); data$prior_nu = c(2, 0.1)
fit_b <- stan(file = "CPF_binomial_robust3.stan", data = data, iter = 2000, pars = c("beta", "sd_population", "sd_clone", "nu"))

data$prior_beta = c(3, 0, 5); data$prior_sd_population = c(3, 0, 5); data$prior_sd_clone = c(3, 0, 5); data$prior_nu = c(2, 0.1)
fit_c <- stan(file = "CPF_binomial_robust3.stan", data = data, iter = 2000, pars = c("beta", "sd_population", "sd_clone", "nu"))

data$prior_beta = c(3, 0, 10); data$prior_sd_population = c(3, 0, 5); data$prior_sd_clone = c(3, 0, 5); data$prior_nu = c(2, 0.1)
fit_d <- stan(file = "CPF_binomial_robust3.stan", data = data, iter = 2000, pars = c("beta", "sd_population", "sd_clone", "nu"))

data$prior_beta = c(3, 0, 25); data$prior_sd_population = c(3, 0, 5); data$prior_sd_clone = c(3, 0, 5); data$prior_nu = c(2, 0.1)
fit_e <- stan(file = "CPF_binomial_robust3.stan", data = data, iter = 2000, pars = c("beta", "sd_population", "sd_clone", "nu"))

pivot_longer(rbind(data.frame(Priors = "Strongly tighter", as.matrix(fit_a)), data.frame(Priors = "Tighter", as.matrix(fit_b)), data.frame(Priors = "Original", as.matrix(fit_c)), data.frame(Priors = "Wider", as.matrix(fit_d)), data.frame(Priors = "Strongly wider", as.matrix(fit_e))), !Priors) %>%
  mutate(Priors = factor(Priors, levels = c("Strongly tighter", "Tighter", "Original", "Wider", "Strongly wider")),
         name = factor(c("Intercept", "Pesticide", "Urbanization",
                         "Temperature", "Pesticide \u00D7\nUrbanization",
                         "Pesticide \u00D7\nTemperature",
                         "Urbanization \u00D7\nTemperature",
                         "Pesticide \u00D7\nUrbanization \u00D7\nTemperature",
                         "log probability", "degrees of freedom",
                         "clone s.d.", "population s.d.")[as.numeric(as.factor(name))],
                       levels = c("Intercept", "Pesticide", "Urbanization",
                                  "Temperature", "Pesticide \u00D7\nUrbanization",
                                  "Pesticide \u00D7\nTemperature",
                                  "Urbanization \u00D7\nTemperature",
                                  "Pesticide \u00D7\nUrbanization \u00D7\nTemperature",
                                  "log probability", "clone s.d.", 
                                  "population s.d.", "degrees of freedom"))) %>%
  ggplot() +
  geom_density(aes(x = value, fill = Priors, color = Priors), alpha = 0.2) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  scale_y_continuous("Density", expand = c(0,0)) +
  scale_color_discrete("Prior choice") +
  scale_fill_discrete("Prior choice") +
  facet_wrap(~ name, ncol = 3, scales = "free", strip.position = "bottom") +
  theme(panel.background = element_blank(),
        panel.grid = element_line(color = "grey95"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 7),
        legend.text = element_text(size = 7),
        legend.key = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", color = "black", size = 7),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 7),
        axis.text = element_text(size = 7),
        axis.line.x = element_line(color = "black"))
ggsave("Prior_sensitivity_analysis_A1.png", width = 16, height = 20, units = "cm", dpi = 600)

# Scenario A2: alternative prior choice for the population-level sd only
data$prior_beta = c(3, 0, 5); data$prior_sd_population = c(3, 0, 1); data$prior_sd_clone = c(3, 0, 5); data$prior_nu = c(2, 0.1)
fit_a <- stan(file = "CPF_binomial_robust3.stan", data = data, iter = 2000, pars = c("beta", "sd_population", "sd_clone", "nu"))

data$prior_beta = c(3, 0, 5); data$prior_sd_population = c(3, 0, 2.5); data$prior_sd_clone = c(3, 0, 5); data$prior_nu = c(2, 0.1)
fit_b <- stan(file = "CPF_binomial_robust3.stan", data = data, iter = 2000, pars = c("beta", "sd_population", "sd_clone", "nu"))

data$prior_beta = c(3, 0, 5); data$prior_sd_population = c(3, 0, 5); data$prior_sd_clone = c(3, 0, 5); data$prior_nu = c(2, 0.1)
fit_c <- stan(file = "CPF_binomial_robust3.stan", data = data, iter = 2000, pars = c("beta", "sd_population", "sd_clone", "nu"))

data$prior_beta = c(3, 0, 5); data$prior_sd_population = c(3, 0, 10); data$prior_sd_clone = c(3, 0, 5); data$prior_nu = c(2, 0.1)
fit_d <- stan(file = "CPF_binomial_robust3.stan", data = data, iter = 2000, pars = c("beta", "sd_population", "sd_clone", "nu"))

data$prior_beta = c(3, 0, 5); data$prior_sd_population = c(3, 0, 25); data$prior_sd_clone = c(3, 0, 5); data$prior_nu = c(2, 0.1)
fit_e <- stan(file = "CPF_binomial_robust3.stan", data = data, iter = 2000, pars = c("beta", "sd_population", "sd_clone", "nu"))

pivot_longer(rbind(data.frame(Priors = "Strongly tighter", as.matrix(fit_a)), data.frame(Priors = "Tighter", as.matrix(fit_b)), data.frame(Priors = "Original", as.matrix(fit_c)), data.frame(Priors = "Wider", as.matrix(fit_d)), data.frame(Priors = "Strongly wider", as.matrix(fit_e))), !Priors) %>%
  mutate(Priors = factor(Priors, levels = c("Strongly tighter", "Tighter", "Original", "Wider", "Strongly wider")),
         name = factor(c("Intercept", "Pesticide", "Urbanization",
                         "Temperature", "Pesticide \u00D7\nUrbanization",
                         "Pesticide \u00D7\nTemperature",
                         "Urbanization \u00D7\nTemperature",
                         "Pesticide \u00D7\nUrbanization \u00D7\nTemperature",
                         "log probability", "degrees of freedom",
                         "clone s.d.", "population s.d.")[as.numeric(as.factor(name))],
                       levels = c("Intercept", "Pesticide", "Urbanization",
                                  "Temperature", "Pesticide \u00D7\nUrbanization",
                                  "Pesticide \u00D7\nTemperature",
                                  "Urbanization \u00D7\nTemperature",
                                  "Pesticide \u00D7\nUrbanization \u00D7\nTemperature",
                                  "log probability", "clone s.d.", 
                                  "population s.d.", "degrees of freedom"))) %>%
  ggplot() +
  geom_density(aes(x = value, fill = Priors, color = Priors), alpha = 0.2) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  scale_y_continuous("Density", expand = c(0,0)) +
  scale_color_discrete("Prior choice") +
  scale_fill_discrete("Prior choice") +
  facet_wrap(~ name, ncol = 3, scales = "free", strip.position = "bottom") +
  theme(panel.background = element_blank(),
        panel.grid = element_line(color = "grey95"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 7),
        legend.text = element_text(size = 7),
        legend.key = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", color = "black", size = 7),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 7),
        axis.text = element_text(size = 7),
        axis.line.x = element_line(color = "black"))
ggsave("Prior_sensitivity_analysis_A2.png", width = 16, height = 20, units = "cm", dpi = 600)

# Scenario A3: alternative prior choice for the clone-level sd only
data$prior_beta = c(3, 0, 5); data$prior_sd_population = c(3, 0, 5); data$prior_sd_clone = c(3, 0, 1); data$prior_nu = c(2, 0.1)
fit_a <- stan(file = "CPF_binomial_robust3.stan", data = data, iter = 2000, pars = c("beta", "sd_population", "sd_clone", "nu"))

data$prior_beta = c(3, 0, 5); data$prior_sd_population = c(3, 0, 5); data$prior_sd_clone = c(3, 0, 2.5); data$prior_nu = c(2, 0.1)
fit_b <- stan(file = "CPF_binomial_robust3.stan", data = data, iter = 2000, pars = c("beta", "sd_population", "sd_clone", "nu"))

data$prior_beta = c(3, 0, 5); data$prior_sd_population = c(3, 0, 5); data$prior_sd_clone = c(3, 0, 5); data$prior_nu = c(2, 0.1)
fit_c <- stan(file = "CPF_binomial_robust3.stan", data = data, iter = 2000, pars = c("beta", "sd_population", "sd_clone", "nu"))

data$prior_beta = c(3, 0, 5); data$prior_sd_population = c(3, 0, 5); data$prior_sd_clone = c(3, 0, 10); data$prior_nu = c(2, 0.1)
fit_d <- stan(file = "CPF_binomial_robust3.stan", data = data, iter = 2000, pars = c("beta", "sd_population", "sd_clone", "nu"))

data$prior_beta = c(3, 0, 5); data$prior_sd_population = c(3, 0, 5); data$prior_sd_clone = c(3, 0, 25); data$prior_nu = c(2, 0.1)
fit_e <- stan(file = "CPF_binomial_robust3.stan", data = data, iter = 2000, pars = c("beta", "sd_population", "sd_clone", "nu"))

pivot_longer(rbind(data.frame(Priors = "Strongly tighter", as.matrix(fit_a)), data.frame(Priors = "Tighter", as.matrix(fit_b)), data.frame(Priors = "Original", as.matrix(fit_c)), data.frame(Priors = "Wider", as.matrix(fit_d)), data.frame(Priors = "Strongly wider", as.matrix(fit_e))), !Priors) %>%
  mutate(Priors = factor(Priors, levels = c("Strongly tighter", "Tighter", "Original", "Wider", "Strongly wider")),
         name = factor(c("Intercept", "Pesticide", "Urbanization",
                         "Temperature", "Pesticide \u00D7\nUrbanization",
                         "Pesticide \u00D7\nTemperature",
                         "Urbanization \u00D7\nTemperature",
                         "Pesticide \u00D7\nUrbanization \u00D7\nTemperature",
                         "log probability", "degrees of freedom",
                         "clone s.d.", "population s.d.")[as.numeric(as.factor(name))],
                       levels = c("Intercept", "Pesticide", "Urbanization",
                                  "Temperature", "Pesticide \u00D7\nUrbanization",
                                  "Pesticide \u00D7\nTemperature",
                                  "Urbanization \u00D7\nTemperature",
                                  "Pesticide \u00D7\nUrbanization \u00D7\nTemperature",
                                  "log probability", "clone s.d.", 
                                  "population s.d.", "degrees of freedom"))) %>%
  ggplot() +
  geom_density(aes(x = value, fill = Priors, color = Priors), alpha = 0.2) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  scale_y_continuous("Density", expand = c(0,0)) +
  scale_color_discrete("Prior choice") +
  scale_fill_discrete("Prior choice") +
  facet_wrap(~ name, ncol = 3, scales = "free", strip.position = "bottom") +
  theme(panel.background = element_blank(),
        panel.grid = element_line(color = "grey95"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 7),
        legend.text = element_text(size = 7),
        legend.key = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", color = "black", size = 7),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 7),
        axis.text = element_text(size = 7),
        axis.line.x = element_line(color = "black"))
ggsave("Prior_sensitivity_analysis_A3.png", width = 16, height = 20, units = "cm", dpi = 600)

# Scenario B: alternative prior choice for df only
data$prior_beta = c(3, 0, 5); data$prior_sd_population = c(3, 0, 5); data$prior_sd_clone = c(3, 0, 5); data$prior_nu = c(2, 0.5)
fit_a <- stan(file = "CPF_binomial_robust3.stan", data = data, iter = 2000, pars = c("beta", "sd_population", "sd_clone", "nu"))

data$prior_beta = c(3, 0, 5); data$prior_sd_population = c(3, 0, 5); data$prior_sd_clone = c(3, 0, 5); data$prior_nu = c(2, 0.1)
fit_b <- stan(file = "CPF_binomial_robust3.stan", data = data, iter = 2000, pars = c("beta", "sd_population", "sd_clone", "nu"))

data$prior_beta = c(3, 0, 5); data$prior_sd_population = c(3, 0, 5); data$prior_sd_clone = c(3, 0, 5); data$prior_nu = c(3, 0.05)
fit_c <- stan(file = "CPF_binomial_robust3.stan", data = data, iter = 2000, pars = c("beta", "sd_population", "sd_clone", "nu"))

pivot_longer(rbind(data.frame(Priors = "Tighter", as.matrix(fit_a)), data.frame(Priors = "Original", as.matrix(fit_b)), data.frame(Priors = "Wider", as.matrix(fit_c))), !Priors) %>%
  mutate(Priors = factor(Priors, levels = c("Strongly tighter", "Tighter", "Original", "Wider", "Strongly wider")),
         name = factor(c("Intercept", "Pesticide", "Urbanization",
                         "Temperature", "Pesticide \u00D7\nUrbanization",
                         "Pesticide \u00D7\nTemperature",
                         "Urbanization \u00D7\nTemperature",
                         "Pesticide \u00D7\nUrbanization \u00D7\nTemperature",
                         "log probability", "degrees of freedom",
                         "clone s.d.", "population s.d.")[as.numeric(as.factor(name))],
                       levels = c("Intercept", "Pesticide", "Urbanization",
                                  "Temperature", "Pesticide \u00D7\nUrbanization",
                                  "Pesticide \u00D7\nTemperature",
                                  "Urbanization \u00D7\nTemperature",
                                  "Pesticide \u00D7\nUrbanization \u00D7\nTemperature",
                                  "log probability", "clone s.d.", 
                                  "population s.d.", "degrees of freedom"))) %>%
  ggplot() +
  geom_density(aes(x = value, fill = Priors, color = Priors), alpha = 0.2) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  scale_y_continuous("Density", expand = c(0,0)) +
  scale_color_discrete("Prior choice") +
  scale_fill_discrete("Prior choice") +
  facet_wrap(~ name, ncol = 3, scales = "free", strip.position = "bottom") +
  theme(panel.background = element_blank(),
        panel.grid = element_line(color = "grey95"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 7),
        legend.text = element_text(size = 7),
        legend.key = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", color = "black", size = 7),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 7),
        axis.text = element_text(size = 7),
        axis.line.x = element_line(color = "black"))
ggsave("Prior_sensitivity_analysis_B.png", width = 16, height = 20, units = "cm", dpi = 600)

# Scenario C: alternative prior choices for all main parameters simultaneously
data$prior_beta = c(3, 0, 1); data$prior_sd_population = c(3, 0, 1); data$prior_sd_clone = c(3, 0, 1); data$prior_nu = c(2, 0.1)
fit_a <- stan(file = "CPF_binomial_robust3.stan", data = data, iter = 2000, pars = c("beta", "sd_population", "sd_clone", "nu"))

data$prior_beta = c(3, 0, 2.5); data$prior_sd_population = c(3, 0, 2.5); data$prior_sd_clone = c(3, 0, 2.5); data$prior_nu = c(2, 0.1)
fit_b <- stan(file = "CPF_binomial_robust3.stan", data = data, iter = 2000, pars = c("beta", "sd_population", "sd_clone", "nu"))

data$prior_beta = c(3, 0, 5); data$prior_sd_population = c(3, 0, 5); data$prior_sd_clone = c(3, 0, 5); data$prior_nu = c(2, 0.1)
fit_c <- stan(file = "CPF_binomial_robust3.stan", data = data, iter = 2000, pars = c("beta", "sd_population", "sd_clone", "nu"))

data$prior_beta = c(3, 0, 10); data$prior_sd_population = c(3, 0, 10); data$prior_sd_clone = c(3, 0, 10); data$prior_nu = c(2, 0.1)
fit_d <- stan(file = "CPF_binomial_robust3.stan", data = data, iter = 2000, pars = c("beta", "sd_population", "sd_clone", "nu"))

data$prior_beta = c(3, 0, 25); data$prior_sd_population = c(3, 0, 25); data$prior_sd_clone = c(3, 0, 25); data$prior_nu = c(2, 0.1)
fit_e <- stan(file = "CPF_binomial_robust3.stan", data = data, iter = 2000, pars = c("beta", "sd_population", "sd_clone", "nu"))

pivot_longer(rbind(data.frame(Priors = "Strongly tighter", as.matrix(fit_a)), data.frame(Priors = "Tighter", as.matrix(fit_b)), data.frame(Priors = "Original", as.matrix(fit_c)), data.frame(Priors = "Wider", as.matrix(fit_d)), data.frame(Priors = "Strongly wider", as.matrix(fit_e))), !Priors) %>%
  mutate(Priors = factor(Priors, levels = c("Strongly tighter", "Tighter", "Original", "Wider", "Strongly wider")),
         name = factor(c("Intercept", "Pesticide", "Urbanization",
                         "Temperature", "Pesticide \u00D7\nUrbanization",
                         "Pesticide \u00D7\nTemperature",
                         "Urbanization \u00D7\nTemperature",
                         "Pesticide \u00D7\nUrbanization \u00D7\nTemperature",
                         "log probability", "degrees of freedom",
                         "clone s.d.", "population s.d.")[as.numeric(as.factor(name))],
                       levels = c("Intercept", "Pesticide", "Urbanization",
                                  "Temperature", "Pesticide \u00D7\nUrbanization",
                                  "Pesticide \u00D7\nTemperature",
                                  "Urbanization \u00D7\nTemperature",
                                  "Pesticide \u00D7\nUrbanization \u00D7\nTemperature",
                                  "log probability", "clone s.d.", 
                                  "population s.d.", "degrees of freedom"))) %>%
  ggplot() +
  geom_density(aes(x = value, fill = Priors, color = Priors), alpha = 0.2) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  scale_y_continuous("Density", expand = c(0,0)) +
  scale_color_discrete("Prior choice") +
  scale_fill_discrete("Prior choice") +
  facet_wrap(~ name, ncol = 3, scales = "free", strip.position = "bottom") +
  theme(panel.background = element_blank(),
        panel.grid = element_line(color = "grey95"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 7),
        legend.text = element_text(size = 7),
        legend.key = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", color = "black", size = 7),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 7),
        axis.text = element_text(size = 7),
        axis.line.x = element_line(color = "black"))
ggsave("Prior_sensitivity_analysis_C.png", width = 16, height = 20, units = "cm", dpi = 600)

################################################################################
# Posterior predictive checks
################################################################################

survived_rep <- as.matrix(fit, pars = "Survived_rep")
survtable <- do.call(rbind, lapply(1:nrow(survived_rep), function(i) table(factor(survived_rep[i,], levels = 0:6))))
ggplot() +
  geom_bar(data = survival_data, aes(x = Survived, y = ..count.., fill = "1")) +
  stat_interval(data = mutate(pivot_longer(data.frame(survtable), cols = everything()), name = as.numeric(substr(name, 2, 2))), aes(x = name, y = value), .width = seq(0.5, 1, by = 0.1)) +
  scale_color_brewer(bquote(bold(y[rep]))) +
  scale_fill_manual("y", values = c("grey50"), labels = c("")) +
  scale_y_continuous("Number of replicates", expand = c(0,0), breaks = seq(0, 250, by = 25)) +
  scale_x_continuous("Number of surviving individuals", breaks = 0:6) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 7, color = "black"),
        axis.title = element_text(face = "bold", size = 7),
        axis.line = element_line(color = "black", size = 0.3),
        legend.title = element_text(face = "bold", size = 7),
        legend.text = element_text(size = 7),
        legend.key = element_blank(),
        legend.key.size = unit(0.3, "cm"))
ggsave("PPC.png", width = 12, height = 6, units = "cm", dpi = 600)
