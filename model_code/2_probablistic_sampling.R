###############################################
# TITLE: Probabilistic Sampling for B360S Model
# AUTHOR: Nichola Naylor (OI Pharma Partners Ltd), aided by GPT-4o,GTP-5 & Github co-pilot
# DATE: September 2025
#
# DESCRIPTION:
# This R script performs probabilistic sensitivity analysis (PSA) sampling
# for a stroke decision model. It:
# - Loads pre-processed parameters
# - Samples utilities, costs, and mortality risks by mRS score
# - Samples treatment probabilities using odds ratios and control probabilities
# - Samples other model parameters from gamma or uniform distributions
#
# INPUTS:
# - inputs/created_inputs/parameters_edited.RData (created by "1_hospital_data_processing.R")
#
# OUTPUTS:
# - inputs/created_inputs/mrs_samples.RData
# - inputs/created_inputs/mrs_samples_mean.RData
# - inputs/created_inputs/params.psa.sample.RData
# - inputs/created_inputs/parameters_edited.RData
# - inputs/created_inputs/parameters_post_psa.csv 
###############################################

#####   1. INITIALISE & LOAD PACKAGES   ####
# clear environment
rm(list=ls())

# Load necessary packages
library(conflicted)
library(truncnorm)
library(tidyverse)
library(data.table)
library(assertthat)
library(stringr)

# Resolve potential function conflicts (tidyverse vs data.table vs base)
conflict_prefer("merge", "data.table")
conflict_prefer("filter", "dplyr")

set.seed(456) # Set seed for reproducibility

########  2. LOAD AND PROCESS INPUTS   ####

# Read in the inputs
load("inputs/created_inputs/parameters_edited.RData")

# Coerce types once (prevents silent mismatches later)
parameters[, `:=`(
  base_case = as.numeric(base_case),
  PSA_low   = as.numeric(PSA_low),
  PSA_high  = as.numeric(PSA_high),
  mrs       = as.integer(mrs)
)]

inputs <- parameters[!is.na(mrs)]

# Separate the data into distribution, utility, and cost
dist <- inputs[model_param != "utility.mrs" & model_param != "cost.mrs" &
                 model_param!="rr.mort"] #distribution of patients across mrs scores
utility <- inputs[model_param == "utility.mrs" & mrs %in% c(0:5)] ## note mrs 6 is dead
cost <- inputs[model_param == "cost.mrs"& mrs %in% c(0:5)] 
mort <- inputs[model_param=="rr.mort" & mrs %in% c(2:5)] 

n.sample <- 1000 # Number of samples for PSA

#' Title: Check Distributions sum
#' Check that each distribution adds to 1 for each model_param group
#' @param data A dataframe with columns 'model_param'(character), 'base_case' (numeric)
#'
#' @return A vector of sums for each group
#' @export
check_distributions <- function(data) {
  # Split the data by the 'model_param' column
  distribution_groups <- split(data, data$model_param)
  
  # Check sums for each group
  sapply(distribution_groups, function(group) {
    sum_value <- sum(group$base_case)
    cat("Sum for", unique(group$model_param), "group:", sum_value, "\n")
    return(sum_value)
  })
}

# Check the distributions in the 'dist' dataframe
stopifnot(all(abs(check_distributions(dist) - 1) < 1e-8))

#######  3. SAMPLING UTILITIES, COSTS, MORT  ######

#' Title : Generate Samples from A Truncated Normal Distribution
#'
#' @param mean a numeric value representing the mean of the distribution
#' @param lower a numeric value representing the low 95% CI
#' @param upper a numeric value representing the high 95% CI
#' @param n a numeric value representing the number of samples to generate
#' @param upperbound a numeric value representing the upper bound used in truncation
#'
#' @return a numeric vector of samples generated from the truncated normal distribution
#' @export
generate_samples_trunc <- function(mean, lower, upper,upperbound, n) {
  
  # Input checks
  assert_that(is.numeric(lower), is.numeric(upper), is.numeric(mean),
              length(lower) == 1, length(upper) == 1, length(mean)==1)
  
  if (is.na(lower) | is.na(upper)|is.na(upperbound)|is.na(mean)) {
    stop("Lower and/or upper bound(s) and/or mean = NA. Cannot proceed.")
  }
  
  if (lower > upper) {
    stop("Lower bound is greater than upper bound.")
  }
  
  if (mean < lower || mean > upper) {
    stop("Inconsistent CI: mean is outside the interval")
  }
  
  # Calculate the standard deviation based on 95% CI
  sd <- (upper - lower) / (2 * 1.96)
  
  # Generate truncated normal samples with stated upperbound
  return(rtruncnorm(n, b = upperbound, mean = mean, sd = sd))
}

#' Title: Generate Gamma Samples Based on Mean and 95% CI (Quantile Matching)
#'
#' @param mean The mean of the distribution
#' @param lower The 2.5% quantile (lower bound of 95% CI)
#' @param upper The 97.5% quantile (upper bound of 95% CI)
#' @param n Number of samples to draw
#'
#' @return A numeric vector of samples from a fitted gamma distribution
#' @export
generate_gamma_samples <- function(mean, lower, upper, n) {
  # Input checks
  assert_that(is.numeric(lower), is.numeric(upper), is.numeric(mean),
              length(lower) == 1, length(upper) == 1, length(mean)==1)
  
  if (is.na(lower) | is.na(upper) | is.na(mean)) {
    stop("One or more of mean/lower/upper is NA.")
  }
  
  if (lower >= upper) {
    stop("Lower bound must be less than upper bound.")
  }
  
  if (mean <= 0) {
    stop("Mean must be positive for gamma distribution.")
  }
  if (mean < lower || mean > upper) {
    stop("Inconsistent CI: mean is outside the interval")
  }
  # Objective function: squared error of quantiles and mean
  obj_fun <- function(par) {
    shape <- par[1]
    scale <- par[2]
    
    diff1 <- qgamma(0.025, shape = shape, scale = scale) - lower
    diff2 <- qgamma(0.975, shape = shape, scale = scale) - upper
    diff3 <- (shape * scale) - mean 
    
    return(diff1^2 + diff2^2 + diff3^2)
  }
  
  # Initial estimates
  sd_guess <- (upper - lower) / (2 * 1.96)
  shape_init <- (mean / sd_guess)^2
  scale_init <- sd_guess^2 / mean
  
  opt <- optim(par = c(shape_init, scale_init),
               fn = obj_fun,
               method = "L-BFGS-B",
               lower = c(1e-5, 1e-5),
               control = list(maxit = 1000))
  
  if (opt$convergence != 0) {
    warning(sprintf("Gamma parameter optimization did not converge. Reverting to SD approximation. [mean=%.2f, lower=%.2f, upper=%.2f]", 
                    mean, lower, upper))
    shape <- shape_init
    scale <- scale_init
  } else {
    shape <- opt$par[1]
    scale <- opt$par[2]
  }
  
  # Draw gamma samples
  return(rgamma(n, shape = shape, scale = scale))
}

#' Title: Generate log-normal samples (based on 95% CIs)
#'  parameterized to hit 2.5%/97.5% exactly; implied mean may differ
#' @param lower 2.5% bound (>0) on original scale
#' @param upper 97.5% bound (>0) on original scale
#' @param n     number of samples
#' @return numeric vector of length n
generate_lognormal_samples <- function(lower, upper, n) {
  # Input validation
  assertthat::assert_that(is.numeric(lower), is.numeric(upper),length(lower) == 1, length(upper) == 1)
  if (isTRUE(anyNA(c(lower, upper)))) stop("NA values detected.")
  if (lower <= 0 || upper <= 0) stop("Lognormal requires positive values.")
  if (lower > upper) stop("Lower bound cannot be greater than the upper bound.")
  
  # Log-transform the CI bounds
  log_lower <- log(lower)
  log_upper <- log(upper)
  
  # Calculate meanlog and sdlog from the log-transformed CI
  # Assumes a 95% CI on the log scale
  z_alpha_2 <- qnorm(0.975)  # 1.96 for a 95% CI
  sdlog <- (log_upper - log_lower) / (2 * z_alpha_2)
  meanlog <- (log_lower + log_upper) / 2
  
  # Generate and return samples
  return(rlnorm(n, meanlog = meanlog, sdlog = sdlog))
}

# Generate utility, cost and RR for mortality samples for each mRS score
util_samples <- list()
cost_samples <- list()
mort_samples <- list()

for (i in 1:nrow(utility)) {
  # Utility samples
  mean_value <- utility$base_case[i]
  utility_lower <- utility$PSA_low[i]
  utility_upper <- utility$PSA_high[i]
  util_samples[[i]] <- generate_samples_trunc(mean = mean_value, lower = utility_lower, upper = utility_upper,
                                        upperbound=1, n = n.sample)
}

for (i in 1:nrow(cost)) {
  # Cost samples
  cost_mean <- cost$base_case[i]
  cost_lower <- cost$PSA_low[i]
  cost_upper <- cost$PSA_high[i]
  cost_samples[[i]] <- generate_gamma_samples(mean = cost_mean, lower = cost_lower, upper = cost_upper, n = n.sample)
}

for (i in 1:nrow(mort)) {
  # Cost samples
  mort_mean <- mort$base_case[i]
  mort_lower <- mort$PSA_low[i]
  mort_upper <- mort$PSA_high[i]
  mort_samples[[i]] <- generate_lognormal_samples(lower = mort_lower, upper = mort_upper, n = n.sample)
}

# Combine all utility samples into a data frame
util_sample_df <- data.frame(
  mrs = rep(utility$mrs, each = n.sample),
  utility_sample = unlist(util_samples),
  sample_id = rep(1:n.sample, times = nrow(utility))
)

# Manually add mRS = 6 utility samples (all zero)
mrs6_util_df <- data.frame(
  mrs = 6,
  utility_sample = rep(0, n.sample),
  sample_id = 1:n.sample
)

# Append to utility sample data
util_sample_df <- bind_rows(util_sample_df, mrs6_util_df)

# Combine all cost samples into a data frame
cost_sample_df <- data.frame(
  mrs = rep(cost$mrs, each = n.sample),
  cost_sample = unlist(cost_samples),
  sample_id = rep(1:n.sample, times = nrow(cost))
)

# Manually add mRS = 6 samples (all zero)
mrs6_cost_df <- data.frame(
  mrs = 6,
  cost_sample = rep(0, n.sample),
  sample_id = 1:n.sample
)

# Append to cost sample data
cost_sample_df <- bind_rows(cost_sample_df, mrs6_cost_df)

# clean
rm(mrs6_cost_df)
rm(mrs6_util_df)

# Combine all RR samples into a data frame
mort_sample_df <- data.frame(
  mrs = rep(mort$mrs, each = n.sample),
  mort_sample = unlist(mort_samples),
  sample_id = rep(1:n.sample, times = nrow(mort))
)

# add equivalent for mrs 0,1 and 6
mort_sample_df <- rbind(mort_sample_df, data.frame(mrs = 0, mort_sample = 1, sample_id = 1:n.sample))
mort_sample_df <- rbind(mort_sample_df, data.frame(mrs = 1, mort_sample = 1, sample_id = 1:n.sample))
mort_sample_df <- rbind(mort_sample_df, data.frame(mrs = 6, mort_sample = NA, sample_id = 1:n.sample))

# checks
assert_that(all(util_sample_df$utility_sample <= 1),
  msg = "Utilities must <= 1 (inclusive).") ## !!! note our input data do have negative utility values
assert_that(all(mort_sample_df$mort_sample > 0 | is.na(mort_sample_df$mort_sample)),
  msg = "Mortality samples must be > 0 or NA.")
assert_that(all(cost_sample_df$cost_sample >= 0),
  msg = "Cost samples must be non-negative.")

# combine mrs sampling data for utility, costs and mortality to read into model
mrs_samples <- mort_sample_df %>%
  left_join(cost_sample_df, by = c("sample_id","mrs")) %>%
  left_join(util_sample_df, by = c("sample_id","mrs"))

# create table of values for deterministic model 

## use this to check values make sense
mrs_samples_mean <- mrs_samples %>%
  group_by(mrs) %>%
  summarise(utility_sample = median(utility_sample, na.rm = TRUE),
            cost_sample = median(cost_sample, na.rm = TRUE),
            mort_sample = median(mort_sample, na.rm = TRUE)) %>%
  as.data.table()

## replace for deterministic values for use in deterministic model
mrs_samples_mean[mrs_samples_mean$mrs == 0, "cost_sample"] <- parameters[model_param=="cost.mrs"& mrs==0,base_case]
mrs_samples_mean[mrs_samples_mean$mrs == 1, "cost_sample"] <- parameters[model_param=="cost.mrs"& mrs==1,base_case]
mrs_samples_mean[mrs_samples_mean$mrs == 2, "cost_sample"] <- parameters[model_param=="cost.mrs"& mrs==2,base_case]
mrs_samples_mean[mrs_samples_mean$mrs == 3, "cost_sample"] <- parameters[model_param=="cost.mrs"& mrs==3,base_case]
mrs_samples_mean[mrs_samples_mean$mrs == 4, "cost_sample"] <- parameters[model_param=="cost.mrs"& mrs==4,base_case]
mrs_samples_mean[mrs_samples_mean$mrs == 5, "cost_sample"] <- parameters[model_param=="cost.mrs"& mrs==5,base_case]
mrs_samples_mean[mrs_samples_mean$mrs == 6, "cost_sample"] <- parameters[model_param=="cost.mrs"& mrs==6,base_case]

# replace with deterministic values for mortality
mrs_samples_mean[mrs_samples_mean$mrs == 0, "mort_sample"] <- parameters[model_param=="rr.mort"& mrs==0,base_case]
mrs_samples_mean[mrs_samples_mean$mrs == 1, "mort_sample"] <- parameters[model_param=="rr.mort"& mrs==1,base_case]
mrs_samples_mean[mrs_samples_mean$mrs == 2, "mort_sample"] <- parameters[model_param=="rr.mort"& mrs==2,base_case]
mrs_samples_mean[mrs_samples_mean$mrs == 3, "mort_sample"] <- parameters[model_param=="rr.mort"& mrs==3,base_case]
mrs_samples_mean[mrs_samples_mean$mrs == 4, "mort_sample"] <- parameters[model_param=="rr.mort"& mrs==4,base_case]
mrs_samples_mean[mrs_samples_mean$mrs == 5, "mort_sample"] <- parameters[model_param=="rr.mort"& mrs==5,base_case]

# replace with deterministic values for utility
mrs_samples_mean[mrs_samples_mean$mrs == 0, "utility_sample"] <- parameters[model_param=="utility.mrs"& mrs==0,base_case]
mrs_samples_mean[mrs_samples_mean$mrs == 1, "utility_sample"] <- parameters[model_param=="utility.mrs"& mrs==1,base_case]
mrs_samples_mean[mrs_samples_mean$mrs == 2, "utility_sample"] <- parameters[model_param=="utility.mrs"& mrs==2,base_case]
mrs_samples_mean[mrs_samples_mean$mrs == 3, "utility_sample"] <- parameters[model_param=="utility.mrs"& mrs==3,base_case]
mrs_samples_mean[mrs_samples_mean$mrs == 4, "utility_sample"] <- parameters[model_param=="utility.mrs"& mrs==4,base_case]
mrs_samples_mean[mrs_samples_mean$mrs == 5, "utility_sample"] <- parameters[model_param=="utility.mrs"& mrs==5,base_case]
mrs_samples_mean[mrs_samples_mean$mrs == 6, "utility_sample"] <- parameters[model_param=="utility.mrs"& mrs==6,base_case]

# save files
mrs_samples_mean <- as.data.table(mrs_samples_mean)
save(mrs_samples_mean, file="inputs/created_inputs/mrs_samples_mean.RData")
save(mrs_samples, file="inputs/created_inputs/mrs_samples.RData")


######   4. SAMPLE OTHER PSA PARAMETERS   ####

## only sample parameters with distributions
params.psa <- parameters[!is.na(PSA_distribution) & PSA_distribution != "" & PSA_distribution != "NA"]

## remove ones already have samples for/will be sampling separately
already.sampled <- c("q.ivt","q.mt","c.ivt.lt",
                     "c.mt.lt","utility.mrs","cost.mrs","or.mt","or.ivt",
                     "rr.mort")

params.psa <- params.psa[!(model_param %in% already.sampled)]

### since same parameter for some ASC and CSC need to rename
params.psa[, param_id := paste(model_param, Presentation.Setting, Intervention, sep = ",")]

params.psa.sample <- params.psa[rep(params.psa[, .I], n.sample)] ## create n.samples of each row
params.psa.sample$sample_id <- rep(1:n.sample, each=nrow(params.psa))

# Split parameters by distribution
unique(params.psa$PSA_distribution)
params_gamma   <- params.psa[PSA_distribution == "gamma"]
params_beta <- params.psa[PSA_distribution == "beta"]
params_uniform <- params.psa[PSA_distribution == "uniform"]
### !!! if add other distributions will have to add here and below

# Initialize an empty list to hold samples
beta_samples_list <- list()

## Checks
### !!! note for the current data set beta sampled values use k (in PSA_low)
### if you have data with known alpha and beta values this code will need adapting
# 1) All Beta rows have valid k
stopifnot(params_beta[, all(is.finite(as.numeric(PSA_low)) & as.numeric(PSA_low) > 0)])

# 2) All base_case for Beta rows are probabilities
stopifnot(params_beta[, all(is.finite(as.numeric(base_case)) &
                              as.numeric(base_case) >= 0 & as.numeric(base_case) <= 1)])

beta_samples <- params_beta[, {
  k <- PSA_low
  if (any(!is.finite(k) | k <= 0)) {
    stop("For Beta PSA rows, PSA_low must be finite and > 0 (k).")
  }
  alpha <- base_case * k
  beta_val <- (1 - base_case) * k
  .(sample_id = 1:n.sample,
    value = rbeta(n.sample, shape1 = alpha, shape2 = beta_val))
}, by = .(model_param, Presentation.Setting, Intervention)]

# Sample from gamma distribution 
# Prepare empty list to store results
gamma_samples_list <- list()

# Loop through each parameter to apply the generate_gamma_samples() function
for (i in seq_len(nrow(params_gamma))) {
  row <- params_gamma[i]
  draws <- generate_gamma_samples(mean = row$base_case,
                                  lower = row$PSA_low,
                                  upper = row$PSA_high,
                                  n = n.sample)
  
  gamma_samples_list[[i]] <- data.table(
    model_param = row$model_param,
    Presentation.Setting = row$Presentation.Setting,
    Intervention = row$Intervention,
    sample_id = 1:n.sample,
    value = draws
  )
}

# Combine the list into one data.table
gamma_samples <- rbindlist(gamma_samples_list)

# Sample from uniform distribution
uniform_samples <- params_uniform[, .(sample_id = 1:n.sample,
                                      value = runif(n.sample, PSA_low, PSA_high)),
                                  by = .(model_param, Presentation.Setting, Intervention)]

# Combine sampled values
sampled_values <- rbind(beta_samples, gamma_samples)
sampled_values <- rbind(sampled_values, uniform_samples)

# Add composite key for merge
sampled_values[, param_id := paste(model_param, Presentation.Setting, Intervention, sep = ",")]

# Merge sampled values back in
params.psa.sample <- merge(params.psa.sample, sampled_values[, .(param_id, sample_id, value)],
                           by = c("param_id", "sample_id"), all.x = TRUE)

### reshape to match other samples 
params.psa.sample <- params.psa.sample %>%
  select(sample_id, param_id, value) %>%
  pivot_wider(names_from = param_id, values_from = value)

## get IQR for sample for DSA for c.lvo and c.mt
parameters[model_param=="c.lvo", DSA_low := quantile(unlist(params.psa.sample$`c.lvo,all,all`), 0.25)]
parameters[model_param=="c.lvo", DSA_high := quantile(unlist(params.psa.sample$`c.lvo,all,all`), 0.75)]
parameters[model_param=="c.mt", DSA_low := quantile(unlist(params.psa.sample$`c.mt,all,all`), 0.25)]
parameters[model_param=="c.mt", DSA_high := quantile(unlist(params.psa.sample$`c.mt,all,all`), 0.75)]

########   5. TREATMENT PROBABILITIES SAMPLING ####

# Extract mean odds ratio (OR) and 95% confidence interval (CI) & convert to log
log_mean_or_mt_asc <- log(parameters[model_param=="or.mt" & 
                                       Presentation.Setting=="ASC", base_case])

log_low_or_mt_asc <- log(parameters[model_param=="or.mt" & 
                                      Presentation.Setting=="ASC", PSA_low])

log_hi_or_mt_asc <- log(parameters[model_param=="or.mt" & 
                                     Presentation.Setting=="ASC", PSA_high])

log_mean_or_mt_csc <- log(parameters[model_param=="or.mt" & 
                                       Presentation.Setting=="CSC", base_case])

log_low_or_mt_csc <- log(parameters[model_param=="or.mt" & 
                                      Presentation.Setting=="CSC", PSA_low])

log_hi_or_mt_csc <- log(parameters[model_param=="or.mt" & 
                                     Presentation.Setting=="CSC", PSA_high])

log_mean_or_ivt <- log(parameters[model_param=="or.ivt", base_case])

log_low_or_ivt <- log(parameters[model_param=="or.ivt", PSA_low])

log_hi_or_ivt <- log(parameters[model_param=="or.ivt", PSA_high])

# Calculate standard error from the 95% CI
# Preferred SE calculation: full CI width / (2 * 1.96)
se_log_or_mt_asc <- (log_hi_or_mt_asc - log_low_or_mt_asc) / (2 * 1.96)
se_log_or_mt_csc <- (log_hi_or_mt_csc - log_low_or_mt_csc) / (2 * 1.96)
se_log_or_ivt    <- (log_hi_or_ivt    - log_low_or_ivt)    / (2 * 1.96)

# Sample log-OR values from a normal distribution
log_smpl_mt_asc <- rnorm(n.sample, mean = log_mean_or_mt_asc, sd = se_log_or_mt_asc)
log_smpl_mt_csc <- rnorm(n.sample, mean = log_mean_or_mt_csc, sd = se_log_or_mt_csc)
log_smpl_ivt <- rnorm(n.sample, mean = log_mean_or_ivt, sd = se_log_or_ivt)

# Convert the log-OR samples back to the OR scale
smpl_mt_asc <- exp(log_smpl_mt_asc)
smpl_mt_csc <- exp(log_smpl_mt_csc)
smpl_ivt <- exp(log_smpl_ivt)

# tidy environment
rm(log_mean_or_mt_asc, log_low_or_mt_asc, log_hi_or_mt_asc, 
   log_mean_or_mt_csc, log_low_or_mt_csc, log_hi_or_mt_csc, log_mean_or_ivt, 
   log_low_or_ivt, log_hi_or_ivt, se_log_or_mt_asc, se_log_or_mt_csc, se_log_or_ivt, 
   log_smpl_mt_asc, log_smpl_mt_csc, log_smpl_ivt)

# Calculate treatment probability from control probability and OR
calculate_ptreatment <- function(p_control, OR) {
  ### standard formulae + clamp to avoid odds explosions if a control prob is 0 or 1
  p_control <- pmin(pmax(p_control, 1e-9), 1 - 1e-9)
  odds_control <- p_control / (1 - p_control)
  odds_treatment <- OR * odds_control
  p <- odds_treatment / (1 + odds_treatment)
  return(pmin(pmax(p, 1e-12), 1 - 1e-12)) ## guards from probabilities being <0, >1 and exactly=0 or 1
}

# Construct p_treatment table using sampled control probabilities
p_treatment <- data.table(sample_id = 1:n.sample)

### check have the columns necessary
ctrl_cols <- c(
  "p.noivt.emt2mt,early;ASC,no intervention",
  "p.ivt.emt2mt,early;ASC,no intervention",
  "p.emt2mt,late;ASC,no intervention",
  "p.noivt.emt2mt,early;CSC,no intervention",
  "p.ivt.emt2mt,early;CSC,no intervention",
  "p.emt2mt,late;CSC,no intervention",
  "p.eivt2ivt,early;ASC,no intervention",
  "p.eivt2ivt,early;CSC,no intervention"
)
stopifnot(all(ctrl_cols %in% names(params.psa.sample)))


# ASC
p_treatment[, `p.noivt.emt2mt,early;ASC,intervention` := calculate_ptreatment(
  params.psa.sample[["p.noivt.emt2mt,early;ASC,no intervention"]], smpl_mt_asc)]
p_treatment[, `p.ivt.emt2mt,early;ASC,intervention` := calculate_ptreatment(
  params.psa.sample[["p.ivt.emt2mt,early;ASC,no intervention"]], smpl_mt_asc)]
p_treatment[, `p.emt2mt,late;ASC,intervention` := calculate_ptreatment(
  params.psa.sample[["p.emt2mt,late;ASC,no intervention"]], smpl_mt_asc)]

# CSC
p_treatment[, `p.noivt.emt2mt,early;CSC,intervention` := calculate_ptreatment(
  params.psa.sample[["p.noivt.emt2mt,early;CSC,no intervention"]], smpl_mt_csc)]
p_treatment[, `p.ivt.emt2mt,early;CSC,intervention` := calculate_ptreatment(
  params.psa.sample[["p.ivt.emt2mt,early;CSC,no intervention"]], smpl_mt_csc)]
p_treatment[, `p.emt2mt,late;CSC,intervention` := calculate_ptreatment(
  params.psa.sample[["p.emt2mt,late;CSC,no intervention"]], smpl_mt_csc)]

# ivt
p_treatment[, `p.eivt2ivt,early;ASC,intervention` := calculate_ptreatment(
  params.psa.sample[["p.eivt2ivt,early;ASC,no intervention"]], smpl_ivt)]
p_treatment[, `p.eivt2ivt,early;CSC,intervention` := calculate_ptreatment(
  params.psa.sample[["p.eivt2ivt,early;CSC,no intervention"]], smpl_ivt)]

# Merge treatment probabilities into PSA sample
params.psa.sample <- params.psa.sample[,which(unlist(lapply(params.psa.sample, function(x)!all(is.na(x))))),with=F]

#### check no duplicates
stopifnot(!anyDuplicated(names(params.psa.sample)))

# Add treatment probabilities to the params.psa.sample
params.psa.sample <- params.psa.sample %>%
  left_join(p_treatment, by = "sample_id") %>%
  as.data.table()

# View the summary statistics for each group as a check
sum_stats_t <- params.psa.sample %>%
  dplyr::select(-sample_id) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), list(
    mean     = ~mean(., na.rm = TRUE),
    lower_CI = ~quantile(., 0.025, na.rm = TRUE),
    upper_CI = ~quantile(., 0.975, na.rm = TRUE)
  )))
sum_stats_t

#### update parameters dt to include the descriptive stats
# Function to update parameter values in the data table

#' Update PSA Parameter Values from Summary Statistics
#'
#' Updates the `parameters_dt` data.table with new `base_case`, `PSA_low`, and `PSA_high`
#' values for a list of model parameters, using summary statistics from a PSA sampling result.
#'
#' @param params_vector A character vector of parameter identifiers, each in the format 
#'   `"model_param,Presentation.Setting, Intervention"` (e.g., `"p.noivt.emt2mt,early;ASC,intervention"`).  
#' @param sum_stats_t A data frame or tibble containing PSA summary statistics.
#'   Column names must follow the pattern `"{param}_mean"`, `"{param}_lower_CI"`, `"{param}_upper_CI"`.
#' @param parameters_dt A `data.table` of model parameters to be updated (modified in-place).
#'
#' @return Invisibly returns the updated `parameters_dt`. The object is modified by reference.
#' @export
update_parameters <- function(params_vector, sum_stats_t, parameters_dt) {
  
  # Input checks
  stopifnot(is.character(params_vector))
  stopifnot(is.data.frame(sum_stats_t))
  stopifnot("data.table" %in% class(parameters_dt))
  
  # Check all required columns exist in parameters_dt
  required_cols <- c("model_param", "Presentation.Setting", "Intervention", "base_case", "PSA_low", "PSA_high")
  missing_cols <- setdiff(required_cols, names(parameters_dt))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns in parameters_dt:", paste(missing_cols, collapse = ", ")))
  }
  
  # Iterate over each parameter string in the vector
  for (param_str in params_vector) {
    
    # Split the string into model_param and Presentation.Setting
    split_param <- strsplit(param_str, ",")[[1]]
    
    model_param_name       <- split_param[1]            # e.g. "p.noivt.emt2mt"
    presentation_setting   <- split_param[2]            # e.g. "early;ASC"
    intervention_arm     <- split_param[3]              # e.g. "no intervention"
    
    # Build column names to pull from sum_stats_t
    mean_col   <- paste0(param_str, "_mean")
    lower_col  <- paste0(param_str, "_lower_CI")
    upper_col  <- paste0(param_str, "_upper_CI")
    
    if (!(mean_col %in% names(sum_stats_t))) {
      stop(paste("Missing summary statistic:", mean_col))
    }
    
    # Do the update in the parameters data.table
    parameters_dt[model_param == model_param_name &
                    Presentation.Setting == presentation_setting &
                    Intervention == intervention_arm, 
                  
                  `:=`(
                    base_case = sum_stats_t[[mean_col]],
                    PSA_low   = sum_stats_t[[lower_col]],
                    PSA_high  = sum_stats_t[[upper_col]]
                  )
    ]
  }
  return(parameters_dt)  # Return the modified data.table 
}



# Build list of all parameter keys from sum_stats_t column names
all_params <- names(sum_stats_t) %>%
  str_remove("_(mean|lower_CI|upper_CI)") %>%
  unique()

# Run the update in one go
update_parameters(all_params, sum_stats_t, parameters)

###### checks 
# Utilities can be negative but must be <= 1
stopifnot(all(util_sample_df$utility_sample <= 1, na.rm = TRUE))

# Probabilities 
# Identify probability columns
prob_cols <- names(params.psa.sample)[grepl("^p\\.", names(params.psa.sample))]
# Flatten to a single numeric vector
prob_values <- unlist(params.psa.sample[, ..prob_cols])
# Assert all finite and within [0,1]
assert_that(
  all(is.finite(prob_values)),
  all(!is.na(prob_values)),
  all(prob_values >= 0 & prob_values <= 1),
  msg = sprintf(
    "Probability values must be finite and within [0,1]. Out-of-range examples: %s",
    paste(head(prob_values[prob_values < 0 | prob_values > 1]), collapse = ", ")
  )
)

# Mortality RR positive (NA allowed for mRS=6)
stopifnot(all(mort_sample_df$mort_sample > 0 | is.na(mort_sample_df$mort_sample)))

######   6. FINAL SAVE   #####

save(parameters, file="inputs/created_inputs/parameters_edited.RData")
save(params.psa.sample, file="inputs/created_inputs/params.psa.sample.RData")

write.csv(parameters, file="inputs/created_inputs/parameters_post_psa.csv")

####### 7. QUICK VISUAL QC PLOTS ( no export) #######

## 7A. params.sample
#### checking plots of distribution of cost for mRS
cols <- setdiff(names(params.psa.sample), "sample_id")
num_cols <- cols[sapply(cols, function(nm) is.numeric(params.psa.sample[[nm]]))]

# adjust rows/cols if you want a different grid
par(mfrow = c(3, 3), mar = c(3, 3, 2, 1))

for (nm in num_cols) {
  x <- params.psa.sample[[nm]]
  x <- x[is.finite(x)]
  if (!length(x)) { plot.new(); title(main = nm); next }
  hist(x, main = gsub(",", " | ", nm), xlab = "")
}

## 7B. by mrs
####### 7. SIMPLE HISTOGRAMS FOR mrs_samples BY mRS #######

metrics <- intersect(c("utility_sample", "cost_sample", "mort_sample"),
                     names(mrs_samples))
mrs_vals <- sort(unique(mrs_samples$mrs))

for (metric in metrics) {
  # grid; tweak if you want a different layout
  par(mfrow = c(3, 3), mar = c(3, 3, 2, 1))
  for (mm in mrs_vals) {
    xx <- mrs_samples[mrs_samples$mrs == mm, metric, drop = TRUE]
    xx <- xx[is.finite(xx)]
    if (!length(xx)) { plot.new(); title(main = paste(metric, "~ mRS", mm)); next }
    hist(xx, main = paste(metric, "~ mRS", mm), xlab = "")
  }
}

## testing
temp_test0 <- subset(mrs_samples,mrs==0)
unique(temp_test0$mort_sample)
temp_test6 <- subset(mrs_samples,mrs==6)
unique(temp_test6$utility_sample)
unique(temp_test6$cost_sample)
unique(temp_test6$mort_sample)