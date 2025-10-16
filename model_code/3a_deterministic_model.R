###############################################
# TITLE: Running the Deterministic Model Runs for B360S Model
# AUTHOR: Nichola Naylor (OI Pharma Partners Ltd), aided by GPT-4o,GTP-5 & Github co-pilot
# DATE: September 2025
#
# DESCRIPTION:
# This script executes the deterministic base case and deterministic sensitivity analysis (DSA)
# for the B360S stroke care cost-effectiveness model. It evaluates the Net Monetary Benefit (NMB)
# under base assumptions and under parameter variations (low/high bounds), and visualizes results
# in a tornado plot.
#
# INPUTS:
# - parameters_edited.RData: a .RData file containing model parameters.
# - mrs_samples_mean.RData: a .RData file containing mean modified Rankin Scale mortality samples.
# - model_functions.R: source file containing the core health economic model (`run_model()`).
#
# OUTPUTS:
# - outputs/base_case_results.RData: results of the base case model run.
# - outputs/population_results.csv: total population results for base case.
# - outputs/NMB.csv: incremental Net Monetary Benefit (NMB) for base case.
# - outputs/proc_results.csv: process results for base case.
# - outputs/mrstrace_results.csv: modified Rankin Scale trace results for base case.
# - outputs/DSA_results.RData: results of the deterministic sensitivity analysis.


#### ======================================= ####
####       1. INITIALISE & LOAD LIBRARIES      ####
#### ======================================= ####

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

#### ======================================= ####
####       2. LOAD DATA & R SCRIPTS         ####
#### ======================================= ####

load("inputs/created_inputs/parameters_edited.RData")
data_main <- parameters
load("inputs/created_inputs/mrs_samples_mean.RData")
source("model_code/model_functions.R")


#### ======================================= ####
####       3. RUN BASE CASE                 ####
#### ======================================= ####

base_case <- run_model(data_main,cycles=10,mrs_samples_mean)


save(base_case, file = "outputs/base_case_results.RData")

write.csv(base_case$summary_data_all, file="outputs/population_results.csv")
write.csv(base_case$incremental_results, file="outputs/NMB.csv")
write.csv(base_case$process_results, file="outputs/proc_results.csv")
write.csv(base_case$mrs.trace, file="outputs/mrstrace_results.csv")

#### =================================================== ####
####       4. RUN DETERMINISTIC SENSITIVTY ANALYSIS     ####
#### =================================================== ####

qaly.threshold <- data_main[model_param=="wtp",base_case]
stopifnot(length(qaly.threshold)==1, is.finite(qaly.threshold), qaly.threshold >= 0)

params <- data_main[`DSA_flag`=="y"]

# Initialize a list to store results
low.results <- list()
high.results <- list()

for (i in 1:nrow(params)) {
  
  variable_name <- params$model_param[i]
  presentation <- params$Presentation.Setting[i]
  intervention <- params$Intervention[i]
  
  
  ### run low
  temp <- copy(data_main)
  
  # Update the parameter value for sensitivity analysis
  temp[model_param==variable_name &
         Presentation.Setting==presentation &
         Intervention==intervention, base_case := DSA_low] 
  
  # Run the health economic model with the updated parameter
  result <- run_model(temp, cycles=10, mrs_samples_mean)
  result <- (qaly.threshold * result$incremental_results$inc.qol) - result$incremental_results$inc.cost
  
  rm(temp)
  
  # Store the result
  low.results[[variable_name]] <- result
  
  ### run high
  temp <- copy(data_main)
  
  # Update the parameter value for sensitivity analysis
  temp[model_param==variable_name&
         Presentation.Setting==presentation&
         Intervention==intervention, base_case := DSA_high]
  
  # Run the health economic model with the updated parameter
  result <- run_model(temp, cycles=10, mrs_samples_mean)
  result <- (qaly.threshold * result$incremental_results$inc.qol) - result$incremental_results$inc.cost
  
  rm(temp)
  
  # Store the result
  high.results[[variable_name]] <- result
  
}

####### tornado plot 

low.results <- data.frame(low.results)

low.results <- low.results %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",      
    values_to = "low.value"         #
  )


high.results <- data.frame(high.results)

high.results <- high.results %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",      
    values_to = "high.value"         
  )

DSA <- merge(low.results, high.results,by="variable",all=TRUE)

### get a column can use in tornado plot 
params[ , low.input := round(DSA_low,2)]
params[ , high.input := round(DSA_high,2)]

params[ , Parameter := paste0(Description," [",low.input," , ", high.input, "]")]
params <- params[,c("Parameter","model_param")]

DSA <- merge(DSA,params,by.x="variable",by.y="model_param")

# Example data frame

# Calculate midpoints for plotting
DSA <- DSA %>%
  mutate(midpoint = (low.value + high.value) / 2,
         ymin = pmin(low.value, high.value),
         ymax = pmax(low.value, high.value))

save(DSA, file = "outputs/DSA_results.RData")
