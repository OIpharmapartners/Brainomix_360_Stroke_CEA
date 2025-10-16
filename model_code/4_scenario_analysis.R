###############################################
# TITLE: Scenario Analyses for B360S Model
# AUTHOR: Nichola Naylor (OI Pharma Partners Ltd), aided by GPT-4o,GTP-5 & Github co-pilot
# DATE: September 2025
#
# DESCRIPTION:
# Runs predefined structural and policy scenarios to test alternative assumptions
# The scenarios are:
# 1. Use 66 as start age
# 2. Remove long term cost savings 
# 3. scenario 2. + Mortality impacts only occur in the first year
# 4. Long term QALY and cost impacts are additive for IVT and MT
# 5. A different IVT MRS distribution
#
# WHAT THIS SCRIPT DOES:
# 1) Defines scenario sets (named lists) with altered parameters/flags.
# 2) For each scenario, updates `data_main` and executes `run_model()`.
# 3) Collates results into a tidy scenario summary and exports figures:
#    - Scenario league table
#    - DSA tornado (if ranges provided)
#    - Optional threshold plots (NMB vs key parameter)
#
# KEY INPUTS (read-only):
# - inputs/created_inputs/parameters_edited.RData   # baseline data_main
# - inputs/created_inputs/mrs_samples_mean.RData
# - config/scenarios.yaml or .csv (optional)        # scenario definitions
#
# KEY OUTPUTS (write):
# - outputs/scenario_table.csv: inc.cost, inc.QALY, ICER, NMB by scenario
# - outputs/tornado.csv / tornado.png               # if DSA specified
# - outputs/threshold.csv / threshold.png           # if thresholding used
#
# DEPENDENCIES:
# - core model functions: mrs_markov(), run_model()
# - PSA outputs optional for overlay (e.g., CEAC under scenario)
# - Packages: data.table, dplyr, ggplot2 (if plotting), assertthat, conflicted
#
#############################################

#### ======================================= ####
####       INITIALISE & LOAD LIBRARIES      ####
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
library(scales)

# Resolve potential function conflicts (tidyverse vs data.table vs base)
conflict_prefer("merge", "data.table")
conflict_prefer("filter", "dplyr")

load("inputs/created_inputs/parameters_edited.RData")
data_main <- parameters
load("inputs/created_inputs/mrs_samples_mean.RData")
source("model_code/model_functions.R")

load("outputs/base_case_results.RData")

base_case$incremental_results$scenario <- "base_case"
base_case$process_results$scenario <- "base_case"

#### ======================================= ####
####       1.CHANGE START AGE                 ####
#### ======================================= ####

sc1_data_main <- copy(data_main)
sc1_data_main[model_param=="age",base_case := 66] 

# Check that mortality tables include ages ≥ start age (66)
max_age_f <- max(sc1_data_main[model_param == "fqx", age], na.rm = TRUE)
max_age_m <- max(sc1_data_main[model_param == "mqx", age], na.rm = TRUE)

assert_that(
  max_age_f >= 66 && max_age_m >= 66,
  msg = sprintf(
    "Mortality tables do not cover age 66. Max female age = %s, max male age = %s.",
    max_age_f, max_age_m
  )
)

sc1 <- run_model(sc1_data_main,cycles=10,mrs_samples_mean)

sc1$incremental_results$scenario <- "sc1"
sc1$process_results$scenario <- "sc1"

#### ======================================= ####
####       2.REMOVE LT COST SAVINGS               ####
#### ======================================= ####
sc2_data_main <- copy(data_main)
sc2_mrs_samples_mean <- copy(mrs_samples_mean)
sc2_mrs_samples_mean$cost_sample <- 0

sc2 <- run_model(sc2_data_main,cycles=10,sc2_mrs_samples_mean)

sc2$incremental_results$scenario <- "sc2"
sc2$process_results$scenario <- "sc2"

#### ======================================= ####
#### 3.REMOVE LT COST SAVINGS + MORTALITY IMPACTS IN FIRST YEAR ####
#### ======================================= ####

source("model_code/model_functions_mortality_sc.R")
sc3_data_main <- copy(data_main)
sc3 <- run_model_Msc(sc3_data_main,cycles=10,sc2_mrs_samples_mean)
sc3$incremental_results$scenario <- "sc3"
sc3$process_results$scenario <- "sc3"

#### ======================================= ####
#### 4. LT costs and QALYs are additive across MT and IVT ####
#### ======================================= ####

source("model_code/model_functions_MTIVT_sc.R")
sc4_data_main <- copy(data_main)
sc4 <- run_model_mtivtsc(sc4_data_main,cycles=10,mrs_samples_mean)
sc4$incremental_results$scenario <- "sc4"
sc4$process_results$scenario <- "sc4"

#### ======================================= ####
#### 5. Different mRS distribution post-IVT ####
#### ======================================= ####

sc5_data_main <- copy(data_main)

sc5_data_main[model_param=="dist.ivt" & mrs=="0",base_case := 0.154]
sc5_data_main[model_param=="dist.ivt" & mrs=="1",base_case := 0.173]
sc5_data_main[model_param=="dist.ivt" & mrs=="2",base_case := 0.061]
sc5_data_main[model_param=="dist.ivt" & mrs=="3",base_case := 0.228]
sc5_data_main[model_param=="dist.ivt" & mrs=="4",base_case := 0.099]
sc5_data_main[model_param=="dist.ivt" & mrs=="5",base_case := 0.067]
sc5_data_main[model_param=="dist.ivt" & mrs=="6",base_case := 0.218]

sc5 <- run_model(sc5_data_main,cycles=10,mrs_samples_mean)
sc5$incremental_results$scenario <- "sc5"
sc5$process_results$scenario <- "sc5"


#### ======================================= ####
#### COMBINE RESULTS                  ####
#### ======================================= ####

incremental.results <- rbind(base_case$incremental_results,sc1$incremental_results,sc2$incremental_results,
                             sc3$incremental_results,sc4$incremental_results,sc5$incremental_results)
process.results <- rbind(base_case$process_results,sc1$process_results,sc2$process_results,
                         sc3$process_results, sc4$process_results,sc5$process_results)
### edit dps
incremental.results <- as.data.table(incremental.results)
incremental.results[, inc.cost := comma(round(inc.cost, 0))]
incremental.results[, inc.qol := comma(round(inc.qol, 2))]
incremental.results[, NMB := comma(round(NMB, 0))]

process.results <- as.data.table(process.results)
process.results[, intervention := comma(round(intervention, 0))]
process.results[, standard := comma(round(standard, 0))]
process.results[, intervention_costs := comma(round(intervention_costs, 0))]
process.results[, standard_costs := comma(round(standard_costs, 0))]
process.results[, intervention_qol := comma(round(intervention_qol, 2))]
process.results[, standard_qol := comma(round(standard_qol, 2))]
process.results[, LT_cost := comma(round(LT_cost, 0))]
process.results[, LT_qol := comma(round(LT_qol, 2))]
process.results[, unit.cost := comma(round(unit.cost, 2))]

### save results
write.csv(incremental.results,"outputs/scenario_incremental_results.csv",row.names=FALSE)
write.csv(process.results,"outputs/scenario_process_results.csv",row.names=FALSE)
