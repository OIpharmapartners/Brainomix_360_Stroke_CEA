###############################################
# TITLE: Scenario Analyses for B360S Model
# AUTHOR: Nichola Naylor (OI Pharma Partners Ltd), aided by GPT-4o,GTP-5 & Github co-pilot
#
# DESCRIPTION:
# Runs predefined structural and policy scenarios to test alternative assumptions
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
# - inputs/created_inputs/parameters_edited.RData    # baseline data_main
# - inputs/created_inputs/mrs_samples_mean.RData
# - outputs/base_case_results.RData                  # base case for comparison
# - outputs/psa_outputs.csv                          # PSA draws (from 3b)
#
# KEY OUTPUTS (write):
# - outputs/scenario_incremental_results.csv         # inc.cost, inc.qol, NMB by scenario
# - outputs/scenario_process_results.csv             # per-procedure results by scenario
# - outputs/psa_summary.csv

#### !!! note scenarios are defined here and need recoding if changed

#### ======================================= ####
####       INITIALISE & LOAD LIBRARIES      ####
#### ======================================= ####
# clear environment
rm(list=ls())

# Load necessary packages
library(conflicted)
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

sc5_data_main[model_param=="dist.ivt" & mrs==0,base_case := 0.14157]
sc5_data_main[model_param=="dist.ivt" & mrs==1,base_case := 0.16769]
sc5_data_main[model_param=="dist.ivt" & mrs==2,base_case := 0.08384]
sc5_data_main[model_param=="dist.ivt" & mrs==3,base_case := 0.19125]
sc5_data_main[model_param=="dist.ivt" & mrs==4,base_case := 0.06035]
sc5_data_main[model_param=="dist.ivt" & mrs==5,base_case := 0.06035]
sc5_data_main[model_param=="dist.ivt" & mrs==6,base_case := 0.29495]


sc5 <- run_model(sc5_data_main,cycles=10,mrs_samples_mean)
sc5$incremental_results$scenario <- "sc5"
sc5$process_results$scenario <- "sc5"

#### ======================================= ####
#### 6. Gradual mortality attenuation (years 1-4, background mortality from year 5) ####
#### ======================================= ####

source("model_code/model_functions_mortality_gradual_sc.R")
sc6_data_main <- copy(data_main)
sc6 <- run_model_Msc_gradual(sc6_data_main, cycles = 10, mrs_samples_mean)
sc6$incremental_results$scenario <- "sc6"
sc6$process_results$scenario <- "sc6"


#### ======================================= ####
#### 7. Half-cycle correction applied        ####
#### ======================================= ####

source("model_code/model_functions_hcc.R")
sc7_data_main <- copy(data_main)
sc7 <- run_model_hcc(sc7_data_main, cycles = 10, mrs_samples_mean)
sc7$incremental_results$scenario <- "sc7"
sc7$process_results$scenario <- "sc7"

#### ======================================= ####
#### 8. Dampened intervention effect (10%)  ####
#### ======================================= ####
# The intervention's effect enters the model through pathway probabilities
# (e.g. p.eivt2ivt, p.ivt.emt2mt), we shrink the incremental effect
# of the intervention towards standard care by 10%, i.e.:
#   dampened = standard + 0.9 * (intervention - standard)

sc8_data_main <- copy(data_main)
damp_factor <- 0.9   # retain 90% of the intervention effect

# paired intervention / standard-care values for probability parameters
int_rows <- sc8_data_main[Intervention == "intervention" & grepl("^p\\.", model_param),
                          .(model_param, Presentation.Setting, mrs, int_val = base_case)]
std_rows <- sc8_data_main[Intervention == "no intervention" & grepl("^p\\.", model_param),
                          .(model_param, Presentation.Setting, mrs, std_val = base_case)]

effect_map <- merge(int_rows, std_rows,
                    by = c("model_param", "Presentation.Setting", "mrs"))

assert_that(nrow(effect_map) > 0,
            msg = "Scenario 8: no paired intervention/standard parameters found")

effect_map[, damp_val := std_val + damp_factor * (int_val - std_val)]

# sanity check: dampened values remain valid probabilities
assert_that(all(effect_map$damp_val >= 0 & effect_map$damp_val <= 1),
            msg = "Scenario 8: dampened probabilities outside [0,1]")

# overwrite the intervention-arm rows only
sc8_data_main[effect_map,
              on = c("model_param", "Presentation.Setting", "mrs"),
              base_case := fifelse(Intervention == "intervention", i.damp_val, base_case)]

print(effect_map[, .(model_param, Presentation.Setting, std_val, int_val, damp_val)])

sc8 <- run_model(sc8_data_main, cycles = 10, mrs_samples_mean)
sc8$incremental_results$scenario <- "sc8"
sc8$process_results$scenario <- "sc8"

#### ======================================= ####
#### 9. Per-AIS patient results               ####
#### ======================================= ####
# Number of AIS patients from base case imaging counts (should equal ~81,565)
n.ais <- sum(base_case$process_results$intervention[
  base_case$process_results$procedure %in% c("NCCT + CTA", "NCCT + CTA + CTP", "NCCT + CTA + CTP + MRI")
])

assert_that(isTRUE(all.equal(n.ais,
                             data_main[model_param=="n.stroke",base_case]*
                               data_main[model_param=="p.ais",base_case])),
            msg = "issue with AIS numbers in scenario 9")

sc9_inc <- data.frame(
  inc.cost = base_case$incremental_results$inc.cost / n.ais,
  inc.qol  = base_case$incremental_results$inc.qol / n.ais,
  NMB      = base_case$incremental_results$NMB / n.ais
)
sc9_inc$scenario <- "sc9"

# also create a placeholder process_results row for consistency
sc9_proc <- copy(base_case$process_results)
sc9_proc$scenario <- "sc9"

#### ======================================= ####
#### COMBINE RESULTS                  ####
#### ======================================= ####

incremental.results <- rbind(base_case$incremental_results, sc1$incremental_results,
                             sc2$incremental_results, sc3$incremental_results,
                             sc4$incremental_results, sc5$incremental_results,
                             sc6$incremental_results,
                             sc7$incremental_results,
                             sc8$incremental_results,
                             sc9_inc)

process.results <- rbind(base_case$process_results, sc1$process_results,
                         sc2$process_results, sc3$process_results,
                         sc4$process_results, sc5$process_results,
                         sc6$process_results, sc7$process_results,
                         sc8$process_results)

### edit dps
incremental.results <- as.data.table(incremental.results)
incremental.results[, inc.cost := comma(round(inc.cost, 0))]
incremental.results[, inc.qol  := comma(round(inc.qol, 2))]
incremental.results[, NMB      := comma(round(NMB, 0))]

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
write.csv(incremental.results, "outputs/scenario_incremental_results.csv", row.names = FALSE)
write.csv(process.results, "outputs/scenario_process_results.csv", row.names = FALSE)

#### ======================================= ####
#### PSA SUMMARY: Mean and 95% CI (percentile method) ####
#### ======================================= ####

psa <- read.csv("outputs/psa_outputs.csv")
psa <- as.data.table(psa)

psa_summary <- data.table(
  metric = c("inc.cost", "inc.qol", "NMB"),
  mean   = c(mean(psa$inc.cost), mean(psa$inc.qol), mean(psa$NMB)),
  lower  = c(quantile(psa$inc.cost, 0.025), quantile(psa$inc.qol, 0.025), quantile(psa$NMB, 0.025)),
  upper  = c(quantile(psa$inc.cost, 0.975), quantile(psa$inc.qol, 0.975), quantile(psa$NMB, 0.975))
)

# Also calculate probability of cost-effectiveness at WTP = 20,000
pCE <- mean(psa$NMB > 0)
cat(sprintf("Probability cost-effective at £20,000/QALY: %.1f%%\n", pCE * 100))

write.csv(psa_summary, "outputs/psa_summary.csv", row.names = FALSE)
