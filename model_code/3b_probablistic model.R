###############################################
# TITLE: Running the Probablistic Model Runs for B360S Model
# AUTHOR: Nichola Naylor (OI Pharma Partners Ltd), aided by GPT-4o,GTP-5 & Github co-pilot
# DATE: September 2025
#
# DESCRIPTION:
# Runs the full decision-tree + mRS Markov model under probabilistic sampling 
# to characterize parameter uncertainty (PSA). Generates distributions for 
# total & incremental costs/QALYs, ICER, NMB, CE plane points, and CEAC. It:
# 1) Draws parameter sets (size = n.sim) from predefined distributions.
# 2) For each draw, runs `run_model()` and stores key outcomes.
# 3) Aggregates PSA results (means, medians, CrIs) and formats for plots/tables.
# 4) Saves tidy outputs for downstream figures: CE plane, CEAC, and tornado/EVPI.
#
# INPUTS:
# - inputs/created_inputs/parameters_edited.RData   # data_main (with PSA columns)
# - inputs/created_inputs/mrs_samples_mean.RData    # mRS-level cost/utility/mortality
# - nputs/psa_draws.RData       # pre-sampled draws 
#
# OUTPUTS:
# - outputs/psa_results.RData (.csv): per-draw costs, QALYs, inc.cost, inc.QALY, NMB
# - outputs/ce_plane.csv: x=ΔQALY, y=ΔCost, label by strategy
# - outputs/ceac.csv: WTP vs Pr(NMB>0)
# - outputs/summary_psa.csv: mean/SD/quantiles for costs, QALYs, ICER, NMB
#
# DEPENDENCIES:
# - 2_probabilistic_sampling.R (defines distributions / draws)   [if used]
# - core model functions: mrs_markov(), run_model()
# - Packages: data.table, dplyr, assertthat, conflicted
###############################################

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

source("model_code/model_functions.R")

### load data
load("inputs/created_inputs/parameters_edited.RData")
data_main <- parameters
load("inputs/created_inputs/params.psa.sample.RData")
load("inputs/created_inputs/mrs_samples.RData")
load("inputs/created_inputs/mrs_samples_mean.RData")

#### ======================================= ####
####       3. FORMAT DATA                   ####
#### ======================================= ####

n.sample <- 1000

###### for params.psa 
### create a new column in params.psa.sample to indicate the sample set
params.psa <- data_main[rep(data_main[, .I], n.sample )] ## create 1000 samples of each row
params.psa$sample_id <- rep(1:n.sample, each=nrow(data_main))

# Reshape sample_set to long format (for easier joining)
sample_set_long <- melt(params.psa.sample, id.vars = "sample_id", 
                        variable.name = "variable", 
                        value.name = "sample_value")

### edit variable to split for commas
# Convert variable names to character
sample_set_long[, variable := as.character(variable)]

# Split into model_param, Presentation.Setting, and Intervention
sample_set_long[, c("model_param", "Presentation.Setting", "Intervention") := {
  split_vars <- strsplit(variable, ",")
  list(
    sapply(split_vars, function(x) x[1]),                              # model_param
    sapply(split_vars, function(x) if (length(x) > 1) x[2] else NA_character_),  # Presentation.Setting
    sapply(split_vars, function(x) if (length(x) > 2) x[3] else NA_character_)   # Intervention
  )
}]

sample_set_long[ , PSAflag := 1]

# Merge all in one step now that columns are structured properly
params.psa <- merge(params.psa,
                    sample_set_long,
                    by = c("sample_id", "model_param", "Presentation.Setting", "Intervention"),
                    all.x = TRUE)

# Update base_case where sampled value exists
params.psa[!is.na(sample_value), base_case := sample_value]

# Clean up
params.psa[, sample_value := NULL]
params.psa[, variable := NULL]
rm(params.psa.sample, sample_set_long)

###### for mrs.sample
mrs_samples <- as.data.table(mrs_samples)
mrs.samples.list <- split(mrs_samples, by = "sample_id", keep.by = FALSE)

### Run model on samples
# Pre-split params.psa by sample_id
params.psa.list <- split(params.psa, by = "sample_id", keep.by = FALSE)

#### ======================================= ####
####     QUICK CHECK: Updated PSA Parameters  ####
#### ======================================= ####

cat("\n--- Sanity check: sampled parameter updates ---\n")

# Pick a few key parameters expected to vary
params_to_check <- c("p.ivt.emt2mt", "p.noivt.emt2mt", "c.mt", "utility.mrs")

# Check that each appears across multiple sample_ids and differs
for (param in params_to_check) {
  cat("\nParameter:", param, "\n")
  print(params.psa[model_param == param,
                   .(sample_id, Presentation.Setting, Intervention, base_case)][1:5, ])
}

# Confirm variation — if variance = 0, something didn’t update or non PSA value
var_summary <- params.psa[, .(variance = var(base_case, na.rm = TRUE)), by = model_param]
cat("\nParameters with zero variance (potential issue):\n")
print(var_summary[variance == 0])

# Overall structural sanity
cat("\nTotal rows:", nrow(params.psa), 
    "\nDistinct sample_ids:", length(unique(params.psa$sample_id)),
    "\nAny missing base_case values?:", sum(is.na(params.psa$base_case)), "\n")

########## RUN MODEL ON SAMPLES ################
psa.outputs <- vector("list", n.sample)
pb <- txtProgressBar(min = 1, max = n.sample, initial = 0, style = 3)

for (i in 1:n.sample) {
  setTxtProgressBar(pb, i)
  temp <- params.psa.list[[i]]
  temp.mrs <- mrs.samples.list[[i]]
  
  temp.outputs <- try(run_model(temp, 10, temp.mrs), silent = TRUE)
  
  ## set up to still run if there is an error on 1 sample
  if (inherits(temp.outputs, "try-error")) {
    warning(sprintf("Model run failed at sample %d", i))
    errors[[i]] <- conditionMessage(attr(temp.outputs, "condition")) 
    next
  }
  psa.outputs[[i]] <- c(temp.outputs[[2]], temp.outputs[[3]]) ### !!! note list order important here for pull through
}
psa.outputs.dt <- rbindlist(psa.outputs, use.names=TRUE, fill=TRUE, idcol=TRUE)
write.csv(psa.outputs.dt, file="outputs/psa_outputs.csv")


