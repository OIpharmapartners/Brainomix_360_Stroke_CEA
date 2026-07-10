###############################################
# TITLE: Running the Probablistic Model Runs for B360S Model
# AUTHOR: Nichola Naylor (OI Pharma Partners Ltd), aided by GPT-4o,GTP-5 & Github co-pilot, Claude Opus 4.6 and Claude Opus 4.8
#
# DESCRIPTION:
# Runs the full decision-tree + mRS Markov model under probabilistic sampling 
# to characterize parameter uncertainty (PSA). 
#
# INPUTS
# - inputs/created_inputs/parameters_edited.RData
#     Contains data_main: model parameters after preprocessing and derived PSA updates.
#
# - inputs/created_inputs/mrs_samples.RData
#     Contains mrs_samples: PSA samples for mRS-specific costs, utilities,
#     mortality relative risks, and mRS treatment distributions.

# Created files:
# - outputs/psa_outputs.csv
#     Combined probabilistic model outputs across PSA samples.
#     Includes total and incremental costs, QALYs, and NMB components
#     returned by run_model().
#################################################

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
library(scales)


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
load("inputs/created_inputs/params_psa_sample.RData")
load("inputs/created_inputs/mrs_samples.RData")
load("inputs/created_inputs/mrs_samples_mean.RData")
load("inputs/created_inputs/dist_sample_df.RData")


#### ======================================= ####
####       3. FORMAT DATA                   ####
#### ======================================= ####

n.sample <- length(unique(mrs_samples$sample_id))

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

### for distribution samples
params.psa <- merge(
  params.psa,
  as.data.table(dist_sample_df)[, .(sample_id, model_param, mrs, dir_value = dist_sample)],
  by = c("sample_id", "model_param", "mrs"), all.x = TRUE)
params.psa[!is.na(dir_value), base_case := dir_value]
params.psa[, dir_value := NULL]

###### for mrs.sample
mrs_samples <- as.data.table(mrs_samples)
mrs.samples.list <- split(mrs_samples, by = "sample_id", keep.by = FALSE)

### Run model on samples
# Pre-split params.psa by sample_id
params.psa.list <- split(params.psa, by = "sample_id", keep.by = FALSE)

####   CHECK: Updated PSA Parameters  ####


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

#### Dirichlet mRS distributions: sampled + integrated correctly ####
dist_params <- c("dist.ivt", "dist.noivt", "dist.mt", "dist.nomt")

# (a) integration intact: every draw still sums to 1 across the 7 mRS rows
dist_sums <- params.psa[model_param %in% dist_params,
                        .(s = sum(base_case)), by = .(model_param, sample_id)]
stopifnot(dist_sums[, all(abs(s - 1) < 1e-8)])

# (b) actually sampled: base_case varies across draws (not left at the fixed base case)
dist_var <- params.psa[model_param %in% dist_params,
                       .(variance = var(base_case)), by = .(model_param, mrs)]
stopifnot(dist_var[, all(variance > 0)])

cat("\nDirichlet mRS distributions: all draws sum to 1 and vary across samples - OK\n")

########## RUN MODEL ON SAMPLES ################
psa.outputs <- vector("list", n.sample)
pb <- txtProgressBar(min = 1, max = n.sample, initial = 0, style = 3)
errors <- list() # set up errors list in case need to store errors

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

assert_that(
  nrow(psa.outputs.dt) == n.sample,
  msg = sprintf(
    "PSA output row count mismatch: expected %s successful samples, got %s. Failed samples: %s",
    n.sample,
    nrow(psa.outputs.dt),
    paste(which(sapply(errors, Negate(is.null))), collapse = ", ")
  )
)


write.csv(psa.outputs.dt, file="outputs/psa_outputs.csv")


