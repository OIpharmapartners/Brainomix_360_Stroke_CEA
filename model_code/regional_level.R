###############################################
# TITLE: ISDN level analysis for the B360S CEA model
# AUTHOR: Nichola Naylor (OI Pharma Partners Ltd), aided by GPT-4o,GTP-5 & Github co-pilot
# DATE: September 2025
#
# DESCRIPTION:
# Adapts the core CEA to a single ISDN region.
#
# WHAT THIS SCRIPT DOES:
# 1) Loads/merges local stroke numbers
# 2) Re-parameterizes `data_main` for the local setting (overrides national means).
# 3) Runs `run_model()` for base
# 4) Produces local outputs
#
#
# INPUTS:
# - inputs/created_inputs/parameters_edited.RData       # national baseline
# - inputs/created_inputs/mrs_samples_mean.RData  # national baseline
# - inputs/created_inputs/IO_isdn_hosp_numbers.csv  # hospital numbers of stroke patients
#
# OUTPUTS:
# - outputs/ISDN_results.RData
# - outputs/ISDN_incremental_results.csv
#
# DEPENDENCIES:
# - core model functions: mrs_markov(), run_model()
# - Packages: see code
######################################################

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

### read in hospital data
hospital_l <- read.csv("inputs/created_inputs/IO_isdn_hosp_numbers.csv")

hospital_l <- as.data.table(hospital_l)

hospital_l[is.na(cscflag), cscflag := 0]


hospital_l <- hospital_l[!is.na(cscflag)&!is.na(ISDN)]
## none should drop

network_results <- list() ## create to put outputs in

hospital_l[ , group_ID := .GRP, by =.(ISDN)]

stopifnot(all(c("cscflag", "ISDN", "n.stroke") %in% names(hospital_l)))

for (i in 1:max(unique(hospital_l$group_ID))){

  
  x <- hospital_l[group_ID==i]
  
  n <- sum(x$n.stroke)
  
  ## might want to add split here of csc flag & impact in pathway in different
  ## iterations !!
  
  ISDN <- x$ISDN[1]
  
  temp_data_main <- data.table::copy(data_main)
  temp_data_main[model_param=="n.asc", base_case := length(which(x$cscflag == 0))]
  temp_data_main[model_param=="n.csc", base_case := length(which(x$cscflag == 1))]
  temp_data_main[model_param=="n.stroke", base_case := n]
    
    y <- run_model(temp_data_main, cycles = 10, mrs_samples_mean = mrs_samples_mean)

    
    network_results[[ISDN]] <- y
 
}

save(network_results,file="outputs/network_results.RData")

## save table of incremental results

incremental_df <- do.call(rbind, lapply(names(network_results), function(region) {
  result <- network_results[[region]]$incremental_results
  if (is.data.frame(result)) {
    result$region <- region
    return(result)
  } else {
    return(NULL)
  }
}))

## formatting
incremental <- as.data.table(incremental_df)
incremental[, inc.cost := comma(round(inc.cost, 0))]
incremental[, inc.qol := comma(round(inc.qol, 2))]
incremental[, NMB := comma(round(NMB, 0))]

write.csv(incremental, file = "outputs/ISDN_incremental_results.csv", row.names = FALSE)
