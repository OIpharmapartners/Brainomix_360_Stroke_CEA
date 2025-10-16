###############################################
# TITLE: Hospital Data Processing for B360S Model
# AUTHOR: Nichola Naylor (OI Pharma Partners Ltd), aided by GPT-4o,GTP-5 & Github co-pilot
# DATE: September 2025
#
# DESCRIPTION:
# This script processes regional hospital data to:
# - Identify CSC-flagged hospitals and assign ISDN numbers
# - Clean and filter stroke and IVT data
# - Derive IVT for those eligible by hospital type (CSC vs non-CSC)
# - Update parameter values used in further modeling R scripts
#
# INPUTS: 
# - inputs/hospital_names_old.csv
# - inputs/csc_key.csv
# - inputs/isdn_key.csv
# - inputs/isdn_ivt.csv
# - inputs/parameters.csv
#
# OUTPUTS:
# - inputs/created_inputs/IO_isdn_hosp_numbers.csv
# - inputs/created_inputs/parameters_edited.RData
# 
# Notes: Ensure all sourced parameters that are used are available and are in the same
# format as assumed in this script/as in the previous data.
###############################################


####### 1. INITIALISE & LOAD LIBRARIES  #######
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

#######  2. LOAD DATA #####

parameters <- read.csv("inputs/parameters.csv")
parameters <- as.data.table(parameters)

### read in hospital data from previous model
hospital_l <- read.csv("inputs/hospital_names_old.csv") ## total list of hospital names

csc_l <- read.csv("inputs/csc_key.csv") ## list of hospitals which are csc
csc_l$cscflag <- 1

## duplicate checks
stopifnot(!any(duplicated(hospital_l$Hospital)))
stopifnot(!any(duplicated(csc_l$Hospital)))

hospital_l <- merge(hospital_l, csc_l, by="Hospital", all=TRUE)
hospital_l <- as.data.table(hospital_l)

hospital_l[is.na(cscflag), cscflag := 0]

####### 3. Merge WITH ISDN DATASET   #####

isdn <- read.csv("inputs/isdn_key.csv")

#### keep latest year
isdn <- as.data.table(isdn)

stopifnot("Apr 2022-Mar 2023" %in% isdn$date) ### !!! different years will need different filters
isdn <- isdn[date=="Apr 2022-Mar 2023"] ### !!! different years will need different filters

stopifnot(!any(duplicated(isdn$Hospital)))

isdn_h <- merge(hospital_l, isdn, by="Hospital", all=TRUE)

isdn_h <- isdn_h  %>% filter(!is.na(ISDN)) %>%  ### REMOVE THOSE WITH NO IDSN or stroke numbers
  filter(!is.na(n.stroke) & n.stroke!="Too few to report" & n.stroke!=".") %>% 
  filter(!stringr::str_detect(Hospital, 'Rehabilitation|Rehab') ) %>% ### REMOVE REHAB UNITS/CENTRES
  filter(ISDN!="Wales" & ISDN!="Northern Ireland") %>% ### REMOVE WALES & NI
  mutate(cscflag=if_else(str_detect(Hospital,'HASU'), true = 1, false = cscflag))%>%
  as.data.table()

### CHANGE CSC FLAGS FOR THOSE WITHOUT 
##!!! if ISDNs/hospitals update this will need checking
isdn_h[is.na(cscflag), cscflag:=0]
isdn_h[Hospital=="Queen's Medical Centre - Nottingham",cscflag:=1]

length(which(isdn_h$cscflag==1))
length(which(isdn_h$cscflag==0))

### some checks for regional analysis
if (anyNA(isdn_h$ISDN) || anyNA(isdn_h$n.stroke))
  stop("Missing ISDN or n.stroke after merges. Check keys and filters before proceeding.")
isdn_h[, n.stroke := suppressWarnings(as.numeric(n.stroke))]
if (any(is.na(isdn_h$n.stroke))) stop("n.stroke has non-numeric entries after cleaning.")


write.csv(isdn_h, file="inputs/created_inputs/IO_isdn_hosp_numbers.csv", row.names = FALSE)

# ##### finding out p.asc
# isdn_h[ , n.stroke := as.numeric(n.stroke)]
# isdn_h%>%
#   group_by(cscflag) %>%
#   summarise(
#     total_value = sum(n.stroke)) %>%
#   as.data.table()

###### 4. IVT ELIGIBILITY PROCESSING   ######

## read in isdn data of ivt out of those eligible
isdn <- read.csv("inputs/isdn_ivt.csv") 

isdn <- as.data.table(isdn)

isdn_ht <- merge(hospital_l, isdn, by="Hospital", all=TRUE)

isdn_ht <- isdn_ht  %>% filter(!is.na(ISDN)) %>%  ### REMOVE THOSE WITH NO IDSN or stroke numbers
  filter(!is.na(p_eivt2ivt) & p_eivt2ivt!="Too few to report" & p_eivt2ivt!=".") %>% 
  filter(!stringr::str_detect(Hospital, 'Rehabilitation|Rehab') ) %>% ### REMOVE REHAB UNITS/CENTRES
  filter(ISDN!="Wales" & ISDN!="Northern Ireland") %>% ### REMOVE WALES & NI
  mutate(cscflag=if_else(str_detect(Hospital,'HASU'), true = 1, false = cscflag))%>%
  as.data.table()

### CHANGE CSC FLAGS FOR THOSE WITHOUT 
##!!! if ISDNs/hospitals update this will need checking
isdn_ht[is.na(cscflag), cscflag:=0]
isdn_ht[Hospital=="Queen's Medical Centre - Nottingham",cscflag:=1]

isdn_ht[ , p_eivt2ivt := as.numeric(p_eivt2ivt)] ## proportion of those eligible for ivt that get ivt
x <- isdn_ht[ ,.(v1average=mean(p_eivt2ivt, na.rm=TRUE)),by=c("cscflag","date")]

### create summary data to parameterise ivt transitions
summary_stats <- x %>%
  group_by(cscflag) %>%
  summarise(
    median_value = median(v1average),
    min_value = min(v1average),
    max_value = max(v1average)) %>%
  as.data.table()

summary_stats

### CREATES EDITED PARAMETER DATA FRAME TO USE
## early ASC
## divided by 100 to make proportion per year
#### !!! note to user: check the inputs are in the correct format (i.e. %s)
stopifnot(max(isdn_ht$p_eivt2ivt, na.rm=TRUE) > 1)  # expect percents, not proportions

parameters[Presentation.Setting == "early;ASC" & 
             model_param == "p.eivt2ivt", 
           base_case := summary_stats[cscflag == 0, median_value]/100]

parameters[Presentation.Setting == "early;ASC" & 
             model_param == "p.eivt2ivt", 
           PSA_low := summary_stats[cscflag == 0, min_value]/100]

parameters[Presentation.Setting == "early;ASC" & 
             model_param == "p.eivt2ivt", 
           PSA_high := summary_stats[cscflag == 0, max_value]/100]


parameters[Presentation.Setting == "early;CSC" & 
             model_param == "p.eivt2ivt", 
           base_case := summary_stats[cscflag == 1, median_value]/100]

parameters[Presentation.Setting == "early;CSC" & 
             model_param == "p.eivt2ivt", 
           PSA_low := summary_stats[cscflag == 1, min_value]/100]

parameters[Presentation.Setting == "early;CSC" & 
             model_param == "p.eivt2ivt", 
           PSA_high := summary_stats[cscflag == 1, max_value]/100]

save(parameters,file="inputs/created_inputs/parameters_edited.RData")



