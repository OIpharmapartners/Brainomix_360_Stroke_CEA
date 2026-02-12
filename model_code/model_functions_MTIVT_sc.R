###############################################
# TITLE: Core Model Function for B360S Model - FOR ALTERNATIVE MT + IVT IMPACT SCENARIO
# AUTHOR: Nichola Naylor (OI Pharma Partners Ltd), aided by GPT-4o,GTP-5 & Github co-pilot
# DATE: September 2025
#
# DESCRIPTION: SEE "model_functions.R" for further descriptions
#
# !!! IMPORTANT NOTE ON CYCLE INTERPRETATION
#
# - The Markov model (mrs_markov function) uses cycles that
#   represent annual transitions (i.e., 1 cycle = 1 year).
#   This aligns with standard practice in health economics
#   for long-term cost and QALY modelling by mRS state.
#
# - In contrast, the decision-tree model (run_model function)
#   uses cycles purely as *transition steps* to clear the cohort
#   through complex pathways (e.g., stroke treatment sequences).
#   These cycles do NOT correspond to actual time units.
#
# Ensure this distinction is maintained when interpreting
# outputs or modifying model logic.
############################################################



###############################################

#### ======================================= ####
####       1. INITIALISE & LOAD LIBRARIES      ####
#### ======================================= ####


# Resolve potential function conflicts (tidyverse vs data.table)
conflict_prefer("select", "dplyr")       # dplyr::select over MASS::select
conflict_prefer("data.table", "data.table")  # Ensure `[.data.table` overrides

# source("model_code/model_functions.R") ## !!! need to source mrs_markov if running standalone

#### ======================================= ####
####       4. COMPLETE MODEL RUN FUNCTION    ####
#### ======================================= ####

# ### testing code (and then # out run_model to run through)
# load("inputs/created_inputs/parameters_edited.RData")
# data_main <- parameters
# load("inputs/created_inputs/mrs_samples_mean.RData")
# cycles <- 10


#' Title: Model run
#'
#' This function executes the full decision-analytic model for one set of
#' input parameters. It combines an acute-phase decision tree (for treatment
#' eligibility, reperfusion, and AI effects) with a long-term Markov model
#' (mrs_markov) that tracks health-state transitions by modified Rankin Scale (mRS).
#'
#' @param data_main_inp A data.table containing all model parameters.
#' @param cycles Integer. Number of Markov cycles to simulate (default = 10).
#'   Each cycle represents one year of follow-up.
#' @param mrs_samples_mean_inp A data.table containing mean mRS-specific
#'   costs, utilities, and mortality rates
#'
#' Internal validation checks ensure:
#' \itemize{
#'   \item Transition matrix rows sum to 1 (within tolerance).
#'   \item All probabilities are within [0,1].
#'   \item No negative or missing state counts.
#'   \item Total stroke cohort (\code{n.stroke}) is conserved.
#' }
#'
#' @return A named \code{list} with the following elements:
#' \describe{
#'   \item{\code{process_results}}{Detailed process-level outcomes and costs (e.g., imaging and reperfusion).}
#'   \item{\code{summary_data_all}}{Summary of population-level costs, QALYs, and totals for each arm.}
#'   \item{\code{incremental_results}}{Incremental costs, QALYs, and Net Monetary Benefit (NMB) between arms.}
#'   \item{\code{mrs.trace}}{Markov trace showing cohort distribution by mRS over cycles.}
#'   \item{\code{trace.standard}}{Acute-phase state trace for the standard-care arm.}
#'   \item{\code{trace.intervention}}{Acute-phase state trace for the AI intervention arm.}
#'   \item{\code{tm.standard}}{Transition matrices (3D array) for the standard-care arm.}
#'   \item{\code{tm.intervention}}{Transition matrices (3D array) for the AI intervention arm.}
#' }
#'
#' @export
run_model_mtivtsc <- function(data_main_inp=data_main, cycles=10,
                              mrs_samples_mean_inp=mrs_samples_mean){
  ##defensive copies to reduce the likelihood of unintended alterations
  data_main       <- data.table::copy(data_main_inp)
  mrs_samples_mean <- data.table::copy(mrs_samples_mean_inp)
  
  tolerance <- 1e-8 ## tolerance for rounding errors in checks
  
  #### ======================================= ####
  ####       4a. SET UP OF DECISION TREE      ####
  #### ======================================= ####
  
  state.names <- c("stroke"       
                   ,"AIS"     
                   , "no_AIS"
                   ,"EARLY"  ## these numbers are used for NCCT + CTA/CTP/MRI splits later on
                   ## !!! note if p.NCCT is not = 1, need to add further splits here or later
                   ## when calculating testing cost
                   ,"LATE"   ## these numbers are used for NCCT + CTA/CTP/MRI splits later on
                   ,"ASC_EARLY"  
                   ,"CSC_EARLY"   
                   ,"EIVT_EARLY_ASC"
                   ,"NOEIVT_EARLY_ASC"
                   ,"EIVT_EARLY_CSC"
                   ,"NOEIVT_EARLY_CSC"
                   ,"IVT_EARLY_ASC"
                   ,"NOIVT_EARLY_ASC" 
                   ,"IVT_EARLY_CSC" 
                   ,"NOIVT_EARLY_CSC" 
                   ,"EMT_EARLY_ASC_IVT" 
                   ,"EMT_EARLY_ASC_NOIVT"
                   ,"EMT_EARLY_CSC_IVT"
                   ,"EMT_EARLY_CSC_NOIVT"
                   ,"NOEMT_EARLY_ASC_IVT"
                   ,"NOEMT_EARLY_ASC_NOIVT"
                   ,"NOEMT_EARLY_CSC_IVT"
                   ,"NOEMT_EARLY_CSC_NOIVT"
                   ,"ASC_LATE"
                   ,"CSC_LATE" 
                   ,"EMT_LATE_ASC" 
                   ,"NOEMT_LATE_ASC"
                   ,"EMT_LATE_CSC"
                   ,"NOEMT_LATE_CSC"
                   ,"MT"
                   ,"NOMT")    
  
  n.states <- length(state.names)
  n.stroke <- data_main[model_param=="n.stroke",base_case]
  assert_that(
    length(n.stroke) == 1,
    !is.na(n.stroke),
    is.numeric(n.stroke),
    n.stroke >= 0,
    abs(n.stroke - round(n.stroke)) < .Machine$double.eps^0.5,
    msg = sprintf("`n.stroke` must be a non-negative integer. Got: %s", n.stroke)
  )
  #  Seed the starting states of the model
  seed <- c(n.stroke, rep(0,(n.states-1))) ## all people start in the first state
  
  #  We start with a three dimensional array in order to capture the time dependencies
  tm.standard <- array(data=0,dim=c(n.states, n.states, cycles),
                       dimnames= list(state.names, state.names, 1:cycles)) 
  ## an empty array of dimensions (number of states, number of states, number of cycles)
  ### create a loop that creates a time dependent transition matrix for each cycle
  
  #### ======================================= ####
  ####       4b. NO INTERVENTION ARM            ####
  #### ======================================= ####
  
  no_int <- data_main[Intervention=="all"|Intervention=="no intervention"]
  
  for (i in 1:cycles) {
    
    ## transitions out of AIS & NCCT
    tm.standard["stroke","AIS",i] <- no_int[model_param=="p.ais",base_case]
    tm.standard["stroke","no_AIS",i] <-1- no_int[model_param=="p.ais",base_case]
    
    tm.standard["AIS","EARLY",i] <- no_int[model_param=="p.early",base_case] 
    tm.standard["AIS","LATE",i] <- 1-no_int[model_param=="p.early",base_case]
    
    ## early pathway:
    tm.standard["EARLY","ASC_EARLY",i] <- no_int[model_param=="p.asc",base_case]
    tm.standard["EARLY","CSC_EARLY",i] <- 1-(no_int[model_param=="p.asc",base_case])
    
    ## early:ASC
    tm.standard["ASC_EARLY","EIVT_EARLY_ASC",i] <- no_int[model_param=="p.eivt"& 
                                                            Presentation.Setting=="early",base_case]
    tm.standard["ASC_EARLY","NOEIVT_EARLY_ASC",i] <- 1-(no_int[model_param=="p.eivt"& 
                                                                 Presentation.Setting=="early",base_case])
    
    tm.standard["EIVT_EARLY_ASC","IVT_EARLY_ASC",i] <- no_int[model_param=="p.eivt2ivt"& 
                                                                Presentation.Setting=="early;ASC",base_case]
    tm.standard["EIVT_EARLY_ASC","NOIVT_EARLY_ASC",i] <- 1-(no_int[model_param=="p.eivt2ivt"& 
                                                                     Presentation.Setting=="early;ASC",base_case])
    
    tm.standard["NOEIVT_EARLY_ASC","NOIVT_EARLY_ASC",i] <- 1
    
    tm.standard["IVT_EARLY_ASC","EMT_EARLY_ASC_IVT",i] <- no_int[model_param=="p.ivt2emt"& 
                                                                   Presentation.Setting=="early",base_case]*
      no_int[model_param=="p.plao.nihss",base_case]
    
    tm.standard["IVT_EARLY_ASC","NOEMT_EARLY_ASC_IVT",i] <- 1- tm.standard["IVT_EARLY_ASC","EMT_EARLY_ASC_IVT",i]
    
    tm.standard["NOIVT_EARLY_ASC","EMT_EARLY_ASC_NOIVT",i] <- no_int[model_param=="p.noivt2emt"& 
                                                                       Presentation.Setting=="early",base_case]*
      no_int[model_param=="p.plao.nihss",base_case]
    
    tm.standard["NOIVT_EARLY_ASC","NOEMT_EARLY_ASC_NOIVT",i] <- 1-  tm.standard["NOIVT_EARLY_ASC","EMT_EARLY_ASC_NOIVT",i] 
    
    tm.standard["EMT_EARLY_ASC_IVT","MT",i] <- no_int[model_param=="p.ivt.emt2mt"& 
                                                        Presentation.Setting=="early;ASC",base_case]
    tm.standard["EMT_EARLY_ASC_IVT","NOMT",i] <- 1-(no_int[model_param=="p.ivt.emt2mt"& 
                                                             Presentation.Setting=="early;ASC",base_case])
    tm.standard["NOEMT_EARLY_ASC_IVT","NOMT",i] <- 1
    
    tm.standard["EMT_EARLY_ASC_NOIVT","MT",i] <- no_int[model_param=="p.noivt.emt2mt"& 
                                                          Presentation.Setting=="early;ASC",base_case]
    tm.standard["EMT_EARLY_ASC_NOIVT","NOMT",i] <- 1-(no_int[model_param=="p.noivt.emt2mt"& 
                                                               Presentation.Setting=="early;ASC",base_case])
    tm.standard["NOEMT_EARLY_ASC_NOIVT","NOMT",i] <- 1
    
    ## early CSC pathway
    tm.standard["CSC_EARLY","EIVT_EARLY_CSC",i] <- no_int[model_param=="p.eivt"& 
                                                            Presentation.Setting=="early",base_case]
    tm.standard["CSC_EARLY","NOEIVT_EARLY_CSC",i] <- 1-(no_int[model_param=="p.eivt"& 
                                                                 Presentation.Setting=="early",base_case])
    
    tm.standard["EIVT_EARLY_CSC","IVT_EARLY_CSC",i] <- no_int[model_param=="p.eivt2ivt"& 
                                                                Presentation.Setting=="early;CSC",base_case]
    tm.standard["EIVT_EARLY_CSC","NOIVT_EARLY_CSC",i] <- 1-(no_int[model_param=="p.eivt2ivt"& 
                                                                     Presentation.Setting=="early;CSC",base_case])
    
    tm.standard["NOEIVT_EARLY_CSC","NOIVT_EARLY_CSC",i] <- 1
    
    tm.standard["IVT_EARLY_CSC","EMT_EARLY_CSC_IVT",i] <- no_int[model_param=="p.ivt2emt"& 
                                                                   Presentation.Setting=="early",base_case]*
      no_int[model_param=="p.plao.nihss",base_case]
    tm.standard["IVT_EARLY_CSC","NOEMT_EARLY_CSC_IVT",i] <- 1-  tm.standard["IVT_EARLY_CSC","EMT_EARLY_CSC_IVT",i]
    
    tm.standard["NOIVT_EARLY_CSC","EMT_EARLY_CSC_NOIVT",i] <- no_int[model_param=="p.noivt2emt"& 
                                                                       Presentation.Setting=="early",base_case]*
      no_int[model_param=="p.plao.nihss",base_case]
    tm.standard["NOIVT_EARLY_CSC","NOEMT_EARLY_CSC_NOIVT",i] <- 1-  tm.standard["NOIVT_EARLY_CSC","EMT_EARLY_CSC_NOIVT",i]
    
    tm.standard["EMT_EARLY_CSC_IVT","MT",i] <- no_int[model_param=="p.ivt.emt2mt"& 
                                                        Presentation.Setting=="early;CSC",base_case]
    tm.standard["EMT_EARLY_CSC_IVT","NOMT",i] <- 1-(no_int[model_param=="p.ivt.emt2mt"& 
                                                             Presentation.Setting=="early;CSC",base_case])
    tm.standard["NOEMT_EARLY_CSC_IVT","NOMT",i] <- 1
    
    tm.standard["EMT_EARLY_CSC_NOIVT","MT",i] <- no_int[model_param=="p.noivt.emt2mt"& 
                                                          Presentation.Setting=="early;CSC",base_case]
    tm.standard["EMT_EARLY_CSC_NOIVT","NOMT",i] <- 1-(no_int[model_param=="p.noivt.emt2mt"& 
                                                               Presentation.Setting=="early;CSC",base_case])
    tm.standard["NOEMT_EARLY_CSC_NOIVT","NOMT",i] <- 1
    
    ## late 
    tm.standard["LATE","ASC_LATE",i] <- no_int[model_param=="p.asc",base_case]  
    tm.standard["LATE","CSC_LATE",i] <- 1-(no_int[model_param=="p.asc",base_case])
    
    ## late ASC pathway
    tm.standard["ASC_LATE","EMT_LATE_ASC",i] <- no_int[model_param=="p.emt" & 
                                                         Presentation.Setting=="late;ASC",base_case]*
      no_int[model_param=="p.plao.nihss",base_case]
    tm.standard["ASC_LATE","NOEMT_LATE_ASC",i] <- 1-  tm.standard["ASC_LATE","EMT_LATE_ASC",i]
    
    tm.standard["EMT_LATE_ASC","MT",i] <- no_int[model_param=="p.emt2mt" & 
                                                   Presentation.Setting=="late;ASC",base_case]
    
    tm.standard["EMT_LATE_ASC","NOMT",i] <- 1-no_int[model_param=="p.emt2mt" & 
                                                       Presentation.Setting=="late;ASC",base_case]
    
    tm.standard["NOEMT_LATE_ASC","NOMT",i] <- 1
    
    ## lates CSC pathway
    tm.standard["CSC_LATE","EMT_LATE_CSC",i] <- no_int[model_param=="p.emt" & 
                                                         Presentation.Setting=="late;CSC",base_case]*
      no_int[model_param=="p.plao.nihss",base_case]
    tm.standard["CSC_LATE","NOEMT_LATE_CSC",i] <- 1-   tm.standard["CSC_LATE","EMT_LATE_CSC",i] 
    
    tm.standard["EMT_LATE_CSC","MT",i] <- no_int[model_param=="p.emt2mt" & 
                                                   Presentation.Setting=="late;CSC",base_case]
    
    tm.standard["EMT_LATE_CSC","NOMT",i] <- 1-no_int[model_param=="p.emt2mt" & 
                                                       Presentation.Setting=="late;CSC",base_case]
    
    tm.standard["NOEMT_LATE_CSC","NOMT",i] <- 1
    
    ### absorbing states
    tm.standard["no_AIS","no_AIS",i] <- 1
    tm.standard["MT","MT",i] <- 1
    tm.standard["NOMT","NOMT",i] <- 1
    
    
  }
  
  # Check row sums for each cycle
  ### will only show when debugging/working within function
  for (j in 1:cycles) {
    row_sums <- rowSums(tm.standard[,,j], na.rm = TRUE)
    
    if (!all(abs(row_sums - 1) < tolerance)) {
      # Find problematic rows (states)
      bad_rows <- which(abs(row_sums - 1) >= tolerance)
      bad_states <- rownames(tm.standard)[bad_rows]
      
      # Create error message
      error_msg <- paste("'No intervention' transition matrix row sums do not equal 1 at cycle", j,
                         "for states:", paste(bad_states, collapse = ", "))
      
      warning(error_msg)
    }
  }
  
  
  #  Create a trace matrix
  trace.standard <- matrix(data=0, nrow=cycles, ncol=n.states)
  colnames(trace.standard) <- state.names
  
  trace.standard[1,] <- seed%*%tm.standard[,,1]
  
  for (i in 2:cycles) {
    trace.standard[i,] <- trace.standard[i-1,]%*%tm.standard[,,i]
  }
  
  # trace.standard
  
  outputs.standard <- list(tm.standard,trace.standard)
  
  #### ======================================= ####
  ####       4c. INTERVENTION ARM               ####
  #### ======================================= ####
  tm.intervention <- tm.standard
  
  data_int <- data_main[Intervention=="all"|Intervention=="intervention"]
  
  ##### adapting impacted pathways
  for (i in 1:cycles){
    ## early:ASC - IVT
    tm.intervention["EIVT_EARLY_ASC","IVT_EARLY_ASC",i] <- data_int[model_param=="p.eivt2ivt"& 
                                                                      Presentation.Setting=="early;ASC",base_case]
    tm.intervention["EIVT_EARLY_ASC","NOIVT_EARLY_ASC",i] <- 1-(data_int[model_param=="p.eivt2ivt"& 
                                                                           Presentation.Setting=="early;ASC",base_case])
    ## early: CSC - IVT
    
    tm.intervention["EIVT_EARLY_CSC","IVT_EARLY_CSC",i]<- data_int[model_param=="p.eivt2ivt"& 
                                                                     Presentation.Setting=="early;CSC",base_case]
    tm.intervention["EIVT_EARLY_CSC","NOIVT_EARLY_CSC",i] <- 1-(data_int[model_param=="p.eivt2ivt"& 
                                                                           Presentation.Setting=="early;CSC",base_case])
    
    ## early; ASC - MT
    tm.intervention["EMT_EARLY_ASC_IVT","MT",i] <- data_int[model_param=="p.ivt.emt2mt"& 
                                                              Presentation.Setting=="early;ASC",base_case]
    tm.intervention["EMT_EARLY_ASC_IVT","NOMT",i] <- 1-(data_int[model_param=="p.ivt.emt2mt"& 
                                                                   Presentation.Setting=="early;ASC",base_case])
    
    tm.intervention["EMT_EARLY_ASC_NOIVT","MT",i] <- data_int[model_param=="p.noivt.emt2mt"& 
                                                                Presentation.Setting=="early;ASC",base_case]
    tm.intervention["EMT_EARLY_ASC_NOIVT","NOMT",i] <- 1-(data_int[model_param=="p.noivt.emt2mt"& 
                                                                     Presentation.Setting=="early;ASC",base_case])
    
    
    ### early; CSC - MT
    tm.intervention["EMT_EARLY_CSC_IVT","MT",i] <- data_int[model_param=="p.ivt.emt2mt"& 
                                                              Presentation.Setting=="early;CSC",base_case]
    tm.intervention["EMT_EARLY_CSC_IVT","NOMT",i] <- 1-(data_int[model_param=="p.ivt.emt2mt"& 
                                                                   Presentation.Setting=="early;CSC",base_case])
    
    tm.intervention["EMT_EARLY_CSC_NOIVT","MT",i] <- data_int[model_param=="p.noivt.emt2mt"& 
                                                                Presentation.Setting=="early;CSC",base_case]
    tm.intervention["EMT_EARLY_CSC_NOIVT","NOMT",i] <- 1-(data_int[model_param=="p.noivt.emt2mt"& 
                                                                     Presentation.Setting=="early;CSC",base_case])
    
    ### late; ASC - MT
    tm.intervention["EMT_LATE_ASC","MT",i] <- data_int[model_param=="p.emt2mt" & 
                                                         Presentation.Setting=="late;ASC",base_case]
    
    tm.intervention["EMT_LATE_ASC","NOMT",i] <- 1-data_int[model_param=="p.emt2mt" & 
                                                             Presentation.Setting=="late;ASC",base_case]
    
    ### late; CSC - MT 
    tm.intervention["EMT_LATE_CSC","MT",i] <- data_int[model_param=="p.emt2mt" & 
                                                         Presentation.Setting=="late;CSC",base_case]
    
    tm.intervention["EMT_LATE_CSC","NOMT",i] <- 1-data_int[model_param=="p.emt2mt" & 
                                                             Presentation.Setting=="late;CSC",base_case]
  }
  
  
  # ## check adds to 1
  # Check row sums for each cycle
  for (j in 1:cycles) {
    row_sums <- rowSums(tm.intervention[,,j], na.rm = TRUE)
    
    if (!all(abs(row_sums - 1) < tolerance)) {
      # Find problematic rows (states)
      bad_rows <- which(abs(row_sums - 1) >= tolerance)
      bad_states <- rownames(tm.intervention)[bad_rows]
      
      # Create error message
      error_msg <- paste("'Intervention' transition matrix row sums do not equal 1 at cycle", j,
                         "for states:", paste(bad_states, collapse = ", "))
      
      warning(error_msg)
      
      # Optional hard stop
      if (exists("stop_on_error") && stop_on_error) stop(error_msg)
    }
  }
  
  ### an additional check to catch smaller negative values outside bounds
  for (j in 1:cycles) {
    bad <- tm.standard[,,j] < -1e-12 | tm.standard[,,j] > 1 + 1e-12
    if (any(bad)) {
      idx <- which(bad, arr.ind = TRUE)[1,]
      stop(sprintf("Standard prob out of [0,1] at cycle %d: %s→%s = %.6f",
                   j, rownames(tm.standard)[idx[1]], colnames(tm.standard)[idx[2]],
                   tm.standard[idx[1], idx[2], j]))
    }
    bad <- tm.intervention[,,j] < -1e-12 | tm.intervention[,,j] > 1 + 1e-12
    if (any(bad)) {
      idx <- which(bad, arr.ind = TRUE)[1,]
      stop(sprintf("Intervention prob out of [0,1] at cycle %d: %s→%s = %.6f",
                   j, rownames(tm.intervention)[idx[1]], colnames(tm.intervention)[idx[2]],
                   tm.intervention[idx[1], idx[2], j]))
    }
  }
  
  #  Create a trace matrix
  trace.intervention <- matrix(data=0, nrow=cycles, ncol=n.states)
  colnames(trace.intervention) <- state.names
  
  trace.intervention[1,] <- seed%*%tm.intervention[,,1]
  
  for (i in 2:cycles) {
    trace.intervention[i,] <- trace.intervention[i-1,]%*%tm.intervention[,,i]
  }
  
  # trace.intervention
  
  outputs.intervention <- list(tm.intervention,trace.intervention)
  
  ### calculate total numbers per state
  ## want the max values (as absorbing states have multiple rows of the same value)
  ## and other states only have 1 row 
  
  # Create a new data.table with the max value of each column
  tmp.s <- as.data.table(trace.standard)
  res.standard <- data.table(t(apply(tmp.s, 2, max)))
  setnames(res.standard, names(tmp.s)) # Set the column names to match the original data.table
  
  tmp.i <- as.data.table(trace.intervention)
  res.int <- data.table(t(apply(tmp.i, 2, max)))
  setnames(res.int, names(tmp.i))
  
  rm(tmp.s)
  rm(tmp.i)
  
  ## check not losing people
  total_standard <- res.standard$no_AIS + res.standard$NOMT + res.standard$MT
  assert_that(
    isTRUE(all.equal(total_standard, n.stroke, tolerance = tolerance)),
    msg = paste0("`No intervention` strategy check failed. Check state trace matrix code")
  )
  message("`No intervention` strategy totals in state trace matrix are ok")
  
  total_int <- res.int$no_AIS + res.int$NOMT + res.int$MT
  assert_that(
    isTRUE(all.equal(total_int, n.stroke, tolerance = tolerance)),
    msg = paste0("Intervention strategy check failed. Check state trace code")
  )
  message("Intervention strategy totals in state trace matrix are ok")
  
  ### combine to one output table
  res.standard$intervention <- 0
  res.int$intervention <- 1
  res.trace <- rbind(res.standard, res.int)
  
  #### =================================================================== ####
  ####       4d. INCORPORATING & COSTING SPLITS ACROSS TYPES OF IMAGING   ####
  #### =================================================================== ####
  
  ####### Calculation of Imaging Procedure Costs and Counts #######
  
  #  Parameters from data_main for costing the tests performed
  p.ctpmri.split   <- data_main[model_param == "p.ctpmri.split", base_case] ## CTP out of CTP and MRI 
  p.ctp.plus.mri.split     <- 1 - (p.ctpmri.split) # CTP plus MRI out of CTP and MRI !!! assume no just MRI
  
  p.ctp.early   <- data_main[model_param == "p.ctp" & Presentation.Setting == "early", base_case]
  p.noCTP.early <- 1 - p.ctp.early
  
  p.ctp.late    <- data_main[model_param == "p.ctp" & Presentation.Setting == "late", base_case]
  p.noCTP.late  <- 1 - p.ctp.late
  
  
  #  Build unit cost lookup 
  unit_cost_lookup <- data_main[grepl("^c\\.", model_param), .(model_param, base_case)]
  unit_costs <- setNames(unit_cost_lookup$base_case, unit_cost_lookup$model_param)
  
  # override for updated MT cost that accounts for difference in LVO and MT
  c.mt <- data_main[model_param=="c.mt",base_case]-
    data_main[model_param=="c.lvo",base_case]
  unit_costs["c.mt"] <- c.mt
  
  
  #' Estimate the number of patients receiving MT after IVT
  #'
  #' This function calculates the total number of patients who ultimately receive
  #' mechanical thrombectomy (MT) after first receiving tissue plasminogen activator (IVT).
  #' It works by tracing the patient flow through specific states in the decision model
  #' across multiple time cycles, using the transition probabilities in a 3D array.
  #'
  #' @param trace_matrix A matrix (cycles × states) showing the number of patients
  #'                     in each state at each cycle of the model.
  #' @param tm_matrix A 3D array (states × states × cycles) of transition probabilities
  #'                  representing the likelihood of moving between states over time.
  #' @param from_state A character string naming the state in which patients receive IVT.
  #' @param to_state A character string naming the state patients enter after IVT if eligible
  #'                 for MT evaluation (e.g., an EMT assessment state).
  #'
  #' @return An integer representing the estimated total number of patients who receive MT
  #'         following IVT administration.
  #'
  #' @export
  get_MT_IVT_count <- function(trace_matrix, tm_matrix, from_state, to_state) {
    # basic checks
    stopifnot(length(dim(tm_matrix)) == 3)
    dn <- dimnames(tm_matrix)
    states <- dn[[1]]
    stopifnot(identical(states, dn[[2]]))  # square slices
    stopifnot(from_state %in% states, to_state %in% states, "MT" %in% states)
    
    # indices
    i_from <- match(from_state, states)
    i_to   <- match(to_state,   states)
    i_MT   <- match("MT",       states)
    
    n_cycles <- dim(tm_matrix)[3]
    if (n_cycles < 2) return(0)  # need at least two cycles for staged path
    
    # Flow from TPA -> EMT in cycle t:
    #   flow_from_to[t] = (# in TPA at t) * P(TPA->EMT at t)
    n_ivt        <- trace_matrix[, from_state]                 # length n_cycles
    p_from_to    <- tm_matrix[i_from, i_to, ]                  # length n_cycles
    flow_from_to <- n_ivt * p_from_to                           # expected arrivals to 'to_state' at t+1
    
    # Those patients arrive in EMT at cycle t+1; on cycle t+1 they can go EMT->MT.
    # Expected MT from these arrivals = flow_from_to[t] * P(EMT->MT at t+1)
    p_to_MT_next <- tm_matrix[i_to, i_MT, ]                    # length n_cycles
    # align (t) with (t+1)
    contrib <- flow_from_to[1:(n_cycles-1)] * p_to_MT_next[2:n_cycles]
    total <- sum(contrib, na.rm = TRUE)
    return(total)
  }
  
  
  # Run for both ASC and CSC arms, intervention and standard
  mt.ivt.asc.I <- get_MT_IVT_count(trace.intervention, tm.intervention, "IVT_EARLY_ASC","EMT_EARLY_ASC_IVT")
  mt.ivt.csc.I <- get_MT_IVT_count(trace.intervention, tm.intervention, "IVT_EARLY_CSC","EMT_EARLY_CSC_IVT")
  mt.ivt.asc.S <- get_MT_IVT_count(trace.standard, tm.standard, "IVT_EARLY_ASC","EMT_EARLY_ASC_IVT")
  mt.ivt.csc.S <- get_MT_IVT_count(trace.standard, tm.standard, "IVT_EARLY_CSC","EMT_EARLY_CSC_IVT")
  
  #  Calculate total MT IVT counts for both intervention and standard arms
  MT.IVT.I <- as.numeric(mt.ivt.asc.I + mt.ivt.csc.I)
  MT.IVT.S <- as.numeric(mt.ivt.asc.S + mt.ivt.csc.S)
  
  
  #  Define cost for each procedure 
  procedure.names <- c("IVT", "MT", "IVT + MT",
                       "NCCT + CTA",
                       "NCCT + CTA + CTP",
                       "NCCT + CTA + CTP + MRI")
  
  procedure_costs <- list(
    "c.ivt",
    "c.mt",
    c("c.ivt", "c.mt"),
    c("c.nct", "c.cta"),
    c("c.nct", "c.cta", "c.ctp"),
    c("c.nct", "c.cta", "c.ctp", "c.mri")
  )
  
  
  #  Helper to compute cost of a combo of components 
  get_unit_cost <- function(params) {
    sum(unit_costs[params], na.rm = TRUE)
  }
  
  unit.cost.vec <- mapply(get_unit_cost, procedure_costs)
  
  #  Helper function to calculate counts from res.trace 
  calculate_procedure_counts <- function(res, MT.IVT, p.ctp.early, p.ctp.late,
                                         p.noCTP.early, p.noCTP.late,
                                         p.ctpmri.split, p.ctp.plus.mri.split ) {
    return(list(
      IVT = (res$IVT_EARLY_ASC + res$IVT_EARLY_CSC) - MT.IVT, # IVT alone = IVT counts minus IVT + MT counts
      MT = res$MT - MT.IVT, # MT counts minus IVT + MT counts
      IVT_MT = MT.IVT,
      NCCT_CTA = (res$EARLY * p.noCTP.early) + (res$LATE * p.noCTP.late),
      NCCT_CTA_CTP = (res$EARLY * p.ctp.early * p.ctpmri.split) +
        (res$LATE * p.ctp.late * p.ctpmri.split),
      NCCT_CTA_CTP_MRI = (res$EARLY * p.ctp.early * p.ctp.plus.mri.split ) +
        (res$LATE * p.ctp.late * p.ctp.plus.mri.split )
    ))
  }
  
  #  Apply to trace data 
  res.int <- res.trace[intervention == 1]
  res.std <- res.trace[intervention == 0]
  
  counts.int <- calculate_procedure_counts(res.int, MT.IVT.I, p.ctp.early, p.ctp.late,
                                           p.noCTP.early, p.noCTP.late,
                                           p.ctpmri.split, p.ctp.plus.mri.split )
  
  counts.std <- calculate_procedure_counts(res.std, MT.IVT.S, p.ctp.early, p.ctp.late,
                                           p.noCTP.early, p.noCTP.late,
                                           p.ctpmri.split, p.ctp.plus.mri.split )
  
  ## check number of tests makes sense
  total_standard <- counts.std$NCCT_CTA +
    counts.std$NCCT_CTA_CTP +
    counts.std$NCCT_CTA_CTP_MRI 
  assert_that(
    isTRUE(all.equal(total_standard, res.standard$AIS, tolerance = tolerance)),
    msg = paste0("No intervention screening state check failed. Check NCCT/CTA/CTP/MRI counts")
  )
  message("`No intervention` screening state check is ok")
  
  total_intervention <- counts.int$NCCT_CTA +
    counts.int$NCCT_CTA_CTP +
    counts.int$NCCT_CTA_CTP_MRI 
  assert_that(
    isTRUE(all.equal(total_intervention, res.int$AIS, tolerance = tolerance)),
    msg = paste0("Intervention screening state check failed. Check NCCT/CTA/CTP/MRI counts")
  )
  message("`Intervention` screening state check is ok")
  
  
  # Final process results table 
  process_results <- data.table(
    procedure = procedure.names,
    intervention = unlist(counts.int),
    standard = unlist(counts.std),
    unit.cost = unit.cost.vec
  )
  
  
  #### =================================================================== ####
  ####       4e. QALY AND COST CALCULSTIONS ACROSS ARMS                     ####
  #### =================================================================== ####
  
  ##### COSTS ######
  ## cost of intervention
  n.asc <- data_main[model_param=="n.asc",base_case]
  n.csc <- data_main[model_param=="n.csc",base_case]
  
  c.B360 <- (n.asc*(data_main[model_param=="c.train",base_case]+
                      (data_main[model_param=="c.360.asc",base_case]))) + 
    (n.csc*(data_main[model_param=="c.train",base_case]+
              (data_main[model_param=="c.360.csc",base_case])))
  
  
  ### calculating LT cost and QoL impacts
  seed_distribution_ivt <- data_main[model_param=="dist.ivt"][order(mrs), base_case]
  seed_distribution_no_ivt <- data_main[model_param=="dist.noivt"][order(mrs), base_case]
  seed_distribution_MT <- data_main[model_param=="dist.mt"][order(mrs), base_case]
  seed_distribution_no_MT <- data_main[model_param=="dist.nomt"][order(mrs), base_case]
  
  out_ivt <- mrs_markov(data_main, mrs_samples_mean, seed_distribution_ivt)
  output_ivt <- out_ivt[[1]]
  mrs_ivt <- as.data.frame(out_ivt[[2]])
  
  mrs_ivt$var <- "ivt"
  mrs_ivt$t <- 1:nrow(mrs_ivt)
  
  out_no_ivt <- mrs_markov(data_main, mrs_samples_mean, seed_distribution_no_ivt)
  output_no_ivt <- out_no_ivt[[1]]
  mrs_no_ivt <- as.data.frame(out_no_ivt[[2]])
  
  mrs_no_ivt$var <- "no_ivt"                         
  mrs_no_ivt$t <- 1:nrow(mrs_no_ivt)                         
  
  out_MT <- mrs_markov(data_main, mrs_samples_mean, seed_distribution_MT)
  output_MT <- out_MT[[1]]
  mrs_MT <- as.data.frame(out_MT[[2]])
  
  mrs_MT$var <- "mt"
  mrs_MT$t <- 1:nrow(mrs_MT)
  
  out_no_MT <- mrs_markov(data_main, mrs_samples_mean, seed_distribution_no_MT)
  output_no_MT <- out_no_MT[[1]]
  mrs_no_MT <- as.data.frame(out_no_MT[[2]])
  
  mrs_no_MT$var <- "no_mt"
  mrs_no_MT$t <- 1:nrow(mrs_no_MT)
  
  ## process results
  mrs.trace <- rbind(mrs_ivt, mrs_no_ivt, mrs_MT, mrs_no_MT)
  
  diff_ivt <- output_ivt - output_no_ivt
  diff_MT <- output_MT - output_no_MT
  
  ## add to process results
  ###!!! change from normal function is in the unit costs here
  process_results$LT_cost <- c(diff_ivt$discounted.costs, diff_MT$discounted.costs,(diff_MT$discounted.costs+diff_ivt$discounted.costs),rep(0,3))
  process_results$LT_qol <- c(diff_ivt$discounted.utility, diff_MT$discounted.utility,(diff_MT$discounted.utility+diff_ivt$discounted.utility),rep(0,3))
  
  ############### RESULTS CALCULATIONS AND FORMATTING ###########
  process_results$intervention_costs <- process_results$intervention * (process_results$unit.cost + 
                                                                          process_results$LT_cost)
  process_results$standard_costs <- process_results$standard * (process_results$unit.cost + 
                                                                  process_results$LT_cost)
  process_results$intervention_qol <- process_results$intervention * process_results$LT_qol
  process_results$standard_qol <- process_results$standard * process_results$LT_qol
  
  
  ## summarise results
  summary_data_all <- process_results %>% 
    select(tail(names(.), 4)) %>% 
    summarise(
      across(where(is.numeric), sum),
      across(where(is.character), ~"Total")
    )
  
  ## add cost of intervention (1-off cost)
  summary_data_all$intervention_costs <- summary_data_all$intervention_costs + c.B360
  
  #### base case results
  incremental_results <- data.frame(inc.cost = summary_data_all$intervention_costs - summary_data_all$standard_costs,
                                    inc.qol = summary_data_all$intervention_qol - summary_data_all$standard_qol)
  incremental_results$NMB <- (incremental_results$inc.qol*data_main[model_param=="wtp",base_case]) - incremental_results$inc.cost
  
  ### return all results
  l <- list(process_results=process_results,
            summary_data_all=summary_data_all,
            incremental_results=incremental_results,
            mrs.trace=mrs.trace,
            trace.standard=trace.standard,
            trace.intervention=trace.intervention,
            tm.standard=tm.standard,
            tm.intervention=tm.intervention)
  
  return(l)
}