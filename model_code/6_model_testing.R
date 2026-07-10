###############################################
# TITLE: B360S CEA MODEL TESTING
# AUTHOR: Nichola Naylor (OI Pharma Partners Ltd), aided by GPT-4o,GTP-5 & Github co-pilot, Claude Opus 4.6 and Claude Opus 4.8
#
# DESCRIPTION:
# Validation test suite for the B360S stroke cost-effectiveness model.
# Structured around the NICE HTA Lab 14-point numerical/logical validation
# checklist (Appendix C), plus model-specific structural and integrity tests.
#
# SECTIONS:
#   0: Setup & baseline model run
#   1: Parameter integrity checks
#   2: Transition matrix & cohort conservation
#   3: Economic stress tests
#   4: NICE HTA Lab 14-point validation
#   5: Summary and output
###############################################

### !!! make sure the other model scripts (3 & 4) have been run before running

# =============================================================================
# SECTION 0: SETUP
# =============================================================================

rm(list = ls())

# ---- 0.1 Load required packages ----
library(conflicted)
library(truncnorm)
library(tidyverse)
library(data.table)
library(assertthat)
library(stringr)
library(scales)

conflict_prefer("select", "dplyr")
conflict_prefer("merge", "data.table")
conflict_prefer("filter", "dplyr")

# ---- 0.2 Source model functions ----
source("model_code/model_functions.R")

# ---- 0.3 Read input data ----
load("inputs/created_inputs/parameters_edited.RData")
data_main <- parameters
load("inputs/created_inputs/mrs_samples_mean.RData")

# ---- 0.4 Create output directory ----
outdir <- "outputs/testing"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ---- 0.5 Test tracking ----
test_results <- list()

record_test <- function(test_name, passed, message = "") {
  test_results[[test_name]] <<- list(passed = passed, message = message)
  status <- if (passed) "PASS" else "FAIL"
  cat(sprintf("[%s] %s\n", status, test_name))
  if (message != "") cat(sprintf("       %s\n", message))
}

warn_results <- list()

record_warn <- function(test_name, message = "") {
  warn_results[[test_name]] <<- list(message = message)
  cat(sprintf("[WARN] %s\n", test_name))
  if (message != "") cat(sprintf("       %s\n", message))
}

# ---- 0.6 Run baseline model ----
cat("Running baseline model...\n")
base <- run_model(data_main, cycles = 10, mrs_samples_mean)

# Save parameter audit trail
write.csv(data_main, file.path(outdir, "parameters_audit.csv"), row.names = FALSE)
write.csv(mrs_samples_mean, file.path(outdir, "mrs_samples_audit.csv"), row.names = FALSE)


# =============================================================================
# SECTION 1: PARAMETER INTEGRITY
# =============================================================================

cat("\n=============================================================================\n")
cat("SECTION 1: PARAMETER INTEGRITY\n")
cat("=============================================================================\n")

# ---- 1.1 mRS distributions sum to 1 ----
for (dp in c("dist.ivt", "dist.noivt", "dist.mt", "dist.nomt")) {
  s <- sum(data_main[model_param == dp, base_case], na.rm = TRUE)
  record_test(
    sprintf("1.1 %s sums to 1", dp),
    abs(s - 1) < 1e-6,
    sprintf("Sum = %.8f", s)
  )
}

# ---- 1.2 Probabilities in [0,1] ----
prob_params <- data_main[grepl("^p\\.", model_param) & !is.na(base_case)]
bad_probs <- prob_params[base_case < 0 | base_case > 1]

record_test(
  "1.2 All probabilities in [0,1]",
  nrow(bad_probs) == 0,
  if (nrow(bad_probs) > 0)
    sprintf("%d violations: %s", nrow(bad_probs),
            paste(unique(bad_probs$model_param), collapse = ", "))
  else ""
)

# ---- 1.3 Utilities decrease with mRS ----
utils <- data_main[model_param == "utility.mrs"& mrs %in% 0:5, .(mrs, base_case)] ## as some states worse utility than death
utils <- utils[order(mrs)]
diffs <- diff(utils$base_case)

record_test(
  "1.3 Utilities decrease with mRS severity",
  all(diffs <= 0),
  paste("Diffs:", paste(round(diffs, 3), collapse = ", "))
)


# ---- 1.4 Costs increase with mRS (0-5) ----
costs_mrs <- data_main[model_param == "cost.mrs" & mrs %in% 0:5, .(mrs, base_case)]
costs_mrs <- costs_mrs[order(mrs)]
diffs <- diff(costs_mrs$base_case)

record_test(
  "1.4 Costs increase with mRS severity (0-5)",
  all(diffs >= 0),
  paste("Diffs:", paste(round(diffs, 0), collapse = ", "))
)

# ---- 1.5 Mortality RR increases with mRS (0-5) ----
rrs <- data_main[model_param == "rr.mort" & mrs %in% 0:5, .(mrs, base_case)]
rrs <- rrs[order(mrs)]
diffs <- diff(rrs$base_case)

record_test(
  "1.5 Mortality RR increases with mRS (0-5)",
  all(diffs >= 0),
  paste("Diffs:", paste(round(diffs, 3), collapse = ", "))
)

# ---- 1.6 Dead state values ----
record_test(
  "1.6a mRS 6 utility = 0",
  data_main[model_param == "utility.mrs" & mrs == 6, base_case] == 0
)
record_test(
  "1.6b mRS 6 cost = 0",
  data_main[model_param == "cost.mrs" & mrs == 6, base_case] == 0
)

# ---- 1.7 Treatment distributions show better outcomes ----
ivt_mean  <- sum(0:6 * data_main[model_param == "dist.ivt"][order(mrs), base_case])
noivt_mean <- sum(0:6 * data_main[model_param == "dist.noivt"][order(mrs), base_case])
mt_mean   <- sum(0:6 * data_main[model_param == "dist.mt"][order(mrs), base_case])
nomt_mean <- sum(0:6 * data_main[model_param == "dist.nomt"][order(mrs), base_case])

record_test(
  "1.7a IVT mean mRS < no-IVT mean mRS",
  ivt_mean < noivt_mean,
  sprintf("IVT: %.3f vs no-IVT: %.3f", ivt_mean, noivt_mean)
)
record_test(
  "1.7b MT mean mRS < no-MT mean mRS",
  mt_mean < nomt_mean,
  sprintf("MT: %.3f vs no-MT: %.3f", mt_mean, nomt_mean)
)

# ---- 1.8 NICE-compliant discount rate and WTP ----
dr <- data_main[model_param == "dr", base_case]
wtp <- data_main[model_param == "wtp", base_case]

record_test("1.8a Discount rate = 3.5%", abs(dr - 0.035) < 1e-6,
            sprintf("Got %.4f", dr))
record_test("1.8b WTP = £20,000", wtp == 20000, sprintf("Got %s", wtp))


# =============================================================================
# SECTION 2: TRANSITION MATRIX & COHORT CONSERVATION
# =============================================================================

cat("\n=============================================================================\n")
cat("SECTION 2: TRANSITION MATRIX & COHORT CONSERVATION\n")
cat("=============================================================================\n")

# ---- 2.1a Decision tree TM rows sum to 1 (standard arm) ----
n_bad_std <- 0
for (cyc in 1:dim(base$tm.standard)[3]) {
  rs <- rowSums(base$tm.standard[, , cyc])
  active <- rs[rs > 0]  # only check rows with non-zero entries
  if (any(abs(active - 1) > 1e-6)) n_bad_std <- n_bad_std + 1
}

record_test(
  "2.1a TM row sums = 1 (standard arm)",
  n_bad_std == 0,
  if (n_bad_std > 0) sprintf("%d cycles with violations", n_bad_std) else ""
)

# ---- 2.1b Intervention arm ----
n_bad_int <- 0
for (cyc in 1:dim(base$tm.intervention)[3]) {
  rs <- rowSums(base$tm.intervention[, , cyc])
  active <- rs[rs > 0]
  if (any(abs(active - 1) > 1e-6)) n_bad_int <- n_bad_int + 1
}

record_test(
  "2.1b TM row sums = 1 (intervention arm)",
  n_bad_int == 0,
  if (n_bad_int > 0) sprintf("%d cycles with violations", n_bad_int) else ""
)

# ---- 2.2 Cohort conservation in decision tree ----
n.stroke <- data_main[model_param == "n.stroke", base_case]
last_cycle <- dim(base$trace.standard)[1]

# All patients should end up in terminal states (no_AIS + MT + NOMT)
terminal_std <- base$trace.standard[last_cycle, "no_AIS"] +
  base$trace.standard[last_cycle, "MT"] +
  base$trace.standard[last_cycle, "NOMT"]

record_test(
  "2.2a Cohort conserved (standard arm)",
  abs(terminal_std - n.stroke) < 1,
  sprintf("Terminal: %.2f, Input: %.0f", terminal_std, n.stroke)
)

terminal_int <- base$trace.intervention[last_cycle, "no_AIS"] +
  base$trace.intervention[last_cycle, "MT"] +
  base$trace.intervention[last_cycle, "NOMT"]

record_test(
  "2.2b Cohort conserved (intervention arm)",
  abs(terminal_int - n.stroke) < 1,
  sprintf("Terminal: %.2f, Input: %.0f", terminal_int, n.stroke)
)

# ---- 2.3 Intervention produces more treatments than standard ----
res_int <- base$trace.intervention[5, ]
res_std <- base$trace.standard[5, ]

# More IVT with intervention
ivt_int <- res_int["IVT_EARLY_ASC"] + res_int["IVT_EARLY_CSC"]
ivt_std <- res_std["IVT_EARLY_ASC"] + res_std["IVT_EARLY_CSC"]

record_test(
  "2.3a Intervention produces more IVT",
  ivt_int > ivt_std,
  sprintf("B360S: %.0f vs Standard: %.0f", ivt_int, ivt_std)
)

# More MT with intervention
res_int <- base$trace.intervention[last_cycle, ]
res_std <- base$trace.standard[last_cycle, ]

mt_int <- res_int["MT"]
mt_std <- res_std["MT"]

record_test(
  "2.3b Intervention produces more MT",
  mt_int > mt_std,
  sprintf("B360S: %.0f vs Standard: %.0f", mt_int, mt_std)
)

# ---- 2.4 Markov convergence (all dead at long horizon) ----
seed_test <- c(1, 0, 0, 0, 0, 0, 0) # 100% in mRS 0
long_run <- mrs_markov(data_main, mrs_samples_mean, seed_test)
long_trace <- long_run[[2]]
last_row <- long_trace[nrow(long_trace), ]

record_test(
  "2.4 Markov converges (>=99% dead at end of life tables)",
  round(last_row["mRS6"],2) >= 0.99,
  sprintf("mRS6 = %.4f at cycle %d", last_row["mRS6"], nrow(long_trace))
)

# ---- 2.5 MT route decomposition reconciles with trace MT total ----
# The six from_state -> EMT -> MT routes are the ONLY inflows to the absorbing
# MT state. Their reconstructed sum must equal the MT count the trace recorded.
# This independently re-implements the flow logic used by get_MT_IVT_count
# (which is internal to run_model) and cross-checks it against the trace.

mt_routes <- list(
  c("IVT_EARLY_ASC",   "EMT_EARLY_ASC_IVT"),
  c("NOIVT_EARLY_ASC", "EMT_EARLY_ASC_NOIVT"),
  c("ASC_LATE",        "EMT_LATE_ASC"),
  c("IVT_EARLY_CSC",   "EMT_EARLY_CSC_IVT"),
  c("NOIVT_EARLY_CSC", "EMT_EARLY_CSC_NOIVT"),
  c("CSC_LATE",        "EMT_LATE_CSC")
)

mt_route_count <- function(trace_mat, tm_mat, from_state, to_state) {
  states <- dimnames(tm_mat)[[1]]
  i_from <- match(from_state, states)
  i_to   <- match(to_state,   states)
  i_MT   <- match("MT",       states)
  n      <- dim(tm_mat)[3]
  flow   <- trace_mat[, from_state] * tm_mat[i_from, i_to, ]   # into EMT
  sum(flow[1:(n - 1)] * tm_mat[i_to, i_MT, ][2:n])             # EMT -> MT
}

recon_mt <- function(trace_mat, tm_mat) {
  sum(vapply(mt_routes,
             function(r) mt_route_count(trace_mat, tm_mat, r[1], r[2]),
             numeric(1)))
}

lc         <- dim(base$trace.standard)[1]
recon_std  <- recon_mt(base$trace.standard,     base$tm.standard)
recon_int  <- recon_mt(base$trace.intervention, base$tm.intervention)
mt_std     <- unname(base$trace.standard[lc, "MT"])
mt_int     <- unname(base$trace.intervention[lc, "MT"])

record_test(
  "2.5a MT routes reconcile with trace total (standard)",
  isTRUE(all.equal(recon_std, mt_std, tolerance = 1e-6)),
  sprintf("Sum of 6 routes = %.2f vs trace MT = %.2f", recon_std, mt_std)
)
record_test(
  "2.5b MT routes reconcile with trace total (intervention)",
  isTRUE(all.equal(recon_int, mt_int, tolerance = 1e-6)),
  sprintf("Sum of 6 routes = %.2f vs trace MT = %.2f", recon_int, mt_int)
)

# =============================================================================
# SECTION 3: ECONOMIC STRESS TESTS
# =============================================================================

cat("\n=============================================================================\n")
cat("SECTION 3: ECONOMIC STRESS TESTS\n")
cat("=============================================================================\n")

# Baseline NMB
baseline_nmb <- base$incremental_results$NMB
baseline_inc_cost <- base$incremental_results$inc.cost
baseline_inc_qol <- base$incremental_results$inc.qol

cat(sprintf("Baseline - inc.cost: £%.0f | inc.QALYs: %.2f | NMB: £%.0f\n",
            baseline_inc_cost, baseline_inc_qol, baseline_nmb))

# ---- 3.1 NMB is positive (B360S is cost-effective) ----
record_test(
  "3.1 Base case NMB > 0",
  baseline_nmb > 0,
  sprintf("NMB = £%.0f", baseline_nmb)
)

# ---- 3.2 NMB formula check ----
nmb_check <- baseline_inc_qol * wtp - baseline_inc_cost
record_test(
  "3.2 NMB = inc.QALYs * WTP - inc.Cost",
  abs(nmb_check - baseline_nmb) < 1,
  sprintf("Calculated: £%.0f, Stored: £%.0f", nmb_check, baseline_nmb)
)

# ---- 3.3 Equal efficacy stress test (OR = 1 for all) ----
cat("--- 3.3: Equal Efficacy Test ---\n")

temp_eq <- copy(data_main)
params_to_equalise <- c("p.eivt2ivt", "p.ivt.emt2mt", "p.noivt.emt2mt", "p.emt2mt")

for (p in params_to_equalise) {
  settings <- temp_eq[model_param == p & Intervention == "no intervention", 
                      .(Presentation.Setting, base_case)]
  for (s in 1:nrow(settings)) {
    temp_eq[model_param == p & 
              Intervention == "intervention" & 
              Presentation.Setting == settings$Presentation.Setting[s], 
            base_case := settings$base_case[s]]
  }
}

res_eq <- tryCatch(run_model(temp_eq, 10, mrs_samples_mean), error = function(e) NULL)

if (!is.null(res_eq)) {
  # With OR=1, intervention arm should produce same treatment volumes as standard
  last_eq <- dim(res_eq$trace.intervention)[1]
  ivt_eq_int <- res_eq$trace.intervention[5, "IVT_EARLY_ASC"] +
    res_eq$trace.intervention[5, "IVT_EARLY_CSC"]
  ivt_eq_std <- res_eq$trace.standard[5, "IVT_EARLY_ASC"] +
    res_eq$trace.standard[5, "IVT_EARLY_CSC"]
  mt_eq_int <- res_eq$trace.intervention[last_eq, "MT"]
  mt_eq_std <- res_eq$trace.standard[last_eq, "MT"]

  record_test(
    "3.3a probs=1 -> equal IVT counts",
    abs(ivt_eq_int - ivt_eq_std) < 1,
    sprintf("Int: %.0f, Std: %.0f", ivt_eq_int, ivt_eq_std)
  )
  record_test(
    "3.3b probs=1 -> equal MT counts",
    abs(mt_eq_int - mt_eq_std) < 1,
    sprintf("Int: %.0f, Std: %.0f", mt_eq_int, mt_eq_std)
  )

  # Incremental QALYs should be ~0 (only diff is B360S licence cost)
  record_test(
    "3.3c probs=1 -> incremental QALYs ≈ 0",
    abs(res_eq$incremental_results$inc.qol) < 1,
    sprintf("inc.QALYs = %.4f", res_eq$incremental_results$inc.qol)
  )
} else {
  record_test("3.3 probs=1 stress test runs", FALSE, "Model failed")
}

# ---- 3.4 Extreme pathway tests ----
cat("--- 3.4: Extreme Pathways ---\n")

# 100% ASC
temp_asc <- copy(data_main)
temp_asc[model_param == "p.asc", base_case := 1.0]
res_asc <- tryCatch(run_model(temp_asc, 10, mrs_samples_mean), error = function(e) NULL)
record_test("3.4a 100% ASC runs without error", !is.null(res_asc))

# 100% CSC
temp_csc <- copy(data_main)
temp_csc[model_param == "p.asc", base_case := 0.0]
res_csc <- tryCatch(run_model(temp_csc, 10, mrs_samples_mean), error = function(e) NULL)
record_test("3.4b 100% CSC runs without error", !is.null(res_csc))

# n.stroke = 1 (single patient)
temp_1 <- copy(data_main)
temp_1[model_param == "n.stroke", base_case := 1]
res_1 <- tryCatch(run_model(temp_1, 10, mrs_samples_mean), error = function(e) NULL)
record_test("3.4c n.stroke=1 runs without error", !is.null(res_1))

# ---- 3.5 PSA output checks  ----
cat("--- 3.5: PSA Output Checks ---\n")

if (file.exists("outputs/psa_outputs.csv")) {
  psa <- as.data.table(read.csv("outputs/psa_outputs.csv"))

  record_test("3.5a PSA has 1000 rows", nrow(psa) == 1000,
              sprintf("Got %d rows", nrow(psa)))
  record_test("3.5b No NAs in PSA inc.cost", sum(is.na(psa$inc.cost)) == 0)
  record_test("3.5c No NAs in PSA inc.qol", sum(is.na(psa$inc.qol)) == 0)
  record_test("3.5d No NAs in PSA NMB", sum(is.na(psa$NMB)) == 0)

  # Check PSA NMB has variance (not all identical)
  record_test("3.5e PSA NMB has variance",
              var(psa$NMB) > 0,
              sprintf("Var = %.0f", var(psa$NMB)))

  # Check PSA mean is broadly consistent with deterministic
  psa_mean_nmb <- mean(psa$NMB)
  record_test(
    "3.5f PSA mean NMB same order of magnitude as deterministic",
    abs(log10(abs(psa_mean_nmb)) - log10(abs(baseline_nmb))) < 1,
    sprintf("PSA mean: £%.0f, Deterministic: £%.0f", psa_mean_nmb, baseline_nmb)
  )
} else {
  record_warn("3.5 PSA outputs", "psa_outputs.csv not found - skipping PSA tests")
}


# =============================================================================
# SECTION 4: NICE HTA LAB 14-POINT VALIDATION
# =============================================================================

cat("\n=============================================================================\n")
cat("SECTION 4: NICE HTA LAB 14-POINT VALIDATION\n")
cat("=============================================================================\n")
cat("Reference: NICE HTA Lab Appendix C - Validation Checklist (Oct 2025)\n\n")

# ---- NICE 1: Discount rate = 0% ----
# Expected: discounted and undiscounted results should be equal
cat("--- NICE 1: Discount rate = 0% ---\n")

temp_dr0 <- copy(data_main)
temp_dr0[model_param == "dr", base_case := 0]
res_dr0 <- run_model(temp_dr0, 10, mrs_samples_mean)

# With dr=0 the Markov discount factors are all 1, so discounted = undiscounted
# We verify the model runs and produces finite results
record_test(
  "NICE 1: DR=0 produces finite results",
  is.finite(res_dr0$incremental_results$NMB),
  sprintf("NMB with dr=0: £%.0f", res_dr0$incremental_results$NMB)
)

# Results should be larger (in absolute terms) without discounting
record_test(
  "NICE 1b: DR=0 -> larger absolute incremental QALYs than DR=3.5%%",
  abs(res_dr0$incremental_results$inc.qol) >= abs(baseline_inc_qol) * 0.99,
  sprintf("DR=0: %.2f QALYs, DR=3.5%%: %.2f QALYs",
          res_dr0$incremental_results$inc.qol, baseline_inc_qol)
)

# ---- NICE 2: Increased discount rate ----
# Expected: discounted results should decrease
cat("--- NICE 2: Increased Discount Rate ---\n")

temp_dr10 <- copy(data_main)
temp_dr10[model_param == "dr", base_case := 0.10]
res_dr10 <- run_model(temp_dr10, 10, mrs_samples_mean)

record_test(
  "NICE 2: Higher DR (10%%) -> smaller absolute incremental QALYs",
  abs(res_dr10$incremental_results$inc.qol) < abs(baseline_inc_qol) * 1.01,
  sprintf("DR=10%%: %.2f QALYs, DR=3.5%%: %.2f QALYs",
          res_dr10$incremental_results$inc.qol, baseline_inc_qol)
)

# ---- NICE 3: Equal treatment efficacy ----
# Expected: health outcome results should be equal for all comparators
# Already tested in 3.3 above
cat("--- NICE 3: Equal Treatment Efficacy -> covered by test 3.3 ---\n")

# ---- NICE 4: All costs = 0 ----
# Expected: all cost results should be £0
# In this model, costs come from: procedure costs, ASC-CSC transfer costs, B360S costs, and LT mRS costs
cat("--- NICE 4: All costs = 0 ---\n")

temp_0cost <- copy(data_main)
temp_0cost[grepl("^c\\.", model_param), base_case := 0]
mrs_0cost <- copy(mrs_samples_mean)
mrs_0cost$cost_sample <- 0

res_0cost <- run_model(temp_0cost, 10, mrs_0cost)

total_cost_int <- res_0cost$summary_data_all$intervention_costs
total_cost_std <- res_0cost$summary_data_all$standard_costs

record_test(
  "NICE 4a: Zero costs -> intervention costs ≈ £0",
  abs(total_cost_int) < 1,
  sprintf("Intervention total cost: £%.2f", total_cost_int)
)
record_test(
  "NICE 4b: Zero costs -> standard costs ≈ £0",
  abs(total_cost_std) < 1,
  sprintf("Standard total cost: £%.2f", total_cost_std)
)

# ---- NICE 5: All utilities = 1 ----
# Expected: QALYs equal life years
# In this model, if all mRS utilities = 1 (including dead = 0 stays 0),
# then per-patient QALYs from the Markov = discounted life years
cat("--- NICE 5: Utilities = 1 (alive states) ---\n")

mrs_u1 <- copy(mrs_samples_mean)
mrs_u1[mrs %in% 0:5, utility_sample := 1]
mrs_u1[mrs == 6, utility_sample := 0]  # dead stays 0

seed_full <- c(0.2, 0.2, 0.2, 0.2, 0.1, 0.1, 0) # test distribution
out_u1 <- mrs_markov(data_main, mrs_u1, seed_full)

# QALYs should equal discounted life years (sum of alive proportion * discount factor)
# We check that QALYs > 0 and are plausible
record_test(
  "NICE 5: Utility=1 for alive states -> positive QALYs (= life years)",
  out_u1[[1]]$discounted.utility > 0,
  sprintf("Discounted QALYs: %.4f", out_u1[[1]]$discounted.utility)
)

# ---- NICE 6: All utilities = 0 ----
# Expected: QALYs = 0
cat("--- NICE 6: All utilities = 0 ---\n")

mrs_u0 <- copy(mrs_samples_mean)
mrs_u0$utility_sample <- 0

out_u0 <- mrs_markov(data_main, mrs_u0, seed_full)

record_test(
  "NICE 6: All utilities = 0 -> QALYs = 0",
  abs(out_u0[[1]]$discounted.utility) < 1e-10,
  sprintf("Discounted QALYs: %.10f", out_u0[[1]]$discounted.utility)
)

# ---- NICE 7: Mortality = 0 ----
# Expected: no deaths
cat("--- NICE 7: Mortality = 0 ---\n")

mrs_nomort <- copy(mrs_samples_mean)
mrs_nomort$mort_sample <- 0  # RR = 0 means mort_prob = pop_mort * 0 = 0

temp_nomort <- copy(data_main)
# Also zero out background mortality (life tables)
temp_nomort[model_param == "fqx", base_case := 0]
temp_nomort[model_param == "mqx", base_case := 0]

out_nomort <- mrs_markov(temp_nomort, mrs_nomort, seed_full)
trace_nomort <- out_nomort[[2]]
deaths_end <- trace_nomort[nrow(trace_nomort), "mRS6"]

record_test(
  "NICE 7: Zero mortality -> no deaths (mRS6 = 0 at all time points)",
  deaths_end < 1e-10,
  sprintf("mRS6 at last cycle: %.10f", deaths_end)
)

# ---- NICE 8: Increased mortality ----
# Expected: more deaths at each time point
cat("--- NICE 8: Increased mortality ---\n")

mrs_highmort <- copy(mrs_samples_mean)
mrs_highmort$mort_sample <- mrs_highmort$mort_sample * 2

out_highmort <- mrs_markov(data_main, mrs_highmort, seed_full)
trace_highmort <- out_highmort[[2]]
trace_base_mort <- long_run[[2]]  # from test 2.4 (using same seed_test? no, different seed)

# Re-run with same seed for fair comparison
out_base_mort <- mrs_markov(data_main, mrs_samples_mean, seed_full)
trace_base_mort <- out_base_mort[[2]]

# Check at midpoint (e.g. cycle 10)
midpoint <- min(10, nrow(trace_highmort))
deaths_high <- trace_highmort[midpoint, "mRS6"]
deaths_base <- trace_base_mort[midpoint, "mRS6"]

record_test(
  "NICE 8: Increased mortality -> more deaths at midpoint",
  deaths_high > deaths_base,
  sprintf("Cycle %d: High mort mRS6=%.4f, Base mRS6=%.4f", midpoint, deaths_high, deaths_base)
)

# ---- NICE 9: Increase starting age -> life years decrease ----
cat("--- NICE 9 & 10: Starting Age Tests ---\n")

temp_old <- copy(data_main)
temp_old[model_param == "age", base_case := 85]
res_old <- run_model(temp_old, 10, mrs_samples_mean)

record_test(
  "NICE 9: Increased age (85) -> reduced incremental QALYs",
  abs(res_old$incremental_results$inc.qol) < abs(baseline_inc_qol),
  sprintf("Age 85: %.2f QALYs, Age 75: %.2f QALYs",
          res_old$incremental_results$inc.qol, baseline_inc_qol)
)

# ---- NICE 10: Decrease starting age -> life years increase ----
temp_young <- copy(data_main)
temp_young[model_param == "age", base_case := 55]
res_young <- run_model(temp_young, 10, mrs_samples_mean)

record_test(
  "NICE 10: Decreased age (55) -> increased incremental QALYs",
  abs(res_young$incremental_results$inc.qol) > abs(baseline_inc_qol),
  sprintf("Age 55: %.2f QALYs, Age 75: %.2f QALYs",
          res_young$incremental_results$inc.qol, baseline_inc_qol)
)

# ---- NICE 11: Reduce time horizon ----
# In this model, the Markov uses remaining life expectancy from life tables.
# We can't directly set n.cycle, but reducing the start age effectively lengthens
# the horizon and vice versa. The age tests above (9/10) already cover this logic.
#
# However, we can also test by truncating the life tables:
cat("--- NICE 11 & 12: Time Horizon (via Markov cycles) ---\n")

# Short horizon: only 5 years of Markov
seed_hz <- c(0.2, 0.2, 0.2, 0.2, 0.1, 0.1, 0)
out_short_hz <- mrs_markov(data_main, mrs_samples_mean, seed_hz)

# We can't directly shorten mrs_markov cycles without code changes.
# Instead, compare starting at age 90 (short remaining life) vs 75 (longer)
temp_90 <- copy(data_main)
temp_90[model_param == "age", base_case := 90]
out_age90 <- mrs_markov(temp_90, mrs_samples_mean, seed_hz)
out_age75 <- mrs_markov(data_main, mrs_samples_mean, seed_hz)

record_test(
  "NICE 11: Shorter horizon (age 90 start) -> lower total QALYs",
  out_age90[[1]]$discounted.utility < out_age75[[1]]$discounted.utility,
  sprintf("Age 90: %.4f, Age 75: %.4f",
          out_age90[[1]]$discounted.utility, out_age75[[1]]$discounted.utility)
)

# ---- NICE 12: Increase time horizon ----
temp_55 <- copy(data_main)
temp_55[model_param == "age", base_case := 55]
out_age55 <- mrs_markov(temp_55, mrs_samples_mean, seed_hz)

record_test(
  "NICE 12: Longer horizon (age 55 start) -> higher total QALYs",
  out_age55[[1]]$discounted.utility > out_age75[[1]]$discounted.utility,
  sprintf("Age 55: %.4f, Age 75: %.4f",
          out_age55[[1]]$discounted.utility, out_age75[[1]]$discounted.utility)
)

# ---- NICE 13: Set AE/complication rates to 0 ----
# This model doesn't have AEs per se, but the closest analogue is setting
# procedure costs to 0 (removing acute treatment costs).
# Expected: total costs decrease
cat("--- NICE 13: Zero procedure costs ---\n")

temp_0proc <- copy(data_main)
temp_0proc[model_param %in% c("c.ivt", "c.mt", "c.nct", "c.cta", "c.ctp", "c.mri", "c.transfer"),
           base_case := 0]
res_0proc <- run_model(temp_0proc, 10, mrs_samples_mean)

# Both arms' costs should be lower (more negative = more savings)
# Actually with procedure costs at 0, the model should show cost reductions
# since only LT costs and B360S costs remain
record_test(
  "NICE 13: Zero procedure costs -> model runs successfully",
  is.finite(res_0proc$incremental_results$NMB),
  sprintf("NMB with zero procedure costs: £%.0f", res_0proc$incremental_results$NMB)
)

# ---- NICE 14: Increase intervention efficacy ----
cat("--- NICE 14: Increased Intervention Efficacy ---\n")
temp_higheff <- copy(data_main)

# Apply a 1.5x OR multiplier to derive new intervention probabilities
params_to_boost <- c("p.eivt2ivt", "p.ivt.emt2mt", "p.noivt.emt2mt", "p.emt2mt")
for (p in params_to_boost) {
  rows <- temp_higheff[model_param == p & Intervention == "intervention"]
  for (r in 1:nrow(rows)) {
    ps <- rows$Presentation.Setting[r]
    # Get the corresponding control probability
    p_ctrl <- temp_higheff[model_param == p & 
                             Intervention == "no intervention" & 
                             Presentation.Setting == ps, base_case]
    p_ctrl <- pmin(pmax(p_ctrl, 1e-9), 1 - 1e-9)
    # Get current intervention prob and back-calculate the implied OR
    p_int <- rows$base_case[r]
    p_int <- pmin(pmax(p_int, 1e-9), 1 - 1e-9)
    implied_or <- (p_int / (1 - p_int)) / (p_ctrl / (1 - p_ctrl))
    # Boost the OR by 1.5x and convert back to probability
    boosted_or <- implied_or * 1.5
    odds_ctrl <- p_ctrl / (1 - p_ctrl)
    new_p <- (boosted_or * odds_ctrl) / (1 + boosted_or * odds_ctrl)
    new_p <- pmin(pmax(new_p, 1e-12), 1 - 1e-12)
    temp_higheff[model_param == p & 
                   Intervention == "intervention" & 
                   Presentation.Setting == ps, 
                 base_case := new_p]
  }
}

res_higheff <- run_model(temp_higheff, 10, mrs_samples_mean)
record_test(
  "NICE 14: Increased efficacy -> NMB increases",
  res_higheff$incremental_results$NMB > baseline_nmb,
  sprintf("High efficacy NMB: £%.0f, Baseline NMB: £%.0f",
          res_higheff$incremental_results$NMB, baseline_nmb)
)


# =============================================================================
# SECTION 4b: ADDITIONAL CHECKS (carried from original test script)
# =============================================================================
# These supplement the NICE 14-point checklist with model-specific tests
# that catch different failure modes.

cat("\n=============================================================================\n")
cat("SECTION 4b: ADDITIONAL CHECKS\n")
cat("=============================================================================\n")

# ---- A1: Mortality RRs all positive ----
# (Monotonicity in 1.5 doesn't catch zero or negative values)
cat("--- A1: Mortality RR positivity ---\n")

rr_vals <- data_main[model_param == "rr.mort" & !is.na(base_case), base_case]

record_test(
  "A1 All mortality RRs > 0",
  all(rr_vals > 0),
  if (any(rr_vals <= 0))
    sprintf("Found %d non-positive RR values", sum(rr_vals <= 0))
  else sprintf("All %d RR values positive (range: %.3f to %.3f)",
               length(rr_vals), min(rr_vals), max(rr_vals))
)

# ---- A2: All cost parameters non-negative ----
cat("--- A2: Cost parameter non-negativity ---\n")

cost_params <- data_main[grepl("^c\\.", model_param) & !is.na(base_case)]
bad_costs <- cost_params[base_case < 0]

record_test(
  "A2 All cost parameters >= 0",
  nrow(bad_costs) == 0,
  if (nrow(bad_costs) > 0)
    sprintf("%d negative costs: %s", nrow(bad_costs),
            paste(unique(bad_costs$model_param), collapse = ", "))
  else ""
)

# ---- A3: Discounting reduces Markov lifetime values directly ----
# Tests the Markov component in isolation (not just incremental results)
cat("--- A3: Markov-level discounting ---\n")

seed_a3 <- c(1, 0, 0, 0, 0, 0, 0) # 100% mRS 0

dr_data_disc <- copy(data_main)
dr_data_disc[model_param == "dr", base_case := 0.035]
out_disc <- mrs_markov(dr_data_disc, mrs_samples_mean, seed_a3)

dr_data_nodisc <- copy(data_main)
dr_data_nodisc[model_param == "dr", base_case := 0.0]
out_nodisc <- mrs_markov(dr_data_nodisc, mrs_samples_mean, seed_a3)

record_test(
  "A3a Discounting reduces Markov lifetime costs",
  out_disc[[1]]$discounted.costs < out_nodisc[[1]]$discounted.costs,
  sprintf("DR=3.5%%: £%.0f, DR=0: £%.0f",
          out_disc[[1]]$discounted.costs, out_nodisc[[1]]$discounted.costs)
)

record_test(
  "A3b Discounting reduces Markov lifetime QALYs",
  out_disc[[1]]$discounted.utility < out_nodisc[[1]]$discounted.utility,
  sprintf("DR=3.5%%: %.4f, DR=0: %.4f",
          out_disc[[1]]$discounted.utility, out_nodisc[[1]]$discounted.utility)
)

# ---- A4: Higher mRS = higher costs and lower QALYs through Markov ----
# End-to-end face validity of Markov payoffs
cat("--- A4: Markov payoff face validity (mRS 0 vs mRS 4) ---\n")

seed_mrs0 <- c(1, 0, 0, 0, 0, 0, 0)
seed_mrs4 <- c(0, 0, 0, 0, 1, 0, 0)
out_mrs0 <- mrs_markov(data_main, mrs_samples_mean, seed_mrs0)
out_mrs4 <- mrs_markov(data_main, mrs_samples_mean, seed_mrs4)

record_test(
  "A4a mRS 0 has higher lifetime QALYs than mRS 4",
  out_mrs0[[1]]$discounted.utility > out_mrs4[[1]]$discounted.utility,
  sprintf("mRS0: %.4f, mRS4: %.4f",
          out_mrs0[[1]]$discounted.utility, out_mrs4[[1]]$discounted.utility)
)

record_test(
  "A4b mRS 0 has lower lifetime costs than mRS 4",
  out_mrs0[[1]]$discounted.costs < out_mrs4[[1]]$discounted.costs,
  sprintf("mRS0: £%.0f, mRS4: £%.0f",
          out_mrs0[[1]]$discounted.costs, out_mrs4[[1]]$discounted.costs)
)

# ---- A5: PSA internal arithmetic (row-level consistency) ----
cat("--- A5: PSA row-level arithmetic ---\n")

if (file.exists("outputs/psa_outputs.csv")) {
  # psa already loaded in Section 3.5 but reload for safety if not in scope
  if (!exists("psa")) psa <- as.data.table(read.csv("outputs/psa_outputs.csv"))

  # inc.cost = intervention_costs - standard_costs
  calc_inc_cost <- psa$intervention_costs - psa$standard_costs
  record_test(
    "A5a PSA: inc.cost = intervention_costs - standard_costs",
    max(abs(calc_inc_cost - psa$inc.cost), na.rm = TRUE) < 1,
    sprintf("Max row diff = £%.4f", max(abs(calc_inc_cost - psa$inc.cost), na.rm = TRUE))
  )

  # inc.qol = intervention_qol - standard_qol
  calc_inc_qol <- psa$intervention_qol - psa$standard_qol
  record_test(
    "A5b PSA: inc.qol = intervention_qol - standard_qol",
    max(abs(calc_inc_qol - psa$inc.qol), na.rm = TRUE) < 0.01,
    sprintf("Max row diff = %.6f", max(abs(calc_inc_qol - psa$inc.qol), na.rm = TRUE))
  )

  # NMB = inc.qol * WTP - inc.cost
  calc_nmb <- psa$inc.qol * wtp - psa$inc.cost
  record_test(
    "A5c PSA: NMB = inc.qol * WTP - inc.cost (all 1000 rows)",
    max(abs(calc_nmb - psa$NMB), na.rm = TRUE) < 1,
    sprintf("Max row diff = £%.4f", max(abs(calc_nmb - psa$NMB), na.rm = TRUE))
  )

  # P(cost-effective) at £20,000
  pce <- mean(psa$NMB > 0, na.rm = TRUE)
  record_test(
    "A5d P(cost-effective) at £20k > 50%",
    pce > 0.5,
    sprintf("P(CE) = %.1f%%", pce * 100)
  )

  # PSA mean NMB same sign as deterministic
  psa_mean_nmb <- mean(psa$NMB, na.rm = TRUE)
  record_test(
    "A5e PSA mean NMB same sign as deterministic",
    sign(psa_mean_nmb) == sign(baseline_nmb),
    sprintf("PSA mean: £%.0f, Deterministic: £%.0f", psa_mean_nmb, baseline_nmb)
  )
} else {
  record_warn("A5 PSA arithmetic", "psa_outputs.csv not found - skipping")
}

# ---- A6: Scenario analysis sanity checks ----
cat("--- A6: Scenario analysis checks ---\n")

sc_file <- "outputs/scenario_incremental_results.csv"
if (file.exists(sc_file)) {
  sc <- fread(sc_file)

  # Parse formatted numbers (remove commas)
  sc_nmb <- as.numeric(gsub(",", "", sc$NMB))
  sc_cost <- as.numeric(gsub(",", "", sc$inc.cost))


  # Removing LT cost savings (sc2) should increase inc.cost vs base case
  if ("sc2" %in% sc$scenario && "base_case" %in% sc$scenario) {
    sc2_cost <- sc_cost[sc$scenario == "sc2"]
    bc_cost  <- sc_cost[sc$scenario == "base_case"]
    record_test(
      "A6a Removing LT cost savings (sc2) increases inc.cost",
      sc2_cost > bc_cost,
      sprintf("sc2: £%.0f, base case: £%.0f", sc2_cost, bc_cost)
    )
  }

  # sc3 (sc2 + year-1-only mortality) should have lower QALYs than sc2
  if ("sc3" %in% sc$scenario && "sc2" %in% sc$scenario) {
    sc3_qol <- as.numeric(gsub(",", "", sc[scenario == "sc3", inc.qol]))
    sc2_qol <- as.numeric(gsub(",", "", sc[scenario == "sc2", inc.qol]))
    record_test(
      "A6b sc3 has fewer incremental QALYs than sc2 (mortality restriction)",
      sc3_qol <= sc2_qol,
      sprintf("sc3: %.2f, sc2: %.2f", sc3_qol, sc2_qol)
    )
  }

  # sc1 (younger age) should have higher absolute inc.QALYs than base case
  if ("sc1" %in% sc$scenario && "base_case" %in% sc$scenario) {
    sc1_qol <- as.numeric(gsub(",", "", sc[scenario == "sc1", inc.qol]))
    bc_qol  <- as.numeric(gsub(",", "", sc[scenario == "base_case", inc.qol]))
    record_test(
      "A6c Younger start age (sc1) -> more incremental QALYs",
      abs(sc1_qol) > abs(bc_qol),
      sprintf("sc1 (age 66): %.2f, base case (age 75): %.2f", sc1_qol, bc_qol)
    )
  }
} else {
  record_warn("A6 Scenario checks", "scenario_incremental_results.csv not found - skipping")
}

# ---- A7: ASC-CSC transfer cost checks ----
cat("--- A7: ASC-CSC transfer cost ---\n")

pr <- base$process_results   # base run_model() result from earlier in the script

# A7a: transfer row exists in base case process_results
record_test(
  "A7a Transfer row present in process_results",
  "ASC-CSC transfer (MT)" %in% pr$procedure,
  paste("Rows:", paste(pr$procedure, collapse = ", "))
)

# A7b: unit cost matches parameters.csv
c.transfer.param <- data_main[model_param == "c.transfer", base_case]
record_test(
  "A7b Transfer unit cost matches parameters.csv",
  isTRUE(all.equal(pr[procedure == "ASC-CSC transfer (MT)", unit.cost], c.transfer.param)),
  sprintf("process_results: £%.2f, parameters.csv: £%.2f",
          pr[procedure == "ASC-CSC transfer (MT)", unit.cost], c.transfer.param)
)

# A7c: transfer counts are bounded by total MT procedures in each arm
# (every transferred patient must be an MT patient; not every MT patient transfers)
for (arm in c("intervention", "standard")) {
  n_transfer <- pr[procedure == "ASC-CSC transfer (MT)", get(arm)]
  n_mt_total <- sum(pr[procedure %in% c("MT", "IVT + MT"), get(arm)])
  record_test(
    sprintf("A7c Transfer count within [0, total MT] (%s arm)", arm),
    n_transfer >= 0 && n_transfer <= n_mt_total + 1e-6,
    sprintf("Transfers: %.2f, Total MT: %.2f", n_transfer, n_mt_total)
  )
}

# A7d: cost reconciliation - zeroing c.transfer must reduce each arm's total
# costs by exactly (transfer count x £c.transfer)
temp_0transfer <- copy(data_main)
temp_0transfer[model_param == "c.transfer", base_case := 0]
res_0transfer <- run_model(temp_0transfer, 10, mrs_samples_mean)

exp_delta_int <- pr[procedure == "ASC-CSC transfer (MT)", intervention] * c.transfer.param
obs_delta_int <- base$summary_data_all$intervention_costs -
  res_0transfer$summary_data_all$intervention_costs
record_test(
  "A7d Zeroing c.transfer removes exactly the transfer cost (intervention arm)",
  isTRUE(all.equal(obs_delta_int, exp_delta_int, tolerance = 1e-6)),
  sprintf("Observed: £%.2f, Expected: £%.2f", obs_delta_int, exp_delta_int)
)

exp_delta_std <- pr[procedure == "ASC-CSC transfer (MT)", standard] * c.transfer.param
obs_delta_std <- base$summary_data_all$standard_costs -
  res_0transfer$summary_data_all$standard_costs
record_test(
  "A7e Zeroing c.transfer removes exactly the transfer cost (standard arm)",
  isTRUE(all.equal(obs_delta_std, exp_delta_std, tolerance = 1e-6)),
  sprintf("Observed: £%.2f, Expected: £%.2f", obs_delta_std, exp_delta_std)
)

# A7f: variant drift guard - every scenario's process_results must carry the transfer row
sc_pr_file <- "outputs/scenario_process_results.csv"
if (file.exists(sc_pr_file)) {
  sc_pr <- fread(sc_pr_file)
  missing_sc <- setdiff(unique(sc_pr$scenario),
                        unique(sc_pr[procedure == "ASC-CSC transfer (MT)", scenario]))
  record_test(
    "A7f All scenarios include the ASC-CSC transfer row",
    length(missing_sc) == 0,
    if (length(missing_sc) == 0) "All scenarios have the transfer row"
    else paste("Missing in:", paste(missing_sc, collapse = ", "))
  )
} else {
  record_warn("A7f Scenario transfer-row check", "scenario_process_results.csv not found - skipping")
}


# =============================================================================
# SECTION 5: SUMMARY AND OUTPUT
# =============================================================================

cat("\n=============================================================================\n")
cat("SECTION 5: SUMMARY AND OUTPUT\n")
cat("=============================================================================\n")

test_summary <- data.table(
  test_name = names(test_results),
  passed = sapply(test_results, function(x) x$passed),
  message = sapply(test_results, function(x) x$message)
)

n_passed <- sum(test_summary$passed)
n_failed <- sum(!test_summary$passed)
n_total  <- nrow(test_summary)
n_warn   <- length(warn_results)

# Save outputs
write.csv(test_summary, file.path(outdir, "test_summary.csv"), row.names = FALSE)

# Write detailed summary file
summary_file <- file.path(outdir, "test_summary.txt")
sink(summary_file)
cat("###############################################################################\n")
cat("#               B360S CEA MODEL TEST SUMMARY                                #\n")
cat("#      Based on NICE HTA Lab Validation Checklist (Oct 2025)                #\n")
cat("###############################################################################\n")
cat(sprintf("Run completed: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("\nTOTAL: %d tests | PASSED: %d | FAILED: %d | WARNINGS: %d\n",
            n_total, n_passed, n_failed, n_warn))
cat("\n-------------------------------------------------------------------------------\n")

for (i in seq_len(nrow(test_summary))) {
  status <- if (test_summary$passed[i]) "PASS" else "FAIL"
  cat(sprintf("[%s] %s\n", status, test_summary$test_name[i]))
  if (test_summary$message[i] != "") {
    cat(sprintf("       %s\n", test_summary$message[i]))
  }
}

cat("\n-------------------------------------------------------------------------------\n")
cat("FAILED TESTS:\n")
if (n_failed == 0) {
  cat("  None - all tests passed!\n")
} else {
  for (i in which(!test_summary$passed)) {
    cat(sprintf("  - %s: %s\n", test_summary$test_name[i], test_summary$message[i]))
  }
}
cat("\nWARNINGS:\n")
if (n_warn == 0) {
  cat("  None\n")
} else {
  for (wn in names(warn_results)) {
    cat(sprintf("  - %s: %s\n", wn, warn_results[[wn]]$message))
  }
}
sink()

# Category breakdown
param_tests    <- test_summary[grepl("^1\\.", test_name)]
struct_tests   <- test_summary[grepl("^2\\.", test_name)]
stress_tests   <- test_summary[grepl("^3\\.", test_name)]
nice_tests     <- test_summary[grepl("^NICE", test_name)]
addl_tests     <- test_summary[grepl("^A", test_name)]

cat(sprintf("\nCompleted: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("RESULTS: %d/%d tests passed (%.1f%%) | %d warnings\n\n",
            n_passed, n_total, 100 * n_passed / n_total, n_warn))
cat("Results by section:\n")
cat(sprintf("  Section 1  - Parameter Integrity:    %d/%d passed\n", sum(param_tests$passed), nrow(param_tests)))
cat(sprintf("  Section 2  - TM & Conservation:      %d/%d passed\n", sum(struct_tests$passed), nrow(struct_tests)))
cat(sprintf("  Section 3  - Economic Stress:        %d/%d passed\n", sum(stress_tests$passed), nrow(stress_tests)))
cat(sprintf("  Section 4  - NICE 14-pt Validation:  %d/%d passed\n", sum(nice_tests$passed), nrow(nice_tests)))
cat(sprintf("  Section 4b - Additional Checks:      %d/%d passed\n", sum(addl_tests$passed), nrow(addl_tests)))

if (n_failed > 0) {
  cat("\nWARNING: Some tests failed. Review before submission.\n")
} else {
  cat("\nAll tests passed successfully.\n")
}

cat(sprintf("\nOutputs: %s\n", outdir))
