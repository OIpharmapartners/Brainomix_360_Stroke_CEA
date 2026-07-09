
###############################################
# TITLE: Clean table creation
# AUTHOR: Nichola Naylor (OI Pharma Partners Ltd), aided by GPT-4o,GTP-5 & Github co-pilot, Claude Opus 4.6 and Claude Opus 4.8

source("model_code/4_scenario_analysis.R")

#### ======================================= ####
#### TABLE 3 for manuscript: National Results  ####
#### ======================================= ####

proc <- read.csv("outputs/proc_results.csv")
proc <- as.data.table(proc)

# Keep only treatment procedures (not imaging)
t3 <- proc[procedure %in% c("IVT", "MT", "IVT + MT")]

# Get n.ais from imaging rows
n.ais.t3 <- sum(proc[procedure %in% c("NCCT + CTA", "NCCT + CTA + CTP", 
                                      "NCCT + CTA + CTP + MRI"), intervention])

# Build Table 3
# Build Table 3
# Procedure counts and percentages are reported together, e.g. 9,109 (11.2%)
table3 <- data.table(
  Procedures = c("IVT only", "MT only", "IVT plus MT"),
  
  `B360S Procedures (% AIS)` = paste0(
    comma(round(t3$intervention, 0)),
    " (",
    round(t3$intervention / n.ais.t3 * 100, 1),
    "%)"
  ),
  
  `Standard Care Procedures (% AIS)` = paste0(
    comma(round(t3$standard, 0)),
    " (",
    round(t3$standard / n.ais.t3 * 100, 1),
    "%)"
  ),
  
  `B360S Procedure Costs (£)` = comma(round(t3$intervention * t3$unit.cost, 0)),
  `Standard Care Procedure Costs (£)` = comma(round(t3$standard * t3$unit.cost, 0)),
  `B360S Long-term Costs (£)` = comma(round(t3$intervention * t3$LT_cost, 0)),
  `Standard Care Long-term Costs (£)` = comma(round(t3$standard * t3$LT_cost, 0)),
  `B360S QALYs` = comma(round(t3$intervention_qol, 0)),
  `Standard Care QALYs` = comma(round(t3$standard_qol, 0))
)


write.csv(table3, "outputs/table3_clean.csv", row.names = FALSE)
print(table3)



#### ======================================= ####
#### TABLE 4: Scenario Analyses + PSA        ####
#### ======================================= ####

sc_raw <- rbind(
  base_case$incremental_results,
  sc1$incremental_results,
  sc2$incremental_results,
  sc3$incremental_results,
  sc4$incremental_results,
  sc5$incremental_results,
  sc6$incremental_results,
  sc7$incremental_results,
  sc8$incremental_results,
  sc9_inc
)
sc_raw <- as.data.table(sc_raw)

# Scenario labels matching V9 manuscript style
sc_labels <- c(
  "base_case" = "Base case",
  "sc1" = "SC1: Start age 66",
  "sc2" = "SC2: Long-term cost savings from IVT and MT removed",
  "sc3" = "SC3: SC2 plus mortality savings from IVT and MT only occur in year 1",
  "sc4" = "SC4: Additive MT and IVT benefits for those who have both",
  "sc5" = "SC5: Alternative IVT mRS distribution",
  "sc6" = "SC6: Gradual mortality risk attenuation (years 1-4)",
  "sc7" = "SC7: Half-cycle correction applied",
  "sc8" = "SC8: Dampened intervention effect",
  "sc9" = "SC9: Per-AIS patient"
)

# Format ICER:
# - Dominant = cost-saving and QALY-gaining
# - Dominated = more costly and less effective
# - Otherwise calculate inc.cost / inc.qol
fmt_icer <- function(cost, qaly, dp = 0) {
  out <- rep(NA_character_, length(cost))
  
  out[cost < 0 & qaly > 0] <- "Dominant"
  out[cost > 0 & qaly < 0] <- "Dominated"
  out[qaly == 0] <- "Not calculated"
  
  calc <- is.na(out) & !is.na(cost) & !is.na(qaly)
  out[calc] <- paste0("£", comma(round(cost[calc] / qaly[calc], dp)))
  
  return(out)
}

# Format deterministic scenario rows
sc_pop <- sc_raw[scenario != "sc9"]  # population-level scenarios
sc_pp  <- sc_raw[scenario == "sc9"]  # per-patient - needs different edit

table4 <- data.table(
  Scenario                    = sc_labels[sc_pop$scenario],
  `Incremental Cost (£)`      = comma(round(sc_pop$inc.cost, 0)),
  `Incremental QALYs`         = comma(round(sc_pop$inc.qol, 0)),
  `Net Monetary Benefit (£)`  = comma(round(sc_pop$NMB, 0)),
  `ICER (£/QALY)`             = fmt_icer(sc_pop$inc.cost, sc_pop$inc.qol)
)

# Per-patient row with appropriate precision
sc9_row <- data.table(
  Scenario                    = sc_labels["sc9"],
  `Incremental Cost (£)`      = comma(round(sc_pp$inc.cost, 0)),
  `Incremental QALYs`         = sprintf("%.3f", sc_pp$inc.qol),
  `Net Monetary Benefit (£)`  = comma(round(sc_pp$NMB, 0)),
  `ICER (£/QALY)`             = fmt_icer(sc_pp$inc.cost, sc_pp$inc.qol)
)

table4 <- rbind(table4, sc9_row)

# Add PSA row with mean (95% CI)
# psa_summary created in scenario analysis
fmt_psa <- function(met, dp = 0) {
  m <- psa_summary$mean[psa_summary$metric == met]
  l <- psa_summary$lower[psa_summary$metric == met]
  u <- psa_summary$upper[psa_summary$metric == met]
  paste0(comma(round(m, dp)), " (", comma(round(l, dp)), " to ", comma(round(u, dp)), ")")
}

psa_inc_cost <- psa_summary$mean[psa_summary$metric == "inc.cost"]
psa_inc_qaly <- psa_summary$mean[psa_summary$metric == "inc.qol"]

psa_row <- data.table(
  Scenario                    = "PSA: Mean (95% CI)",
  `Incremental Cost (£)`      = fmt_psa("inc.cost"),
  `Incremental QALYs`         = fmt_psa("inc.qol"),
  `Net Monetary Benefit (£)`  = fmt_psa("NMB"),
  `ICER (£/QALY)`             = fmt_icer(psa_inc_cost, psa_inc_qaly)
)

table4 <- rbind(table4, psa_row)

write.csv(table4, "outputs/table4_clean.csv", row.names = FALSE)
print(table4)


#### information for narrative text
# % increase in procedures: B360S vs Standard Care
t3_pct <- data.table(
  Procedure = t3$procedure,
  pct_increase = round((t3$intervention - t3$standard) / t3$standard * 100, 1)
)
print(t3_pct)

# Total IVT = IVT only + IVT+MT
total_ivt_int <- t3[procedure == "IVT", intervention] + t3[procedure == "IVT + MT", intervention]
total_ivt_std <- t3[procedure == "IVT", standard] + t3[procedure == "IVT + MT", standard]
pct_ivt_total <- round((total_ivt_int - total_ivt_std) / total_ivt_std * 100, 1)

# Total MT = MT only + IVT+MT
total_mt_int <- t3[procedure == "MT", intervention] + t3[procedure == "IVT + MT", intervention]
total_mt_std <- t3[procedure == "MT", standard] + t3[procedure == "IVT + MT", standard]
pct_mt_total <- round((total_mt_int - total_mt_std) / total_mt_std * 100, 1)

cat(sprintf("Total IVT increase: %.1f%%\n", pct_ivt_total))
cat(sprintf("Total MT increase: %.1f%%\n", pct_mt_total))

icer_sc2 <- sc2$incremental_results$inc.cost / sc2$incremental_results$inc.qol
icer_sc3 <- sc3$incremental_results$inc.cost / sc3$incremental_results$inc.qol
sprintf("SC2 ICER: £%.0f, SC3 ICER: £%.0f", icer_sc2, icer_sc3)

#### ================================================================= ####
#### SUPPLEMENTARY TABLE: Full cost breakdown + incremental             ####
####   All procedures + scans + ASC-CSC transfer. Immediate (count x    ####
####   unit.cost), long-term (count x LT_cost), and total costs, plus   ####
####   the incremental (B360S - Standard) so the table reconciles to    ####
####   both the population totals AND the ICER / NMB.                    ####
#### ================================================================= ####


supp_num <- proc[, .(
  Procedure                      = procedure,
  `B360S N`                      = intervention,
  `Standard Care N`              = standard,
  `B360S Immediate (£)`          = intervention * unit.cost,
  `Standard Care Immediate (£)`  = standard     * unit.cost,
  `B360S Long-term (£)`          = intervention * LT_cost,
  `Standard Care Long-term (£)`  = standard     * LT_cost,
  `B360S Total (£)`              = intervention * (unit.cost + LT_cost),
  `Standard Care Total (£)`      = standard     * (unit.cost + LT_cost),
  `Incremental Immediate (£)`    = (intervention - standard) * unit.cost,
  `Incremental Long-term (£)`    = (intervention - standard) * LT_cost,
  `Incremental Total (£)`        = intervention * (unit.cost + LT_cost) -
    standard     * (unit.cost + LT_cost)
)]

# one-off B360S software + training cost (intervention arm only)

  n.asc <- data_main[model_param == "n.asc", base_case]
  n.csc <- data_main[model_param == "n.csc", base_case]
  c.B360 <- n.asc * (data_main[model_param == "c.train",   base_case] +
             data_main[model_param == "c.360.asc", base_case]) +
    n.csc * (data_main[model_param == "c.train",   base_case] +
               data_main[model_param == "c.360.csc", base_case])

cost_cols <- c("B360S Immediate (£)", "Standard Care Immediate (£)",
               "B360S Long-term (£)", "Standard Care Long-term (£)",
               "B360S Total (£)",     "Standard Care Total (£)",
               "Incremental Immediate (£)", "Incremental Long-term (£)",
               "Incremental Total (£)")

subtotal <- c(list(Procedure = "Subtotal (procedures)",
                   `B360S N` = NA_real_, `Standard Care N` = NA_real_),
              as.list(colSums(supp_num[, ..cost_cols])))

extra <- list(subtotal)

  oneoff <- list(Procedure = "B360S software & training (one-off)",
                 `B360S N` = NA_real_, `Standard Care N` = NA_real_,
                 `B360S Immediate (£)` = c.B360, `Standard Care Immediate (£)` = 0,
                 `B360S Long-term (£)` = 0,      `Standard Care Long-term (£)` = 0,
                 `B360S Total (£)` = c.B360,     `Standard Care Total (£)` = 0,
                 `Incremental Immediate (£)` = c.B360,
                 `Incremental Long-term (£)` = 0,
                 `Incremental Total (£)` = c.B360)
  grand <- list(Procedure = "Grand total",
                `B360S N` = NA_real_, `Standard Care N` = NA_real_,
                `B360S Immediate (£)` = subtotal$`B360S Immediate (£)` + c.B360,
                `Standard Care Immediate (£)` = subtotal$`Standard Care Immediate (£)`,
                `B360S Long-term (£)` = subtotal$`B360S Long-term (£)`,
                `Standard Care Long-term (£)` = subtotal$`Standard Care Long-term (£)`,
                `B360S Total (£)` = subtotal$`B360S Total (£)` + c.B360,
                `Standard Care Total (£)` = subtotal$`Standard Care Total (£)`,
                `Incremental Immediate (£)` = subtotal$`Incremental Immediate (£)` + c.B360,
                `Incremental Long-term (£)` = subtotal$`Incremental Long-term (£)`,
                `Incremental Total (£)` = subtotal$`Incremental Total (£)` + c.B360)
  extra <- list(subtotal, oneoff, grand)

supp_num <- rbind(supp_num, rbindlist(extra, use.names = TRUE), use.names = TRUE)

# format: thousands separators, signed, blank where not applicable
fmt <- function(x) ifelse(is.na(x), "", comma(round(x, 0)))
tableS <- copy(supp_num)
for (cc in setdiff(names(tableS), "Procedure")) tableS[[cc]] <- fmt(tableS[[cc]])

write.csv(tableS, "outputs/tableS_cost_breakdown.csv", row.names = FALSE)
print(tableS)

# --- CEA reconciliation: the table's grand-total incremental IS inc.cost ---
  ir       <- base_case$incremental_results
  inc_tab  <- supp_num[Procedure == "Grand total", `Incremental Total (£)`]
  wtp      <- if (exists("data_main")) data_main[model_param == "wtp", base_case] else 20000
  stopifnot(isTRUE(all.equal(inc_tab, ir$inc.cost, tolerance = 1)))
  cat(sprintf(paste0(
    "\n--- CEA reconciliation (traceable from Table S) ---\n",
    "Incremental cost (table grand total): \u00a3%s\n",
    "Incremental cost (model):             \u00a3%s\n",
    "Incremental QALYs:                     %s\n",
    "ICER (inc.cost / inc.QALYs):          \u00a3%s per QALY%s\n",
    "NMB at WTP \u00a3%s:                     \u00a3%s\n"),
    comma(round(inc_tab, 0)), comma(round(ir$inc.cost, 0)),
    comma(round(ir$inc.qol, 1)),
    comma(round(ir$inc.cost / ir$inc.qol, 0)),
    ifelse(ir$inc.cost < 0 & ir$inc.qol > 0, " (dominant)", ""),
    comma(wtp), comma(round(ir$NMB, 0))))
