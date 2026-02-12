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
# Note: costs are negative (savings) in the model output, so we take abs()
# and the header says "Costs Saved"
table3 <- data.table(
  Procedures           = c("IVT only", "MT only", "IVT plus MT"),
  `B360S Procedures`   = comma(round(t3$intervention, 0)),
  `B360S (% AIS)`      = paste0(round(t3$intervention / n.ais.t3 * 100, 1), "%"),
  `Standard Care Procedures` = comma(round(t3$standard, 0)),
  `Standard Care (% AIS)` = paste0(round(t3$standard / n.ais.t3 * 100, 1), "%"),
  `B360S Costs Saved (£)` = comma(round(abs(t3$intervention_costs), 0)),
  `Standard Care Costs Saved (£)` = comma(round(abs(t3$standard_costs), 0)),
  `B360S QALYs`        = comma(round(t3$intervention_qol, 0)),
  `Standard Care QALYs` = comma(round(t3$standard_qol, 0))
)

write.csv(table3, "outputs/table3_clean.csv", row.names = FALSE)
print(table3)


#### ============================================================ ####
#### CHUNK 6: Clean Table 4 — Scenario Results + PSA              ####
####          Run after CHUNK 3 and CHUNK 4                       ####
#### ============================================================ ####

#### ======================================= ####
#### TABLE 4: Scenario Analyses + PSA        ####
#### ======================================= ####

# Read back the scenario results (already formatted as strings from Chunk 3)
# OR if running in same session, use the raw incremental.results before formatting
# Here we rebuild from the raw model outputs for cleaner control:

sc_raw <- rbind(
  base_case$incremental_results,
  sc1$incremental_results,
  sc2$incremental_results,
  sc3$incremental_results,
  sc4$incremental_results,
  sc5$incremental_results,
  sc6$incremental_results,
  sc7$incremental_results,
  sc8_inc
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
  "sc8" = "SC8: Per-AIS patient"
  
)

# Format deterministic scenario rows
# SC8 (per-patient) needs different precision: costs ~£111, QALYs ~0.019
# So handle SC8 separately

sc_pop <- sc_raw[scenario != "sc8"]  # population-level scenarios
sc_pp  <- sc_raw[scenario == "sc8"]  # per-patient

table4 <- data.table(
  Scenario               = sc_labels[sc_pop$scenario],
  `Incremental Cost (£)` = comma(round(sc_pop$inc.cost, 0)),
  `Incremental QALYs`    = comma(round(sc_pop$inc.qol, 0)),
  `Net Monetary Benefit (£)` = comma(round(sc_pop$NMB, 0))
)

# Per-patient row with appropriate precision
sc8_row <- data.table(
  Scenario               = sc_labels["sc8"],
  `Incremental Cost (£)` = comma(round(sc_pp$inc.cost, 0)),
  `Incremental QALYs`    = sprintf("%.3f", sc_pp$inc.qol),
  `Net Monetary Benefit (£)` = comma(round(sc_pp$NMB, 0))
)
table4 <- rbind(table4, sc8_row)

# Add PSA row with mean (95% CI)
# psa_summary created in Chunk 4
fmt_psa <- function(met, dp = 0) {
  m <- psa_summary$mean[psa_summary$metric == met]
  l <- psa_summary$lower[psa_summary$metric == met]
  u <- psa_summary$upper[psa_summary$metric == met]
  paste0(comma(round(m, dp)), " (", comma(round(l, dp)), " to ", comma(round(u, dp)), ")")
}

psa_row <- data.table(
  Scenario               = "PSA: Mean (95% CI)",
  `Incremental Cost (£)` = fmt_psa("inc.cost"),
  `Incremental QALYs`    = fmt_psa("inc.qol"),
  `Net Monetary Benefit (£)` = fmt_psa("NMB")
)

table4 <- rbind(table4, psa_row)

write.csv(table4, "outputs/table4_clean.csv", row.names = FALSE)
print(table4)

