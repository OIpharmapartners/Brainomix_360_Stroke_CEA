###############################################
# TITLE: Journal-ready DSA tornado plot
# DESCRIPTION:
# - Keeps the same logic as the original fig_DSA.R
# - Uses all DSA parameters available in outputs/DSA_results.RData
# - Keeps manual relabelling
# - Removes duplicate parameter labels using the largest DSA range, as before
# - Fixes encoding issue by using "GBP" instead of the pound symbol
# - Saves PNG, TIFF, PDF, and plot-data CSV
###############################################

rm(list = ls())

#### Load required packages ####

library(data.table)
library(ggplot2)
library(scales)
library(stringr)

#### Load model outputs ####

load("outputs/DSA_results.RData")          # loads object: DSA
load("outputs/base_case_results.RData")    # loads object: base_case

if (!dir.exists("outputs")) dir.create("outputs", recursive = TRUE)

#### Prepare DSA data ####

DSA <- as.data.table(DSA)

# Required columns check
required_cols <- c(
  "variable", "Parameter",
  "low.value", "high.value",
  "midpoint", "ymin", "ymax"
)

missing_cols <- setdiff(required_cols, names(DSA))
if (length(missing_cols) > 0) {
  stop(
    "DSA is missing required column(s): ",
    paste(missing_cols, collapse = ", ")
  )
}

#### Manual relabelling retained from original script ####

DSA[variable == "c.360.asc",
    Parameter := "Cost per B360S for ASC [15000 , 30000]"]

DSA[variable == "c.360.csc",
    Parameter := "Cost per B360S for CSC [30000 , 60000]"]

DSA[variable == "p.eivt",
    Parameter := "% eligible for IVT [0.21 , 0.88]"]

DSA[variable == "p.ivt2emt",
    Parameter := "% of IVT patients that are eligible for MT [0.41 , 0.8]"]

DSA[variable == "p.noivt2emt",
    Parameter := "% of patients who didn't get IVT who are eligible for MT [0.41 , 0.8]"]

DSA[variable == "p.emt",
    Parameter := "% of late stroke patients in CSCs that are eligible for MT [0.41 , 0.8]"]

#### Clean text encoding ####
# This avoids the "invalid input ... utf8towcs" error when saving PDF/TIFF.
# It also prevents hidden bad characters in parameter labels from breaking export.

DSA[, Parameter := enc2utf8(as.character(Parameter))]
DSA[, Parameter := gsub("\uFFFD", "GBP", Parameter)]          # replacement character
DSA[, Parameter := gsub("\u00A3", "GBP", Parameter, fixed = TRUE)]  # pound sign
DSA[, Parameter := gsub("£", "GBP", Parameter, fixed = TRUE)]       # visible pound sign

#### Remove duplicate rows, as in the original script ####
# This is NOT a top-10 filter.
# It keeps one row per displayed Parameter label, choosing the widest tornado range.

DSA[, impact := abs(high.value - low.value)]
DSA <- DSA[, .SD[which.max(impact)], by = Parameter]

#### Order and wrap labels ####

DSA <- DSA[order(impact)]   # smallest first so largest appears at top after coord_flip()

DSA[, Parameter_plot := stringr::str_wrap(Parameter, width = 42)]

# Use a plotting ID rather than making labels unique with ".1", ".2", etc.
# This keeps wrapped labels clean.
DSA[, plot_id := factor(seq_len(.N), levels = seq_len(.N))]
axis_labels <- setNames(DSA$Parameter_plot, levels(DSA$plot_id))

#### Base-case NMB line ####

base_case_value <- base_case$incremental_results$NMB

if (length(base_case_value) != 1 || is.na(base_case_value)) {
  stop("base_case$incremental_results$NMB must be a single non-missing value.")
}

#### Dynamic figure height ####
# Makes the figure taller when all parameters are shown, avoiding label overlap.

height_mm <- max(180, 7.5 * nrow(DSA) + 45)

#### Create tornado plot ####

options(scipen = 999)

dsa_plot <- ggplot(
  DSA,
  aes(x = plot_id, y = midpoint)
) +
  geom_linerange(
    aes(ymin = ymin, ymax = ymax),
    colour = "lightskyblue",
    linewidth = 3.5,
    lineend = "round"
  ) +
  geom_point(
    aes(y = low.value),
    colour = "#E69F00",
    size = 2.2
  ) +
  geom_point(
    aes(y = high.value),
    colour = "#009E73",
    size = 2.2
  ) +
  geom_hline(
    yintercept = base_case_value,
    colour = "black",
    linewidth = 0.4
  ) +
  coord_flip(clip = "off") +
  scale_x_discrete(
    labels = axis_labels
  ) +
  scale_y_continuous(
    name = "Net monetary benefit (GBP million)",
    labels = scales::label_number(
      scale = 1e-6,
      suffix = "m",
      accuracy = 1
    ),
    breaks = scales::pretty_breaks(n = 5),
    expand = expansion(mult = c(0.03, 0.08))
  ) +
  labs(
    x = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    
    axis.text.y = element_text(
      size = 7,
      colour = "black",
      lineheight = 0.9
    ),
    axis.text.x = element_text(
      size = 8,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 9,
      face = "bold",
      margin = margin(t = 8)
    ),
    
    axis.line = element_line(
      colour = "black",
      linewidth = 0.3
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.3
    ),
    
    panel.grid.major.x = element_line(
      colour = "grey85",
      linewidth = 0.3,
      linetype = "dotted"
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    
    plot.margin = margin(t = 5, r = 12, b = 5, l = 12)
  )

print(dsa_plot)

#### Save outputs ####

# Keep original output name for compatibility with report/user guide
ggsave(
  filename = "outputs/dsa_plot.png",
  plot = dsa_plot,
  width = 190,
  height = height_mm,
  units = "mm",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

# Journal-ready TIFF
ggsave(
  filename = "outputs/dsa_plot_journal.tiff",
  plot = dsa_plot,
  width = 190,
  height = height_mm,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white",
  limitsize = FALSE
)

# Vector PDF
ggsave(
  filename = "outputs/dsa_plot_journal.pdf",
  plot = dsa_plot,
  width = 190,
  height = height_mm,
  units = "mm",
  device = grDevices::cairo_pdf,
  bg = "white",
  limitsize = FALSE
)

# Optional journal PNG copy
ggsave(
  filename = "outputs/dsa_plot_journal.png",
  plot = dsa_plot,
  width = 190,
  height = height_mm,
  units = "mm",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

#### Export tidy CSV of the values used in the plot ####

dsa_export <- DSA[, .(
  Parameter,
  low.value,
  high.value,
  midpoint
)]

# Order parameters by absolute swing, largest first, for readability
dsa_export <- dsa_export[
  order(abs(high.value - low.value), decreasing = TRUE)
]

write.csv(
  dsa_export,
  file = "outputs/dsa_plot_data.csv",
  row.names = FALSE
)