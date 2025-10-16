###############################################
# TITLE: Tornado plotting
## This script runs the DSA tornado plot code
## Inputs: outputs/DSA_results.RData
## Outputs: outputs/dsa_plot.png

### Clear environment
rm(list=ls())

## Load required packages
library(data.table)
library(ggplot2)

load("outputs/DSA_results.RData")
load("outputs/base_case_results.RData")

DSA <- as.data.table(DSA)
## relabel and 
DSA[variable=="c.360.asc", Parameter := "Cost per B360S for ASC [15000 , 30000]"]
DSA[variable=="c.360.csc", Parameter := "Cost per B360S for CSC  [30000 , 60000]"]
DSA[variable=="p.eivt", Parameter := "% eligible for IVT [0.21 , 0.88]"]
DSA[variable=="p.ivt2emt", Parameter := "% of IVT patients that are eligible for MT [0.41 , 0.8]"]
DSA[variable=="p.noivt2emt", Parameter := "% of patients who didn't get IVT who are eligble for MT [0.41 , 0.8]"]


## remove duplicate rows
DSA <- DSA[, .SD[which.max(abs(high.value - low.value))], by = Parameter]

# Set the base case value
base_case_value <- base_case$incremental_results$NMB

# Create the tornado plot
options(scipen = 999)
ggplot(DSA, aes(x = reorder(Parameter, abs(high.value - low.value)), y = midpoint)) +
  geom_linerange(aes(ymin = ymin, ymax = ymax), color = "lightskyblue", size = 5) +
  geom_point(aes(y = low.value), color = "red", size = 3) +
  geom_point(aes(y = high.value), color = "green", size = 3) +
  geom_hline(yintercept = base_case_value, color = "black")+
  coord_flip() +scale_y_continuous(labels = scales::comma)+
  labs(x = "Variable", y = "Net Monetary Benefit (NMB)",
       title = "Tornado Plot") 



dsa_plot <- ggplot(DSA, aes(x = reorder(Parameter, abs(high.value - low.value)), y = midpoint)) +
  geom_linerange(aes(ymin = ymin, ymax = ymax), color = "lightskyblue", size = 4) +       # blue
  geom_point(aes(y = low.value), color = "#E69F00", size = 2.8) +                     # orange
  geom_point(aes(y = high.value), color = "#009E73", size = 2.8) +                    # green
  geom_hline(yintercept = base_case_value, color = "black") +
  coord_flip() +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x = "Variable",
    y = "Net Monetary Benefit (£)"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.y = element_text(size = 20, face = "bold", color = "black"),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    panel.grid.major.x = element_line(color = "gray80", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 15)
  )

dsa_plot
ggsave("outputs/dsa_plot.png", 
       plot = dsa_plot, 
       width = 15, 
       height = 14, 
       units = "in", 
       dpi = 300, 
       limitsize = FALSE)


## export a tidy CSV of the values used in the plot (for supplement)
dsa_export <- DSA[, c("Parameter", "low.value", "high.value", "midpoint")]
# Order parameters by absolute swing (largest first) for readability
dsa_export <- dsa_export[order(abs(dsa_export$high.value - dsa_export$low.value), decreasing = TRUE), ]
write.csv(dsa_export, file = "outputs/dsa_plot_data.csv", row.names = FALSE)