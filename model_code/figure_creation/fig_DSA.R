###############################################
# TITLE: Tornado plotting
## This script runs the DSA tornado plot code
## Inputs: outputs/DSA_results.RData
## Outputs: outputs/dsa_plot.png

### Clear environment
rm(list=ls())

## Load required packages
library(ggplot2)

load("outputs/DSA_results.RData")
load("outputs/base_case_results.RData")

DSA <- as.data.table(DSA)
## relabel and remove duplicate rows
DSA[variable=="c.360.asc", Parameter := "Cost per B360S for ASC [15000 , 30000]"]
DSA[variable=="c.360.csc", Parameter := "Cost per B360S for CSC  [30000 , 60000]"]
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
  geom_linerange(aes(ymin = ymin, ymax = ymax), color = "lightskyblue", size = 5) +
  geom_point(aes(y = low.value), color = "red", size = 3) +
  geom_point(aes(y = high.value), color = "green", size = 3) +
  geom_hline(yintercept = base_case_value, color = "black") +
  coord_flip() +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Variable", y = "Net Monetary Benefit (NMB, £)",
       title = "Tornado Plot") +
  theme_minimal(base_size = 14) +  # Set a larger base font size
  theme(
    axis.text = element_text(size = 12),  # Larger axis text
    axis.title = element_text(size = 14),  # Larger axis title text
    plot.title = element_text(size = 16, face = "bold"),  # Larger plot title
    plot.margin = margin(10, 10, 10, 20)  # Extra space on the left margin for y-axis labels
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