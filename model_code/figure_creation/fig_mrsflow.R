###############################################
# TITLE: mRS state flow
### this script generates a plot of the flow of patients across mRS states over time
### inputs: outputs/mrstrace_results.csv
### outputs: outputs/mrsflow.png


### Clear environment
rm(list=ls())

### load libraries & data
library(tidyverse)

# Load your trace data
trace <- read.csv("outputs/mrstrace_results.csv")
# Ensure correct column types
trace <- trace %>%
  mutate(var = as.factor(var))  # convert 'var' to a factor for plotting

# Pivot data to long format for mRS states
trace_long <- trace %>%
  pivot_longer(cols = starts_with("mRS"), names_to = "State", values_to = "Proportion")

# Plot with facets for each value of `var`
mrs_plot <- ggplot(trace_long, aes(x = t, y = Proportion, color = State)) +
  geom_line(size = 1) +
  facet_wrap(~ var, ncol = 1) +
  labs(
    title = "Flow of Patients Across mRS States Over Time",
    x = "Cycle (t)",
    y = "Proportion of Cohort",
    color = "mRS State"
  ) +
  theme_minimal()

ggsave("outputs/mrs_plot.png", 
       plot = mrs_plot, 
       width = 15, 
       height = 14, 
       units = "in", 
       dpi = 300, 
       limitsize = FALSE)