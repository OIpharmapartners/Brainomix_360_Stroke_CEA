###############################################
# TITLE: transition state plot function
### This script is for plotting the states after AIS
### inputs: outputs/base_case_results.RData
### outputs: sankey_plot.html


### Clear environment
rm(list=ls())

## Load required packages

library(data.table)
library(networkD3)
library(htmlwidgets)
  
make_static_max_sankey <- function(trace_matrix, transition_matrix) {
  state_names <- colnames(trace_matrix)
  flows <- data.table()
  
  n_cycles <- dim(transition_matrix)[3]
  
  for (cycle in 1:n_cycles) {
    for (from_state in state_names) {
      for (to_state in state_names) {
        count <- trace_matrix[cycle, from_state] * transition_matrix[from_state, to_state, cycle]
        if (count > 0) {
          flows <- rbind(flows, data.table(source = from_state,
                                           target = to_state,
                                           value = round(count)))
        }
      }
    }
  }
  
  # Aggregate and take max for each transition
  max_flows <- flows[, .(value = max(value)), by = .(source, target)]
  
  # Remove flows from or to 'stroke' and 'no_AIS' or from MT/noMT to MT/noMT
  max_flows <- max_flows[!source %in% c("stroke", "no_AIS") & !target %in% c("stroke", "no_AIS")]
  max_flows <- max_flows[!source %in% c("MT", "NOMT")]
  
  # Create nodes
  nodes <- data.table(name = unique(c(max_flows$source, max_flows$target)))
  max_flows[, source_id := match(source, nodes$name) - 1]
  max_flows[, target_id := match(target, nodes$name) - 1]
  
  # Plot Sankey
  sankeyNetwork(Links = max_flows[, .(source = source_id, target = target_id, value)],
                Nodes = nodes,
                Source = "source",
                Target = "target",
                Value = "value",
                NodeID = "name",
                fontSize = 12,
                nodeWidth = 30)
}

load("outputs/base_case_results.RData")

standard <- make_static_max_sankey(base_case$trace.standard,
                       base_case$tm.standard)

standard

intervention <- make_static_max_sankey(base_case$trace.intervention,
                       base_case$tm.intervention)
intervention

# Save the plots as HTML files

saveWidget(standard, "outputs/standard_sankey_plot.html", selfcontained = TRUE)
saveWidget(intervention, "outputs/intervention_sankey_plot.html", selfcontained = TRUE)