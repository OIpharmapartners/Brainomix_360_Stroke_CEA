#######code thats used to check the distributions of PSA parameters looks reasonable
#### not used for manuscript/reports


library(data.table)

load("inputs/created_inputs/params.psa.sample.RData")
# Remove 'sample_id' if it's the first column
param_data <- params.psa.sample[, !("sample_id"), with = FALSE]

# Loop over each column and plot
for (param in names(param_data)) {
  hist(param_data[[param]],
       main = param,
       xlab = "",
       col = "grey",
       border = "white")
  
  # Optional: pause between plots until user presses a key
  readline(prompt = "Press [enter] to continue to next plot...")
}