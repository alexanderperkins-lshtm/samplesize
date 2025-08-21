# Sample Size Calculation for Cluster RCTs using Design Effect Method

# Prompt user for input values
cat("Enter the range of effect sizes (e.g., 0.1,0.5): ")
effect_sizes <- as.numeric(unlist(strsplit(readline(), ",")))

cat("Enter the range of standard event rates (e.g., 0.1,0.5): ")
standard_event_rates <- as.numeric(unlist(strsplit(readline(), ",")))

cat("Enter the cluster size: ")
cluster_size <- as.numeric(readline())

cat("Enter the intra-cluster correlation coefficient (ICC): ")
icc <- as.numeric(readline())

cat("Enter the significance level (alpha, e.g., 0.05): ")
alpha <- as.numeric(readline())

cat("Enter the power (1-beta, e.g., 0.8): ")
power <- as.numeric(readline())

# Function to calculate sample size for a given effect size and event rate
calculate_sample_size <- function(effect_size, event_rate, cluster_size, icc, alpha, power) {
  # Design effect
  design_effect <- 1 + (cluster_size - 1) * icc

  # Calculate sample size using power.prop.test
  result <- power.prop.test(p1 = event_rate, p2 = event_rate + effect_size, sig.level = alpha, power = power)

  # Adjust for design effect
  adjusted_sample_size <- ceiling(result$n * design_effect)

  return(adjusted_sample_size)
}

# Generate table of sample sizes
sample_size_table <- matrix(nrow = length(effect_sizes), ncol = length(standard_event_rates))
rownames(sample_size_table) <- paste("Effect Size:", effect_sizes)
colnames(sample_size_table) <- paste("Event Rate:", standard_event_rates)

for (i in seq_along(effect_sizes)) {
  for (j in seq_along(standard_event_rates)) {
    sample_size_table[i, j] <- calculate_sample_size(effect_sizes[i], standard_event_rates[j], cluster_size, icc, alpha, power)
  }
}

# Print the table as HTML
cat("\nSample Size Table (HTML):\n")
html_table <- paste0(
  "<table border='1'>",
  "<tr><th>Effect Size/Event Rate</th>",
  paste0("<th>", colnames(sample_size_table), "</th>", collapse = ""),
  "</tr>",
  paste(apply(sample_size_table, 1, function(row, effect_size) {
    paste0(
      "<tr><td>", effect_size, "</td>",
      paste0("<td>", row, "</td>", collapse = ""),
      "</tr>"
    )
  }, effect_size = rownames(sample_size_table)), collapse = ""),
  "</table>"
)
cat(html_table)

# Optionally, save the HTML table to a file
writeLines(html_table, "sample_size_table.html")
cat("\nHTML table saved to 'sample_size_table.html'\n")
