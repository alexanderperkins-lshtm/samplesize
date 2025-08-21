library(shiny)
setwd("/Users/alexanderperkins/Library/CloudStorage/OneDrive-LondonSchoolofHygieneandTropicalMedicine/samplesize")

# Define UI for the app
ui <- fluidPage(
  titlePanel("Sample Size Calculation for Cluster RCTs"),

  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        tabPanel("Categorical Proportions",
          textInput("effect_sizes", "Enter range of effect sizes (e.g., 0.1,0.5):", "0.1,0.5"),
          textInput("event_rates", "Enter range of standard event rates (e.g., 0.1,0.5):", "0.1,0.5"),
          numericInput("cluster_size_cat", "Enter cluster size:", value = 10, min = 1),
          numericInput("icc_cat", "Enter intra-cluster correlation coefficient (ICC):", value = 0.01, min = 0, max = 1),
          numericInput("alpha_cat", "Enter significance level (alpha):", value = 0.05, min = 0, max = 1),
          numericInput("power_cat", "Enter power (1-beta):", value = 0.8, min = 0, max = 1),
          numericInput("loss_cat", "Enter expected loss to follow-up (%):", value = 10, min = 0, max = 100),
          actionButton("calculate_cat", "Calculate")
        ),
        tabPanel("Continuous Variables",
          textInput("mean_diff_range", "Enter range of mean differences (e.g., 5,10):", "5,10"),
          textInput("sd_range", "Enter range of standard deviations (e.g., 10,15):", "10,15"),
          numericInput("cluster_size_cont", "Enter cluster size:", value = 10, min = 1),
          numericInput("icc_cont", "Enter intra-cluster correlation coefficient (ICC):", value = 0.01, min = 0, max = 1),
          numericInput("alpha_cont", "Enter significance level (alpha):", value = 0.05, min = 0, max = 1),
          numericInput("power_cont", "Enter power (1-beta):", value = 0.8, min = 0, max = 1),
          numericInput("loss_cont", "Enter expected loss to follow-up (%):", value = 10, min = 0, max = 100),
          actionButton("calculate_cont", "Calculate")
        )
      )
    ),

    mainPanel(
      h3("Sample Size Table"),
      tableOutput("sample_size_table"),
      h3("Number of Clusters Table"),
      tableOutput("cluster_table"),
      h3("Rounded Sample Size Table"),
      tableOutput("rounded_table"),
      downloadButton("download_table", "Download CSV Table")
    )
  )
)

# Define server logic
server <- function(input, output) {
  calculate_sample_size <- function(effect_size, event_rate, cluster_size, icc, alpha, power) {
    design_effect <- 1 + (cluster_size - 1) * icc
    result <- power.prop.test(p1 = event_rate, p2 = event_rate + effect_size, sig.level = alpha, power = power)
    adjusted_sample_size <- ceiling(result$n * design_effect)
    return(adjusted_sample_size)
  }

  calculate_sample_size_continuous <- function(mean_diff, sd, cluster_size, icc, alpha, power) {
    design_effect <- 1 + (cluster_size - 1) * icc
    result <- power.t.test(delta = mean_diff, sd = sd, sig.level = alpha, power = power, type = "two.sample")
    adjusted_sample_size <- ceiling(result$n * design_effect)
    return(adjusted_sample_size)
  }

  observeEvent(input$calculate_cat, {
    effect_sizes <- as.numeric(unlist(strsplit(input$effect_sizes, ",")))
    event_rates <- as.numeric(unlist(strsplit(input$event_rates, ",")))

    sample_size_table <- matrix(nrow = length(effect_sizes), ncol = length(event_rates))
    cluster_table <- matrix(nrow = length(effect_sizes), ncol = length(event_rates))
    rounded_table <- matrix(nrow = length(effect_sizes), ncol = length(event_rates))
    rownames(sample_size_table) <- paste("Effect Size:", effect_sizes)
    colnames(sample_size_table) <- paste("Event Rate:", event_rates)
    rownames(cluster_table) <- paste("Effect Size:", effect_sizes)
    colnames(cluster_table) <- paste("Event Rate:", event_rates)
    rownames(rounded_table) <- paste("Effect Size:", effect_sizes)
    colnames(rounded_table) <- paste("Event Rate:", event_rates)

    for (i in seq_along(effect_sizes)) {
      for (j in seq_along(event_rates)) {
        base_sample_size <- 2 * calculate_sample_size(effect_sizes[i], event_rates[j], input$cluster_size_cat, input$icc_cat, input$alpha_cat, input$power_cat)
        adjusted_sample_size <- ceiling(base_sample_size * (1 + input$loss_cat / 100))
        sample_size_table[i, j] <- adjusted_sample_size
        cluster_table[i, j] <- ceiling(adjusted_sample_size / input$cluster_size_cat)
        rounded_table[i, j] <- ifelse(input$cluster_size_cat > 0, 
                                      ceiling(adjusted_sample_size / input$cluster_size_cat) * input$cluster_size_cat, 
                                      NA)
      }
    }

    output$sample_size_table <- renderTable({
      sample_size_table
    }, rownames = TRUE)

    output$cluster_table <- renderTable({
      cluster_table
    }, rownames = TRUE)

    output$rounded_table <- renderTable({
      rounded_table
    }, rownames = TRUE)

    output$download_table <- downloadHandler(
      filename = function() {
        "sample_size_and_clusters.csv"
      },
      content = function(file) {
        combined_table <- cbind(Sample_Size = sample_size_table, Number_of_Clusters = cluster_table, Rounded_Sample_Size = rounded_table)
        combined_table <- rbind(combined_table, Note = "Figures are for a two-arm trial with adjustments included.")
        write.csv(combined_table, file, row.names = TRUE)
      }
    )
  })

  observeEvent(input$calculate_cont, {
    mean_diffs <- as.numeric(unlist(strsplit(input$mean_diff_range, ",")))
    sds <- as.numeric(unlist(strsplit(input$sd_range, ",")))

    sample_size_table <- matrix(nrow = length(mean_diffs), ncol = length(sds))
    cluster_table <- matrix(nrow = length(mean_diffs), ncol = length(sds))
    rownames(sample_size_table) <- paste("Mean Difference:", mean_diffs)
    colnames(sample_size_table) <- paste("SD:", sds)
    rownames(cluster_table) <- paste("Mean Difference:", mean_diffs)
    colnames(cluster_table) <- paste("SD:", sds)

    for (i in seq_along(mean_diffs)) {
      for (j in seq_along(sds)) {
        base_sample_size <- 2 * calculate_sample_size_continuous(mean_diffs[i], sds[j], input$cluster_size_cont, input$icc_cont, input$alpha_cont, input$power_cont)
        adjusted_sample_size <- ceiling(base_sample_size * (1 + input$loss_cont / 100))
        sample_size_table[i, j] <- adjusted_sample_size
        cluster_table[i, j] <- ceiling(adjusted_sample_size / input$cluster_size_cont)
      }
    }

    output$sample_size_table <- renderTable({
      sample_size_table
    }, rownames = TRUE)

    output$cluster_table <- renderTable({
      cluster_table
    }, rownames = TRUE)

    output$download_table <- downloadHandler(
      filename = function() {
        "sample_size_and_clusters.csv"
      },
      content = function(file) {
        combined_table <- cbind(Sample_Size = sample_size_table, Number_of_Clusters = cluster_table)
        combined_table <- rbind(combined_table, Note = "Figures are for a two-arm trial with adjustments included.")
        write.csv(combined_table, file, row.names = TRUE)
      }
    )
  })
}

# Run the application
shinyApp(ui = ui, server = server)
