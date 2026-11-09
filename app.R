library(shiny)
library(bslib)
library(rmarkdown)

# 1. Modernized UI using bslib
ui <- page_sidebar(
  title = "Sample Size Calculation for Cluster RCTs",
  theme = bs_theme(version = 5, bootswatch = "zephyr"),
  
  sidebar = sidebar(
    width = 350,
    tabsetPanel(
      id = "calc_tabs",
      
      tabPanel("Proportional",
               br(),
               textInput("effect_sizes", "Effect sizes (comma-separated):", "0.1, 0.5"),
               textInput("event_rates", "Standard event rates (comma-separated):", "0.1, 0.5"),
               numericInput("cluster_size_cat", "Cluster size:", value = 10, min = 1),
               numericInput("icc_cat", "Intra-cluster correlation (ICC):", value = 0.01, min = 0, max = 1, step = 0.01),
               numericInput("alpha_cat", "Significance level (alpha):", value = 0.05, min = 0, max = 1, step = 0.01),
               numericInput("power_cat", "Power (1-beta):", value = 0.8, min = 0, max = 1, step = 0.05),
               numericInput("loss_cat", "Expected loss to follow-up (%):", value = 10, min = 0, max = 100),
               actionButton("calculate_cat", "Calculate", class = "btn-primary w-100")
      ),
      
      tabPanel("Continuous",
               br(),
               textInput("mean_diff_range", "Mean differences (comma-separated):", "5, 10"),
               textInput("sd_range", "Standard deviations (comma-separated):", "10, 15"),
               numericInput("cluster_size_cont", "Cluster size:", value = 10, min = 1),
               numericInput("icc_cont", "Intra-cluster correlation (ICC):", value = 0.05, min = 0, max = 1, step = 0.01),
               numericInput("alpha_cont", "Significance level (alpha):", value = 0.05, min = 0, max = 1, step = 0.01),
               numericInput("power_cont", "Power (1-beta):", value = 0.8, min = 0, max = 1, step = 0.05),
               numericInput("loss_cont", "Expected loss to follow-up (%):", value = 10, min = 0, max = 100),
               actionButton("calculate_cont", "Calculate", class = "btn-primary w-100")
      )
    )
  ),
  
  # Introduction text
  card(
    card_body(
      p("This application calculates a range of sample sizes for proportional or continuous outcomes for a cluster RCT. It also calculates the corresponding number of clusters required."),
      p("The output shows the total sample size for a two-arm trial and includes the loss to follow-up and clustering adjustments. For questions or comments, please contact Alexander.Perkins@LSHTM.ac.uk. Figures are advisory and should be checked with a qualified statistician before use. No responsibility is taken for their use.")
    )
  ),
  
  # Results cards
  layout_columns(
    col_widths = c(6, 6),
    card(
      card_header("Total Sample Size (Adjusted)"),
      tableOutput("sample_size_table")
    ),
    card(
      card_header("Total Number of Clusters"),
      tableOutput("cluster_table")
    )
  ),
  
  uiOutput("download_ui")
)


# 2. Refactored Server Logic
server <- function(input, output, session) {
  
  # Centralized reactive values to hold state
  res <- reactiveValues(
    sample_table = NULL,
    cluster_table = NULL,
    conditions = NULL,
    type = NULL
  )
  
  # Helper: Parse comma-separated inputs safely
  parse_inputs <- function(text_input) {
    as.numeric(trimws(unlist(strsplit(text_input, ","))))
  }
  
  # Vectorized calculation functions
  calc_prop <- Vectorize(function(eff, rate, m, icc, alpha, pwr, loss) {
    de <- 1 + (m - 1) * icc
    n_per_arm_ind <- power.prop.test(p1 = rate, p2 = rate + eff, sig.level = alpha, power = pwr)$n
    n_per_arm_adj <- (n_per_arm_ind * de) / (1 - loss / 100)
    
    clusters_per_arm <- ceiling(n_per_arm_adj / m)
    total_clusters <- clusters_per_arm * 2
    total_n <- total_clusters * m 
    
    c(total_n, total_clusters)
  })
  
  calc_cont <- Vectorize(function(diff, sd, m, icc, alpha, pwr, loss) {
    de <- 1 + (m - 1) * icc
    n_per_arm_ind <- power.t.test(delta = diff, sd = sd, sig.level = alpha, power = pwr, type = "two.sample")$n
    n_per_arm_adj <- (n_per_arm_ind * de) / (1 - loss / 100)
    
    clusters_per_arm <- ceiling(n_per_arm_adj / m)
    total_clusters <- clusters_per_arm * 2
    total_n <- total_clusters * m
    
    c(total_n, total_clusters)
  })
  
  # Proportional Trigger
  observeEvent(input$calculate_cat, {
    effs <- parse_inputs(input$effect_sizes)
    rates <- parse_inputs(input$event_rates)
    
    grid <- outer(effs, rates, function(e, r) calc_prop(e, r, input$cluster_size_cat, input$icc_cat, input$alpha_cat, input$power_cat, input$loss_cat))
    
    # Extract matrices from the 3D array returned by outer + Vectorize
    res$sample_table <- matrix(grid[1,,], nrow = length(effs), dimnames = list(paste("Effect:", effs), paste("Rate:", rates)))
    res$cluster_table <- matrix(grid[2,,], nrow = length(effs), dimnames = list(paste("Effect:", effs), paste("Rate:", rates)))
    
    res$type <- "Proportional"
    res$conditions <- list(
      "Effect Sizes" = input$effect_sizes,
      "Event Rates" = input$event_rates,
      "Cluster Size" = input$cluster_size_cat,
      "ICC" = input$icc_cat,
      "Alpha" = input$alpha_cat,
      "Power" = input$power_cat,
      "Loss to Follow-up (%)" = input$loss_cat
    )
  })
  
  # Continuous Trigger
  observeEvent(input$calculate_cont, {
    diffs <- parse_inputs(input$mean_diff_range)
    sds <- parse_inputs(input$sd_range)
    
    grid <- outer(diffs, sds, function(d, s) calc_cont(d, s, input$cluster_size_cont, input$icc_cont, input$alpha_cont, input$power_cont, input$loss_cont))
    
    res$sample_table <- matrix(grid[1,,], nrow = length(diffs), dimnames = list(paste("Diff:", diffs), paste("SD:", sds)))
    res$cluster_table <- matrix(grid[2,,], nrow = length(diffs), dimnames = list(paste("Diff:", diffs), paste("SD:", sds)))
    
    res$type <- "Continuous"
    res$conditions <- list(
      "Mean Differences" = input$mean_diff_range,
      "Standard Deviations" = input$sd_range,
      "Cluster Size" = input$cluster_size_cont,
      "ICC" = input$icc_cont,
      "Alpha" = input$alpha_cont,
      "Power" = input$power_cont,
      "Loss to Follow-up (%)" = input$loss_cont
    )
  })
  
  # Render Tables (only when data exists)
  output$sample_size_table <- renderTable({ req(res$sample_table); res$sample_table }, rownames = TRUE, digits = 0)
  output$cluster_table <- renderTable({ req(res$cluster_table); res$cluster_table }, rownames = TRUE, digits = 0)
  
  # Render Download Button conditionally
  output$download_ui <- renderUI({
    req(res$sample_table)
    div(class = "d-flex justify-content-end",
        downloadButton("download_report", "Download Report", class = "btn-success")
    )
  })
  
  # Extracted top-level download handler
  output$download_report <- downloadHandler(
    filename = function() { paste0("sample_size_report_", Sys.Date(), ".html") },
    content = function(file) {
      req(res$sample_table)
      rmarkdown::render("report.Rmd", 
                        output_file = file, 
                        output_format = "html_document", 
                        params = list(
                          sample_size_table = res$sample_table,
                          cluster_table = res$cluster_table,
                          conditions = res$conditions,
                          calculation_type = res$type
                        ))
    }
  )
}

shinyApp(ui = ui, server = server)