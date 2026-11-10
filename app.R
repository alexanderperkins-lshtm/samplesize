library(shiny)
library(bslib)
library(rmarkdown)
library(tinytex)

# 1. UI Definition
ui <- page_navbar(
  title = "Sample Size Calculator",
  theme = bs_theme(version = 5, bootswatch = "zephyr"),
  
  # ---------- INDIVIDUAL RCT TAB ----------
  nav_panel("Individual RCT",
    layout_sidebar(
      sidebar = sidebar(
        width = 350,
        textInput("study_name_ind", "Study Name (Max 20 chars):", placeholder = "e.g., CARisSA Trial"),
        textAreaInput("notes_ind", "Study Notes (Max 200 chars):", placeholder = "Brief trial description or methodology notes...", rows = 3),
        hr(),
        tabsetPanel(
          tabPanel("Proportional",
                   br(),
                   textInput("eff_ind_prop", "Effect sizes (comma-separated):", "0.1, 0.2"),
                   textInput("rate_ind_prop", "Standard event rates (comma-separated):", "0.3, 0.5"),
                   numericInput("alpha_ind_prop", "Significance level (alpha):", value = 0.05, min = 0, max = 1, step = 0.01),
                   numericInput("power_ind_prop", "Power (1-beta):", value = 0.8, min = 0, max = 1, step = 0.05),
                   numericInput("loss_ind_prop", "Expected loss to follow-up (%):", value = 10, min = 0, max = 100),
                   actionButton("calc_ind_prop", "Calculate", class = "btn-primary w-100")
          ),
          tabPanel("Continuous",
                   br(),
                   textInput("diff_ind_cont", "Mean differences (comma-separated):", "5, 10"),
                   textInput("sd_ind_cont", "Standard deviations (comma-separated):", "10, 15"),
                   numericInput("alpha_ind_cont", "Significance level (alpha):", value = 0.05, min = 0, max = 1, step = 0.01),
                   numericInput("power_ind_cont", "Power (1-beta):", value = 0.8, min = 0, max = 1, step = 0.05),
                   numericInput("loss_ind_cont", "Expected loss to follow-up (%):", value = 10, min = 0, max = 100),
                   actionButton("calc_ind_cont", "Calculate", class = "btn-primary w-100")
          )
        )
      ),
      
      card(
        card_body(
          p("Calculates sample sizes for individually randomized controlled trials (1:1 allocation)."),
          p("Outputs represent the total sample size required across both arms, adjusted for loss to follow-up.")
        )
      ),
      
      card(
        card_header("Total Sample Size (Adjusted)"),
        tableOutput("table_ind_sample")
      ),
      
      uiOutput("dl_ui_ind")
    )
  ),
  
  # ---------- CLUSTER RCT TAB ----------
  nav_panel("Cluster RCT",
    layout_sidebar(
      sidebar = sidebar(
        width = 350,
        textInput("study_name_clu", "Study Name (Max 20 chars):", placeholder = "e.g., CARisSA Trial"),
        textAreaInput("notes_clu", "Study Notes (Max 200 chars):", placeholder = "Brief trial description or methodology notes...", rows = 3),
        hr(),
        tabsetPanel(
          tabPanel("Proportional",
                   br(),
                   textInput("eff_clu_prop", "Effect sizes (comma-separated):", "0.1, 0.5"),
                   textInput("rate_clu_prop", "Standard event rates (comma-separated):", "0.1, 0.5"),
                   numericInput("m_clu_prop", "Cluster size:", value = 10, min = 1),
                   numericInput("icc_clu_prop", "Intra-cluster correlation (ICC):", value = 0.01, min = 0, max = 1, step = 0.01),
                   numericInput("alpha_clu_prop", "Significance level (alpha):", value = 0.05, min = 0, max = 1, step = 0.01),
                   numericInput("power_clu_prop", "Power (1-beta):", value = 0.8, min = 0, max = 1, step = 0.05),
                   numericInput("loss_clu_prop", "Expected loss to follow-up (%):", value = 10, min = 0, max = 100),
                   actionButton("calc_clu_prop", "Calculate", class = "btn-primary w-100")
          ),
          tabPanel("Continuous",
                   br(),
                   textInput("diff_clu_cont", "Mean differences (comma-separated):", "5, 10"),
                   textInput("sd_clu_cont", "Standard deviations (comma-separated):", "10, 15"),
                   numericInput("m_clu_cont", "Cluster size:", value = 10, min = 1),
                   numericInput("icc_clu_cont", "Intra-cluster correlation (ICC):", value = 0.05, min = 0, max = 1, step = 0.01),
                   numericInput("alpha_clu_cont", "Significance level (alpha):", value = 0.05, min = 0, max = 1, step = 0.01),
                   numericInput("power_clu_cont", "Power (1-beta):", value = 0.8, min = 0, max = 1, step = 0.05),
                   numericInput("loss_clu_cont", "Expected loss to follow-up (%):", value = 10, min = 0, max = 100),
                   actionButton("calc_clu_cont", "Calculate", class = "btn-primary w-100")
          )
        )
      ),
      
      card(
        card_body(
          p("Calculates sample sizes for cluster randomized controlled trials (1:1 allocation)."),
          p("Outputs are adjusted for the design effect (ICC and cluster size) and expected loss to follow-up, ensuring an equal, whole number of clusters per arm.")
        )
      ),
      
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("Total Sample Size (Adjusted)"),
          tableOutput("table_clu_sample")
        ),
        card(
          card_header("Total Number of Clusters"),
          tableOutput("table_clu_cluster")
        )
      ),
      
      uiOutput("dl_ui_clu")
    )
  ),
  
  # ---------- METHODOLOGY TAB ----------
  nav_panel("Methodology",
    card(
      card_header("Statistical Methods & Variables"),
      card_body(
        markdown("
        ### Base Sample Size Calculation
        *   **Proportional Outcomes:** Calculated using the standard normal approximation for the difference between two proportions (`power.prop.test`).
        *   **Continuous Outcomes:** Calculated using a two-sample t-test framework (`power.t.test`).
        
        ### Cluster Adjustments
        For cluster RCTs, individuals within the same cluster are likely to have more similar outcomes than individuals in different clusters. This is quantified by the **Intra-cluster Correlation Coefficient (ICC)**. 
        
        To account for this loss of statistical efficiency, the base sample size is multiplied by the **Design Effect (DE)**:
        <br>*DE = 1 + (m - 1) &times; ICC* 
        <br>*(where **m** is the cluster size).*
        
        The required number of clusters is derived by dividing the adjusted sample size by the cluster size, rounding up (ceiling) to the nearest whole cluster per arm, and multiplying by two to ensure 1:1 trial balance.
        
        ### Loss to Follow-up
        The final sample size is inflated to account for expected dropout using the formula: 
        <br>*N_final = N_adjusted / (1 - Loss%)*
        
        ### Variable Definitions
        *   **Significance level (&alpha;):** The probability of a Type I error (rejecting the null hypothesis when it is true). Typically set at 0.05.
        *   **Power (1 - &beta;):** The probability of correctly rejecting the null hypothesis when a true effect exists. Typically set at 0.80 or 0.90.
        *   **Effect Size / Mean Difference:** The minimum clinically important difference you wish to detect between the two trial arms.
        *   **Event Rate / SD:** The expected baseline proportion of events (for binary outcomes) or the standard deviation (for continuous outcomes) in the control group.
        ")
      )
    )
  )
)

# 2. Server Logic
server <- function(input, output, session) {
  
  res <- reactiveValues(
    sample_table = NULL,
    cluster_table = NULL,
    conditions = NULL,
    trial_type = NULL,
    outcome_type = NULL,
    study_name = "Sample_Size",
    study_notes = ""
  )
  
  parse_inputs <- function(text_input) { as.numeric(trimws(unlist(strsplit(text_input, ",")))) }
  
  # --- Math Functions ---
  
  calc_ind_prop_fn <- Vectorize(function(eff, rate, alpha, pwr, loss) {
    n_per_arm <- power.prop.test(p1 = rate, p2 = rate + eff, sig.level = alpha, power = pwr)$n
    (n_per_arm) / (1 - loss / 100)
  })
  
  calc_ind_cont_fn <- Vectorize(function(diff, sd, alpha, pwr, loss) {
    n_per_arm <- power.t.test(delta = diff, sd = sd, sig.level = alpha, power = pwr, type = "two.sample")$n
    (n_per_arm) / (1 - loss / 100)
  })
  
  calc_clu_prop_fn <- Vectorize(function(eff, rate, m, icc, alpha, pwr, loss) {
    de <- 1 + (m - 1) * icc
    n_per_arm <- power.prop.test(p1 = rate, p2 = rate + eff, sig.level = alpha, power = pwr)$n
    (n_per_arm * de) / (1 - loss / 100)
  })
  
  calc_clu_cont_fn <- Vectorize(function(diff, sd, m, icc, alpha, pwr, loss) {
    de <- 1 + (m - 1) * icc
    n_per_arm <- power.t.test(delta = diff, sd = sd, sig.level = alpha, power = pwr, type = "two.sample")$n
    (n_per_arm * de) / (1 - loss / 100)
  })
  
  # --- Observers: Individual RCT ---
  
  observeEvent(input$calc_ind_prop, {
    effs <- parse_inputs(input$eff_ind_prop)
    rates <- parse_inputs(input$rate_ind_prop)
    n_adj_matrix <- outer(effs, rates, function(e, r) calc_ind_prop_fn(e, r, input$alpha_ind_prop, input$power_ind_prop, input$loss_ind_prop))
    samples <- ceiling(n_adj_matrix) * 2 
    
    # Store text fields securely (truncated to limits)
    raw_name <- substr(trimws(input$study_name_ind), 1, 20)
    res$study_name <- if (nchar(raw_name) > 0) raw_name else "Sample_Size"
    res$study_notes <- substr(trimws(input$notes_ind), 1, 200)
    
    res$sample_table <- matrix(samples, nrow = length(effs), dimnames = list(paste("Effect:", effs), paste("Rate:", rates)))
    res$cluster_table <- NULL
    res$trial_type <- "Individual"
    res$outcome_type <- "Proportional"
    res$conditions <- list("Effect Sizes" = input$eff_ind_prop, "Event Rates" = input$rate_ind_prop, "Alpha" = input$alpha_ind_prop, "Power" = input$power_ind_prop, "Loss to Follow-up (%)" = input$loss_ind_prop)
  })
  
  observeEvent(input$calc_ind_cont, {
    diffs <- parse_inputs(input$diff_ind_cont)
    sds <- parse_inputs(input$sd_ind_cont)
    n_adj_matrix <- outer(diffs, sds, function(d, s) calc_ind_cont_fn(d, s, input$alpha_ind_cont, input$power_ind_cont, input$loss_ind_cont))
    samples <- ceiling(n_adj_matrix) * 2 
    
    raw_name <- substr(trimws(input$study_name_ind), 1, 20)
    res$study_name <- if (nchar(raw_name) > 0) raw_name else "Sample_Size"
    res$study_notes <- substr(trimws(input$notes_ind), 1, 200)
    
    res$sample_table <- matrix(samples, nrow = length(diffs), dimnames = list(paste("Diff:", diffs), paste("SD:", sds)))
    res$cluster_table <- NULL
    res$trial_type <- "Individual"
    res$outcome_type <- "Continuous"
    res$conditions <- list("Mean Differences" = input$diff_ind_cont, "Standard Deviations" = input$sd_ind_cont, "Alpha" = input$alpha_ind_cont, "Power" = input$power_ind_cont, "Loss to Follow-up (%)" = input$loss_ind_cont)
  })
  
  # --- Observers: Cluster RCT ---
  
  observeEvent(input$calc_clu_prop, {
    effs <- parse_inputs(input$eff_clu_prop)
    rates <- parse_inputs(input$rate_clu_prop)
    n_adj_matrix <- outer(effs, rates, function(e, r) calc_clu_prop_fn(e, r, input$m_clu_prop, input$icc_clu_prop, input$alpha_clu_prop, input$power_clu_prop, input$loss_clu_prop))
    clusters <- ceiling(n_adj_matrix / input$m_clu_prop) * 2
    samples <- clusters * input$m_clu_prop
    
    raw_name <- substr(trimws(input$study_name_clu), 1, 20)
    res$study_name <- if (nchar(raw_name) > 0) raw_name else "Sample_Size"
    res$study_notes <- substr(trimws(input$notes_clu), 1, 200)
    
    res$sample_table <- matrix(samples, nrow = length(effs), dimnames = list(paste("Effect:", effs), paste("Rate:", rates)))
    res$cluster_table <- matrix(clusters, nrow = length(effs), dimnames = list(paste("Effect:", effs), paste("Rate:", rates)))
    res$trial_type <- "Cluster"
    res$outcome_type <- "Proportional"
    res$conditions <- list("Effect Sizes" = input$eff_clu_prop, "Event Rates" = input$rate_clu_prop, "Cluster Size" = input$m_clu_prop, "ICC" = input$icc_clu_prop, "Alpha" = input$alpha_clu_prop, "Power" = input$power_clu_prop, "Loss to Follow-up (%)" = input$loss_clu_prop)
  })
  
  observeEvent(input$calc_clu_cont, {
    diffs <- parse_inputs(input$diff_clu_cont)
    sds <- parse_inputs(input$sd_clu_cont)
    n_adj_matrix <- outer(diffs, sds, function(d, s) calc_clu_cont_fn(d, s, input$m_clu_cont, input$icc_clu_cont, input$alpha_clu_cont, input$power_clu_cont, input$loss_clu_cont))
    clusters <- ceiling(n_adj_matrix / input$m_clu_cont) * 2
    samples <- clusters * input$m_clu_cont
    
    raw_name <- substr(trimws(input$study_name_clu), 1, 20)
    res$study_name <- if (nchar(raw_name) > 0) raw_name else "Sample_Size"
    res$study_notes <- substr(trimws(input$notes_clu), 1, 200)
    
    res$sample_table <- matrix(samples, nrow = length(diffs), dimnames = list(paste("Diff:", diffs), paste("SD:", sds)))
    res$cluster_table <- matrix(clusters, nrow = length(diffs), dimnames = list(paste("Diff:", diffs), paste("SD:", sds)))
    res$trial_type <- "Cluster"
    res$outcome_type <- "Continuous"
    res$conditions <- list("Mean Differences" = input$diff_clu_cont, "Standard Deviations" = input$sd_clu_cont, "Cluster Size" = input$m_clu_cont, "ICC" = input$icc_clu_cont, "Alpha" = input$alpha_clu_cont, "Power" = input$power_clu_cont, "Loss to Follow-up (%)" = input$loss_clu_cont)
  })
  
  # --- Render Outputs ---
  
  output$table_ind_sample <- renderTable({ req(res$trial_type == "Individual"); res$sample_table }, rownames = TRUE, digits = 0)
  output$table_clu_sample <- renderTable({ req(res$trial_type == "Cluster"); res$sample_table }, rownames = TRUE, digits = 0)
  output$table_clu_cluster <- renderTable({ req(res$trial_type == "Cluster"); res$cluster_table }, rownames = TRUE, digits = 0)
  
  # --- Downloads ---
  
  output$dl_ui_ind <- renderUI({ req(res$trial_type == "Individual"); div(class = "d-flex justify-content-end", downloadButton("dl_pdf_ind", "Download PDF Report", class = "btn-success")) })
  output$dl_ui_clu <- renderUI({ req(res$trial_type == "Cluster"); div(class = "d-flex justify-content-end", downloadButton("dl_pdf_clu", "Download PDF Report", class = "btn-success")) })
  
  dl_handler_logic <- function(file) {
    req(res$sample_table)
    tempReport <- file.path(tempdir(), "report.Rmd")
    file.copy("report.Rmd", tempReport, overwrite = TRUE)
    
    rmarkdown::render(tempReport, 
                      output_file = file, 
                      output_format = "pdf_document", 
                      params = list(
                        study_name = res$study_name,
                        study_notes = res$study_notes,
                        sample_size_table = res$sample_table,
                        cluster_table = res$cluster_table,
                        conditions = res$conditions,
                        trial_type = res$trial_type,
                        outcome_type = res$outcome_type
                      ),
                      envir = new.env(parent = globalenv())
    )
  }
  
  # Inject the cleaned study name directly into the downloaded file name
  output$dl_pdf_ind <- downloadHandler(
    filename = function() { paste0(gsub("[^A-Za-z0-9_-]", "_", res$study_name), "_Report_", Sys.Date(), ".pdf") }, 
    content = dl_handler_logic
  )
  
  output$dl_pdf_clu <- downloadHandler(
    filename = function() { paste0(gsub("[^A-Za-z0-9_-]", "_", res$study_name), "_Report_", Sys.Date(), ".pdf") }, 
    content = dl_handler_logic
  )
}

shinyApp(ui = ui, server = server)