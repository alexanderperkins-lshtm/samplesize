# Deploy to shinyapps.io
rsconnect::deployApp(
  appDir = ".", 
  appName = "samplesize",
  appPrimaryDoc = "app.R"
)