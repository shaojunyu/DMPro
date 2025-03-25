library(shiny)
library(shinythemes)
library(tidyverse)
library(DT)

# probe <- head(epic, 100)
# saveRDS(probe, "probe.rds")
# Function to generate mock methylation data
# probe <- readRDS("probe.rds")
print(getwd())

# Define UI with icons and custom CSS
ui <- fluidPage(
  theme = shinytheme("flatly"),
  titlePanel("DNA Methylation Imputation", windowTitle = "DNA Methylation Imputation"),
  
  tags$head(
    tags$style(HTML("
            .sidebar .form-group {
                margin-bottom: 15px;
            }
            .flex-container {
                display: flex;
                flex-direction: column;
                align-items: flex-start;
            }
            .flex-container .download-button {
                align-self: flex-end;
                margin-bottom: 10px;
            }
        "))
  ),
  
  sidebarLayout(
    sidebarPanel(
      fileInput('file1', 'Upload DNA Data', 
                buttonLabel = tags$span(icon("upload"), " Upload")),
      selectInput("tissueType", "Select Tissue Type", 
                  choices = c("Whole blood", "Peripheral blood", "Cortex")),
      selectInput("diseaseStatus", "Disease Status", 
                  choices = c("NA" = "Healthy", "Diseased" = "Cancer")),
      # start imputation button, enabled when file is uploaded
      shinyjs::useShinyjs(),
      actionButton("Start imputation", "Start Imputation",
                   class = "btn-primary", disabled = T),
      # placeholder for imputation button
      
      
      
      # add vertical space
      tags$hr(),
      downloadButton('downloadSample', 'Download Sample Input', class = "btn-info",
                     icon = icon("download"))
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Status", icon = icon("info-circle"), textOutput("status")),
        tabPanel("Uploaded Data", icon = icon("file"), 
                 DT::dataTableOutput("dataPreview")),
        tabPanel("Results Preview", icon = icon("table"),
                 div(class = "flex-container",
                     div(class = "download-button",
                         downloadButton('downloadData', 'Download Results', class = "btn-primary",
                                        icon = icon("file-download"))),
                     tableOutput("resultPreview")
                 )
        ),
        id = "mainTabset"
      )
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  # check if file is uploaded, enable imputation button
  observe({
    shinyjs::toggleState("Start imputation", !is.null(input$file1))
  })
  
  # event handler for imputation button
  observeEvent(input[["Start imputation"]], {
    # Placeholder for imputation logic
    print("Imputation started")
  })
  
  output$status <- renderText({
    if (is.null(input$file1)) {
      return("Awaiting data...")
    }
    paste("Data uploaded. Processing for", 
          input$tissueType, "with", 
          input$diseaseStatus, "status.",
          input$file1$datapath)
  })
  
  output$dataPreview <- renderDataTable({
    if (is.null(input$file1)) {
      return(NULL)
    }
    read.table(input$file1$datapath, sep = ",", header = T)
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("imputed_data_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      # Placeholder for data export logic
    }
  )
  
  output$downloadSample <- downloadHandler(
    filename = function() {
      "HM450K_chr20_test_data.csv"
    },
    content = function(file) {
      file.copy("test_data_for_shiny.csv", file)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
