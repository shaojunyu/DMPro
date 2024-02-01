library(shiny)
library(shinythemes)
library(tidyverse)
# probe <- head(epic, 100)
# saveRDS(probe, "probe.rds")
# Function to generate mock methylation data
# probe <- readRDS("probe.rds")
print(getwd())
generateMockData <- function(rows = 100) {
  set.seed(123)  # For reproducibility
  data.frame(
    SampleID = paste("Sample1"),
    ID = str_split(probe$ID, "_", simplify = T)[,1],
    # CpG_Island = sample(c("Island", "Shore", "Shelf", "OpenSea"), rows, replace = TRUE),
    CpG_Island = probe$Relation_to_UCSC_CpG_Island,
    Methylation_Level = rbeta(rows, 1, 2)
  )
}

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
      downloadButton('downloadSample', 'Download Sample Input', class = "btn-info",
                     icon = icon("download"))
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Status", icon = icon("info-circle"), textOutput("status")),
        tabPanel("Results Preview", icon = icon("table"),
                 div(class = "flex-container",
                     div(class = "download-button",
                         downloadButton('downloadData', 'Download Results', class = "btn-primary",
                                        icon = icon("file-download"))),
                     tableOutput("dataPreview")
                 )
        ),
        id = "mainTabset"
      )
    )
  )
)

# Define server logic
server <- function(input, output) {
  output$status <- renderText({
    if (is.null(input$file1)) {
      return("Awaiting data...")
    }
    paste("Data uploaded. Processing for", input$tissueType, "with", input$diseaseStatus, "status.")
  })
  
  output$dataPreview <- renderTable({
    # if (is.null(input$file1)) {
    #   return(generateMockData())
    # }
    return(generateMockData())
    # Placeholder for actual file processing
  }, digits = 4, hover = T)
  
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
      paste("sample_data_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(generateMockData(), file, row.names = FALSE)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
