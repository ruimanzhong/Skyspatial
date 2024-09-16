library(shiny)
library(ggplot2)

# Define UI for application
ui <- fluidPage(
  titlePanel("CSV File Upload and Plotting"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Choose CSV File",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")),
      tags$hr(),
      checkboxInput("header", "Header", TRUE),
      radioButtons("sep", "Separator",
                   choices = c(Comma = ",",
                               Semicolon = ";",
                               Tab = "\t"),
                   selected = ",")
    ),
    
    mainPanel(
      plotOutput("plot")
    )
  )
)

# Define server logic
server <- function(input, output) {
  output$plot <- renderPlot({
    req(input$file1)
    
    # Read the uploaded file
    inFile <- input$file1
    
    # Read the CSV file
    data <- read.csv(inFile$datapath,
                     header = input$header,
                     sep = input$sep)
    
    # Check if the required columns exist
    if (!all(c("x", "y") %in% names(data))) {
      stop("CSV must contain 'x' and 'y' columns.")
    }
    
    # Plot the data
    ggplot(data, aes(x = x, y = y)) +
      geom_point() +
      theme_minimal() +
      labs(title = "Scatter Plot of CSV Data", x = "X-axis", y = "Y-axis")
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
