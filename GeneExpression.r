install.packages("pheatmap")
# Load necessary libraries
library(shiny)
library(ggplot2)
library(pheatmap)
library(DT)

# Define UI for application
ui <- fluidPage(
  
  # Application title
  titlePanel("Gene Expression Data Explorer"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    
    # Sidebar panel for inputs
    sidebarPanel(
      
      # Input: Upload a file
      fileInput("file", "Upload gene expression data (CSV file)",
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      
      # Input: Filter genes by expression level
      sliderInput("expr_filter", "Filter genes by expression level:",
                  min = 0, max = 100, value = c(0, 100)),
      
      # Input: Select samples for heatmap
      selectInput("samples", "Select samples for heatmap:",
                  choices = NULL, multiple = TRUE),
      
      # Input: Select gene for boxplot
      selectInput("gene_boxplot", "Select gene for boxplot:",
                  choices = NULL)
    ),
    
    # Main panel for displaying outputs
    mainPanel(
      
      # Output: Heatmap
      plotOutput("heatmap"),
      
      # Output: Boxplot
      plotOutput("boxplot"),
      
      # Output: Filtered data table
      DTOutput("filtered_table")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Reactive expression to read the uploaded file
  data <- reactive({
    req(input$file)
    read.csv(input$file$datapath, row.names = 1)
  })
  
  # Update selectInput choices based on the uploaded dataset
  observeEvent(data(), {
    updateSelectInput(session, "samples", choices = colnames(data()))
    updateSelectInput(session, "gene_boxplot", choices = rownames(data()))
  })
  
  # Output: Heatmap
  output$heatmap <- renderPlot({
    req(input$samples)
    filtered_data <- data()[rowMeans(data()) >= input$expr_filter[1] & rowMeans(data()) <= input$expr_filter[2], input$samples]
    pheatmap(filtered_data, scale = "row", clustering_distance_rows = "euclidean")
  })
  
  # Output: Boxplot
  output$boxplot <- renderPlot({
    req(input$gene_boxplot)
    gene_data <- data()[input$gene_boxplot, ]
    ggplot(data.frame(Sample = colnames(data()), Expression = as.numeric(gene_data)), 
           aes(x = Sample, y = Expression)) +
      geom_boxplot(fill = "orange") +
      theme_minimal() +
      labs(title = paste("Expression of", input$gene_boxplot))
  })
  
  # Output: Filtered data table
  output$filtered_table <- renderDT({
    req(input$expr_filter)
    filtered_data <- data()[rowMeans(data()) >= input$expr_filter[1] & rowMeans(data()) <= input$expr_filter[2], ]
    datatable(filtered_data)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)