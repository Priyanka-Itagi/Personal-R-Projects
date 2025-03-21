# Load necessary libraries
library(shiny)
library(ggplot2)
library(dplyr)

# Define UI for application
ui <- fluidPage(
  
  # Application title
  titlePanel("Biostatistics Data Analysis"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    
    # Sidebar panel for inputs
    sidebarPanel(
      
      # Input: Upload a file
      fileInput("file", "Upload your dataset (CSV file)",
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      
      # Input: Select variable for summary statistics
      selectInput("var_summary", "Select variable for summary statistics:",
                  choices = NULL),
      
      # Input: Select variables for scatter plot
      selectInput("var_x", "Select X variable for scatter plot:",
                  choices = NULL),
      selectInput("var_y", "Select Y variable for scatter plot:",
                  choices = NULL),
      
      # Input: Checkbox for linear regression line
      checkboxInput("reg_line", "Add regression line", value = FALSE),
      
      # Input: Select variable for histogram
      selectInput("var_hist", "Select variable for histogram:",
                  choices = NULL),
      
      # Input: Select variable for boxplot
      selectInput("var_box", "Select variable for boxplot:",
                  choices = NULL),
      
      # Input: Grouping variable for boxplot
      selectInput("var_group", "Select grouping variable for boxplot (optional):",
                  choices = NULL)
    ),
    
    # Main panel for displaying outputs
    mainPanel(
      
      # Output: Summary statistics
      verbatimTextOutput("summary"),
      
      # Output: Scatter plot
      plotOutput("scatter_plot"),
      
      # Output: Histogram
      plotOutput("histogram"),
      
      # Output: Boxplot
      plotOutput("boxplot"),
      
      # Output: Error messages
      textOutput("error_message")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Reactive expression to read the uploaded file
  data <- reactive({
    req(input$file)
    tryCatch({
      read.csv(input$file$datapath)
    }, error = function(e) {
      showNotification("Error: Invalid file format. Please upload a valid CSV file.", type = "error")
      return(NULL)
    })
  })
  
  # Update selectInput choices based on the uploaded dataset
  observeEvent(data(), {
    updateSelectInput(session, "var_summary", choices = names(data()))
    updateSelectInput(session, "var_x", choices = names(data()))
    updateSelectInput(session, "var_y", choices = names(data()))
    updateSelectInput(session, "var_hist", choices = names(data()))
    updateSelectInput(session, "var_box", choices = names(data()))
    updateSelectInput(session, "var_group", choices = c("None", names(data())))
  })
  
  # Output: Summary statistics
  output$summary <- renderPrint({
    req(input$var_summary)
    summary(data()[[input$var_summary]])
  })
  
  # Output: Scatter plot
  output$scatter_plot <- renderPlot({
    req(input$var_x, input$var_y)
    tryCatch({
      p <- ggplot(data(), aes_string(x = input$var_x, y = input$var_y)) +
        geom_point() +
        theme_minimal()
      
      if (input$reg_line) {
        p <- p + geom_smooth(method = "lm", se = FALSE, color = "red")
      }
      
      p
    }, error = function(e) {
      output$error_message <- renderText("Error: Invalid variables selected for scatter plot.")
      return(NULL)
    })
  })
  
  # Output: Histogram
  output$histogram <- renderPlot({
    req(input$var_hist)
    tryCatch({
      ggplot(data(), aes_string(x = input$var_hist)) +
        geom_histogram(binwidth = 30, fill = "yellow", color = "black") +
        theme_minimal() +
        labs(title = paste("Histogram of", input$var_hist))
    }, error = function(e) {
      output$error_message <- renderText("Error: Invalid variable selected for histogram.")
      return(NULL)
    })
  })
  
  # Output: Boxplot
  output$boxplot <- renderPlot({
    req(input$var_box)
    tryCatch({
      if (input$var_group == "None") {
        ggplot(data(), aes_string(y = input$var_box)) +
          geom_boxplot(fill = "orange") +
          theme_minimal() +
          labs(title = paste("Boxplot of", input$var_box))
      } else {
        ggplot(data(), aes_string(x = input$var_group, y = input$var_box)) +
          geom_boxplot(fill = "orange") +
          theme_minimal() +
          labs(title = paste("Boxplot of", input$var_box, "by", input$var_group))
      }
    }, error = function(e) {
      output$error_message <- renderText("Error: Invalid variables selected for boxplot.")
      return(NULL)
    })
  })
  
  # Output: Error messages
  output$error_message <- renderText({
    NULL
  })
}

# Run the application 
shinyApp(ui = ui, server = server)