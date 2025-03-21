library(shiny)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(DT)
library(shinythemes)  # For better themes

# Define UI for application
ui <- fluidPage(
  
  # Apply a modern theme
  theme = shinytheme("flatly"),
  
  # Application title
  titlePanel("Pathway Enrichment Analysis", windowTitle = "Pathway Enrichment"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    
    # Sidebar panel for inputs
    sidebarPanel(
      h4("Input Parameters", style = "color: #2C3E50;"),
      
      # Input: Upload a file or paste gene list
      textAreaInput("gene_list", "Paste gene list (one gene per line):", 
                    rows = 10, placeholder = "Enter gene IDs (e.g., 7157, 7422)"),
      
      # Input: Select organism database
      selectInput("organism_db", "Select organism database:",
                  choices = c("Human" = "org.Hs.eg.db",
                              "Mouse" = "org.Mm.eg.db"),
                  selected = "org.Hs.eg.db"),
      
      # Input: Select ontology
      selectInput("ontology", "Select ontology:",
                  choices = c("Biological Process" = "BP",
                              "Molecular Function" = "MF",
                              "Cellular Component" = "CC"),
                  selected = "BP"),
      
      # Action button to run analysis
      actionButton("run_analysis", "Run Enrichment Analysis", 
                   class = "btn-primary", icon = icon("play")),
      
      # Add some spacing
      br(),
      br(),
      
      # Add a note for users
      p("Note: Enter gene IDs in Entrez format (e.g., 7157 for TP53).", 
        style = "font-size: 12px; color: #7F8C8D;")
    ),
    
    # Main panel for displaying outputs
    mainPanel(
      tabsetPanel(
        # Tab 1: Enrichment Plot
        tabPanel("Enrichment Plot", 
                 plotOutput("enrichment_plot", height = "600px")),
        
        # Tab 2: Enrichment Table
        tabPanel("Enrichment Table", 
                 DTOutput("enrichment_table"))
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Reactive expression to perform enrichment analysis
  enrichment_result <- eventReactive(input$run_analysis, {
    req(input$gene_list)
    genes <- unlist(strsplit(input$gene_list, "\n"))
    genes <- as.character(genes)  # Ensure genes are in character format
    enrichGO(genes, OrgDb = input$organism_db, ont = input$ontology)
  })
  
  # Output: Enrichment plot
  output$enrichment_plot <- renderPlot({
    req(enrichment_result())
    dotplot(enrichment_result(), showCategory = 20, font.size = 12) +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # Output: Enrichment results table
  output$enrichment_table <- renderDT({
    req(enrichment_result())
    datatable(as.data.frame(enrichment_result()), 
              options = list(pageLength = 10, scrollX = TRUE),
              class = "display nowrap")
  })
}

# Run the application 
shinyApp(ui = ui, server = server)