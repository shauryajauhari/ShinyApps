if(interactive()){

  ui <- fluidPage(
  titlePanel(title = "Welcome to the Benchmarking of GSA tools."),
  br(),
  h5("This application allows analyzing a user defined dataset for ranking various enrichment tools for genomic regions.
     The current suite of tools include GREAT, Enrichr, Chipenrich, Broadenrich, Seq2pathway. However, only the latter 3 are available
     as functions in R and hence form a part of this interface. The user is required to select the dataset, tool, and the metric of 
     comparison. The application returns the plot and values as a result."),
  br(),
  sidebarLayout(
    
    sidebarPanel(
      selectInput("m", "Select a comparison metric",
                   choices = c("Sensitivity"='sn', "Specificity"='sp',  "Prioritization"='pn', "Precision"='pr')),
      br(),
      radioButtons("t", "Select a tool",
                   choices = c("Chipenrich"='ce', "Broadenrich"='be',  "Seq2pathway"='sy')),
      br(),
      fileInput("bds", "Upload benchmark dataset",
              accept = c(
                "text/csv",
                "tab/comma-separated-values,text/plain",
                ".csv",
                ".txt")),
      tags$hr(),
      checkboxInput("header", "Header", TRUE)),
    
    # Show file contents.
    mainPanel(
      tabsetPanel(
        tabPanel("Preview", tableOutput("contents")))
#        tabPanel("Plot"))
      )
    ))
      

# Server logic 
server <- function(input, output) {

  ## Display the contents of the table in the main panel.  
  output$contents <- renderTable({
    inFile <- input$bds
    if (is.null(inFile))
      return(NULL)
      head(read.table(inFile$datapath, header = input$header, sep = "\t", quote = ""))
  })
  
  
  ## Visualization plot of the data.
#  source_url("https://github.com/mora-lab/benchmarks/blob/master/genomic_range/R/Plotting_Comparison_Metrics.R")
}

# Complete app with UI and server components
shinyApp(ui, server)
}
