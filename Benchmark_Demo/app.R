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
      helpText("These four metrics are necessary and sufficient to highlight efficacy of a tool in classification task. However, as highlighted 
               by Tarca et al. (https://doi.org/10.1371/journal.pone.0079217) we dodge the case of surrogate sensitivity, not finding it
               befitting in our scenario."),
      selectInput(inputId = "m", label = "Select a comparison metric",
                   choices = c("Sensitivity"='sn', "Specificity"='sp',  "Prioritization"='pn', "Precision"='pr'),
                  selected = 'sn'),
      br(),
      radioButtons(inputId = "t", label = "Select a tool",
                   choices = c("Chipenrich"='ce', "Broadenrich"='be',  "Seq2pathway"='sy'),
                   selected = 'ce'),
      br(),
      fileInput(inputId = "bds", label = "Upload benchmark dataset",
              accept = c(
                "text/csv",
                "tab/comma-separated-values,text/plain",
                ".csv",
                ".txt")),
      tags$hr(), # Horizontal rule
      checkboxInput(inputId = "header", label = "Header", TRUE)),
    
    
    mainPanel(
      tabsetPanel( type = "tab",
        tabPanel("Preview", tableOutput(outputId = "contents")), # Show file contents.
        tabPanel("Summary", verbatimTextOutput(outputId = "summary")), # Summary
        tabPanel("Plots", 
                 conditionalPanel(condition = "input.m=='sn'", tags$img(src="Sensitivity.jpeg")),
                 conditionalPanel(condition = "input.m=='sp'", tags$img(src="Specificity.jpeg")),
                 conditionalPanel(condition = "input.m=='pr'", tags$img(src="Prioritization.jpeg")),
                 conditionalPanel(condition = "input.m=='pn'", tags$img(src="Precision.jpeg"))     # Display plot
      )
    ),
  actionButton(inputId = "submit", label="Submit", icon("fas fa-magic"))
  )))

# Server logic 
server <- function(input, output, session) {

  ## Display the contents of the table in the main panel.  
    output$contents <- renderTable({
      if (is.null(input$bds))
        return(NULL)
        head(read.table(input$bds$datapath, header = input$header, sep = "\t", quote = ""))
    })
    
  ## Display the summary of the table.
    output$summary <- renderPrint({
      summary(input$bds)
    })
  
  ## Visualization plot of the data.
#  source_url("https://github.com/mora-lab/benchmarks/blob/master/genomic_range/R/Plotting_Comparison_Metrics.R")
}

# Complete app with UI and server components
shinyApp(ui = ui, server = server)
}
