if(interactive()){

  ui <- fluidPage(
  titlePanel(title = "Welcome to the Benchmarking of Gene Set Analysis tools."),
  br(),
  h5('This application allows analyzing a user defined dataset for ranking various enrichment tools for genomic regions.
   The current suite of tools include GREAT, Enrichr, Chipenrich, Broadenrich, Seq2pathway. However, only the latter 3 are available
   as functions in R and hence form a part of this interface. The user is required to select the dataset, tool, and the metric of
   comparison. The application returns the plot and values as a result.'),
  h5('This interactive module allows for reviewing the results from our current study of benchmark data as we formulated.'),  
  br(),
  sidebarLayout(
    
    sidebarPanel(
      helpText('These four metrics are necessary and sufficient to highlight efficacy of a tool in classification task. However, as highlighted 
               by Tarca et al. (https://doi.org/10.1371/journal.pone.0079217) we dodge the case of surrogate sensitivity, not finding it
               befitting in our scenario.'),
      selectInput(inputId = "m", label = "Select a comparison metric for viewing plots from our benchmark dataset",
                   choices = c("Sensitivity"='sn', "Specificity"='sp',  "Prioritization"='pn', "Precision"='pr'),
                  selected = 'sn'),
      br(),
      radioButtons(inputId = "t", label = "Select a tool",
                   choices = c("GREAT"='gt', "Chipenrich"='ce', "Broadenrich"='be',  "Seq2pathway"='sy', "Enrichr"='er'),
                   selected = 'gt'),
      br(),
      radioButtons(inputId = "d", label = "Select a disease",
                   choices = c("Colorectal Cancer"='cc', "Prostate Cancer"='pc',  "Gastric Cancer"='gc', 
                               "Alzheimer's Disease"='ad'),
                   selected = 'cc'),
      br(),
      helpText('Please upload the metadata for your benchmark dataset here. The file may be structured with the 
               following variables (for reference): S.No., 	Experimental Method, 	Organism , Cell Type/Tissue, 	
                Genome Version,Disease Target Pathway, Samples,	Replicates,	PMID, Publication Journal, 	
                Publication Year, GSE, GSM, TF or histone mark, Cistrome ID, KEGG ID,	Metadata.'),
      br(),
      helpText('However, the mandatory attributes include, Disease Target Pathway and GSM, as the samples will be
               loaded and processed against the disease terms as available with the KEGG and GO ids (BP,CC,MF).'),
      fileInput(inputId = "bds", label = "Upload your benchmark dataset",
              accept = c(
                "text/csv",
                "tab/comma-separated-values,text/plain",
                ".csv",
                ".txt")),
      tags$hr(), # Horizontal rule
      checkboxInput(inputId = "header", label = "Header", TRUE),
      actionButton(inputId = "submit", label="Submit", icon("fas fa-magic"))
      ),
    
    
    mainPanel(
      tabsetPanel( type = "tab",
                   tabPanel("Dataset Reference",
                            tags$iframe(style="height:1000px; width:100%; scrolling=yes", 
                                        src="GSA_ChIP_Seq_Master_Table.pdf")), # Summary
                   tabPanel("Preview", tableOutput(outputId = "contents")), # Show file contents.
                   tabPanel("Benchmark Plots",
                            conditionalPanel(condition = "input.m=='sn'", tags$img(src="Sensitivity.jpeg", 
                                                                                   height="700", 
                                                                                   width="700",
                                                                                   align="center")),
                            conditionalPanel(condition = "input.m=='sp'", tags$img(src="Specificity.jpeg",
                                                                                   height="700", 
                                                                                   width="700",
                                                                                   align="center")),
                            conditionalPanel(condition = "input.m=='pr'", tags$img(src="Prioritization.jpeg",
                                                                                   height="700", 
                                                                                   width="700",
                                                                                   align="center")),
                            conditionalPanel(condition = "input.m=='pn'", tags$img(src="Precision.jpeg",
                                                                                   height="700", 
                                                                                   width="700",
                                                                                   align="center"))),     # Display plot
                   tabPanel("Plots", plotOutput(outputId = "userplot")) # User selected plots         
      
    )
  
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
  
  ## Aggregating plots from the radio button inputs from user- tool and disease.
    for_plot <- function(tool,disease){

      ## Initialize local variables for input from user
      tool <- input$t
      disease <- input$d


      ## Import results for each metric from external, manually curated files.

      precision_data <- read.table("./Benchmark_Demo/www/Precision_table.txt", sep = "\t", header = TRUE, quote = "")
      prioritization_data <- read.table("./Benchmark_Demo/www/Prioritization_table.txt", sep = "\t", header = TRUE, quote = "")
      sensitivity_data <- read.table("./Benchmark_Demo/www/Sensitivity_table.txt", sep = "\t", header = TRUE, quote = "")
      specificity_data <- read.table("./Benchmark_Demo/www/Specificity_table.txt", sep = "\t", header = TRUE, quote = "")
      
      ## Creating all possible combinations for the tool and disease inputs from the user, 4 diseases * 5 tools, i.e. 20 precisely. 
      ## Colorectal Cancer and Chipenrich
      if(input$d <- 'cc' && input$t <- 'ce')
      {
            db_sn <- sensitivity_data[,which(grepl("Colorectal", colnames(sensitivity_data[,which(grepl("Chipenrich", colnames(sensitivity_data)))])))]
            db_sp <- specificity_data[,which(grepl("Colorectal", colnames(specificity_data[,which(grepl("Chipenrich", colnames(specificity_data)))])))]
            db_pr <- prioritization_data[,which(grepl("Colorectal", colnames(prioritization_data[,which(grepl("Chipenrich", colnames(prioritization_data)))])))]
            db_pn <- precision_data[,which(grepl("Colorectal", colnames(precision_data[,which(grepl("Chipenrich", colnames(precision_data)))])))]
            db <- cbind(db_sn,db_sp,db_pr,db_pn)
            boxplot(db,
                  boxwex = 0.1,
                  names = c("Sensitivity","Specificity","Prioritization","Precision"),
                  ylab = "Colorectal Cancer",
                  xlab = "Chipenrich",
                  ylim = c(0,100),
                  col = c("salmon","tan","khaki","lavender"))
            abline(h=50, col="red")
      }
      
      ## Colorectal Cancer and Broadenrich
      if(input$d <- 'cc' && input$t <- 'be')
      {
        db_sn <- sensitivity_data[,which(grepl("Colorectal", colnames(sensitivity_data[,which(grepl("Broadenrich", colnames(sensitivity_data)))])))]
        db_sp <- specificity_data[,which(grepl("Colorectal", colnames(specificity_data[,which(grepl("Broadenrich", colnames(specificity_data)))])))]
        db_pr <- prioritization_data[,which(grepl("Colorectal", colnames(prioritization_data[,which(grepl("Broadenrich", colnames(prioritization_data)))])))]
        db_pn <- precision_data[,which(grepl("Colorectal", colnames(precision_data[,which(grepl("Broadenrich", colnames(precision_data)))])))]
        db <- cbind(db_sn,db_sp,db_pr,db_pn)
        boxplot(db,
                boxwex = 0.1,
                names = c("Sensitivity","Specificity","Prioritization","Precision"),
                ylab = "Colorectal Cancer",
                xlab = "Broadenrich",
                ylim = c(0,100),
                col = c("salmon","tan","khaki","lavender"))
        abline(h=50, col="red")
      }
      
      ## Colorectal Cancer and Seq2pathway
      if(input$d <- 'cc' && input$t <- 'sy')
      {
        db_sn <- sensitivity_data[,which(grepl("Colorectal", colnames(sensitivity_data[,which(grepl("Seq2pathway", colnames(sensitivity_data)))])))]
        db_sp <- specificity_data[,which(grepl("Colorectal", colnames(specificity_data[,which(grepl("Seq2pathway", colnames(specificity_data)))])))]
        db_pr <- prioritization_data[,which(grepl("Colorectal", colnames(prioritization_data[,which(grepl("Seq2pathway", colnames(prioritization_data)))])))]
        db_pn <- precision_data[,which(grepl("Colorectal", colnames(precision_data[,which(grepl("Seq2pathway", colnames(precision_data)))])))]
        db <- cbind(db_sn,db_sp,db_pr,db_pn)
        boxplot(db,
                boxwex = 0.1,
                names = c("Sensitivity","Specificity","Prioritization","Precision"),
                ylab = "Colorectal Cancer",
                xlab = "Seq2pathway",
                ylim = c(0,100),
                col = c("salmon","tan","khaki","lavender"))
        abline(h=50, col="red")
      }
      
      ## Colorectal Cancer and Enrichr
      if(input$d <- 'cc' && input$t <- 'er')
      {
        db_sn <- sensitivity_data[,which(grepl("Colorectal", colnames(sensitivity_data[,which(grepl("Enrichr", colnames(sensitivity_data)))])))]
        db_sp <- specificity_data[,which(grepl("Colorectal", colnames(specificity_data[,which(grepl("Enrichr", colnames(specificity_data)))])))]
        db_pr <- prioritization_data[,which(grepl("Colorectal", colnames(prioritization_data[,which(grepl("Enrichr", colnames(prioritization_data)))])))]
        db_pn <- precision_data[,which(grepl("Colorectal", colnames(precision_data[,which(grepl("Enrichr", colnames(precision_data)))])))]
        db <- cbind(db_sn,db_sp,db_pr,db_pn)
        boxplot(db,
                boxwex = 0.1,
                names = c("Sensitivity","Specificity","Prioritization","Precision"),
                ylab = "Colorectal Cancer",
                xlab = "Enrichr",
                ylim = c(0,100),
                col = c("salmon","tan","khaki","lavender"))
        abline(h=50, col="red")
      }
      
      ## Colorectal Cancer and GREAT
      if(input$d <- 'cc' && input$t <- 'gt')
      {
        db_sn <- sensitivity_data[,which(grepl("Colorectal", colnames(sensitivity_data[,which(grepl("GREAT", colnames(sensitivity_data)))])))]
        db_sp <- specificity_data[,which(grepl("Colorectal", colnames(specificity_data[,which(grepl("GREAT", colnames(specificity_data)))])))]
        db_pr <- prioritization_data[,which(grepl("Colorectal", colnames(prioritization_data[,which(grepl("GREAT", colnames(prioritization_data)))])))]
        db_pn <- precision_data[,which(grepl("Colorectal", colnames(precision_data[,which(grepl("GREAT", colnames(precision_data)))])))]
        db <- cbind(db_sn,db_sp,db_pr,db_pn)
        boxplot(db,
                boxwex = 0.1,
                names = c("Sensitivity","Specificity","Prioritization","Precision"),
                ylab = "Colorectal Cancer",
                xlab = "GREAT",
                ylim = c(0,100),
                col = c("salmon","tan","khaki","lavender"))
        abline(h=50, col="red")
      }
  
    }
    
  ## Visualization plot of the data.
#  source_url("https://github.com/mora-lab/benchmarks/blob/master/genomic_range/R/Plotting_Comparison_Metrics.R")
}

# Complete app with UI and server components
shinyApp(ui = ui, server = server)
}
