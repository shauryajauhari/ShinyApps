## Loading libraries.
library(shiny)

## Main module.

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
      br(),
      selectInput(inputId = "m", label = "Select a comparison metric for viewing plots from our benchmark dataset",
                  choices = c("Sensitivity"='sn', "Specificity"='sp',  "Prioritization"='pn', "Precision"='pr'),
                  selected = 'sn'),
      br(),
      radioButtons(inputId = "t", label = "Select a tool",
                   choices = c("Chipenrich"='ce', "Broadenrich"='be',  "Seq2pathway"='sy', "Enrichr"='er', "GREAT"='gt'),
                   selected = 'ce'),
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
               loaded and processed against the disease terms as available with the KEGG and GO ids (BP, CC, MF).'),
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
                   tabPanel("Dataset Plots",
                            conditionalPanel(condition = "input.m=='sn'", tags$img(src="Sensitivity_ggplot.jpeg", 
                                                                                   height="500", 
                                                                                   width="1000",
                                                                                   align="center")),
                            conditionalPanel(condition = "input.m=='sp'", tags$img(src="Specificity_ggplot.jpeg",
                                                                                   height="500", 
                                                                                   width="1000",
                                                                                   align="center")),
                            conditionalPanel(condition = "input.m=='pr'", tags$img(src="Prioritization_ggplot.jpeg",
                                                                                   height="500", 
                                                                                   width="1000",
                                                                                   align="center")),
                            conditionalPanel(condition = "input.m=='pn'", tags$img(src="Precision_ggplot.jpeg",
                                                                                   height="500", 
                                                                                   width="1000",
                                                                                   align="center"))),     # Display plot
                   tabPanel("User Plots", plotOutput(outputId = "userplot", width = "100%")), # User selected plots         
                   tabPanel("ROC Plot",
                            tags$img(src="ROC_Plot.jpeg",
                                     height="500",
                                     width="1000",
                                     align="center")) # ROC Plot
                   
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
  for_plot <- function(){
    
    
    if (!is.null(input$t) && !is.null(input$t)) ## Check for valid inputs
    {
      
      ## Import results for each metric from external, manually curated files.
      
      precision_data <- read.table("www/Precision_Table.txt", sep = "\t", header = TRUE, quote = "")
      prioritization_data <- read.table("www/Prioritization_Table.txt", sep = "\t", header = TRUE, quote = "")
      sensitivity_data <- read.table("www/Sensitivity_Table.txt", sep = "\t", header = TRUE, quote = "")
      specificity_data <- read.table("www/Specificity_Table.txt", sep = "\t", header = TRUE, quote = "")
      
      ## Creating all possible combinations for the tool and disease inputs from the user, 4 diseases * 5 tools, i.e. 20 precisely. 
      ## Colorectal Cancer and Chipenrich
      
      if(input$d == 'cc' && input$t == 'ce')
      {
        db_sn <- sensitivity_data[,which(grepl("Chipenrich_Colorectal", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Chipenrich_Colorectal", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Chipenrich_Colorectal", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Chipenrich_Colorectal", colnames(precision_data)))]
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
      if(input$d == 'cc' && input$t == 'be')
      {
        db_sn <- sensitivity_data[,which(grepl("Broadenrich_Colorectal", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Broadenrich_Colorectal", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Broadenrich_Colorectal", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Broadenrich_Colorectal", colnames(precision_data)))]
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
      if(input$d == 'cc' && input$t == 'sy')
      {
        db_sn <- sensitivity_data[,which(grepl("Seq2pathway_Colorectal", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Seq2pathway_Colorectal", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Seq2pathway_Colorectal", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Seq2pathway_Colorectal", colnames(precision_data)))]
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
      if(input$d == 'cc' && input$t == 'er')
      {
        db_sn <- sensitivity_data[,which(grepl("Enrichr_Colorectal", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Enrichr_Colorectal", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Enrichr_Colorectal", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Enrichr_Colorectal", colnames(precision_data)))]
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
      if(input$d == 'cc' && input$t == 'gt')
      {
        db_sn <- sensitivity_data[,which(grepl("GREAT_Colorectal", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("GREAT_Colorectal", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("GREAT_Colorectal", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("GREAT_Colorectal", colnames(precision_data)))]
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
      
      ## Gastric Cancer and Chipenrich
      if(input$d == 'gc' && input$t == 'ce')
      {
        db_sn <- sensitivity_data[,which(grepl("Chipenrich_Gastric", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Chipenrich_Gastric", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Chipenrich_Gastric", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Chipenrich_Gastric", colnames(precision_data)))]
        db <- cbind(db_sn,db_sp,db_pr,db_pn)
        boxplot(db,
                boxwex = 0.1,
                names = c("Sensitivity","Specificity","Prioritization","Precision"),
                ylab = "Gastric Cancer",
                xlab = "Chipenrich",
                ylim = c(0,100),
                col = c("salmon","tan","khaki","lavender"))
        abline(h=50, col="red")
      }
      
      ## Gastric Cancer and Broadenrich
      if(input$d == 'gc' && input$t == 'be')
      {
        db_sn <- sensitivity_data[,which(grepl("Broadenrich_Gastric", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Broadenrich_Gastric", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Broadenrich_Gastric", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Broadenrich_Gastric", colnames(precision_data)))]
        db <- cbind(db_sn,db_sp,db_pr,db_pn)
        boxplot(db,
                boxwex = 0.1,
                names = c("Sensitivity","Specificity","Prioritization","Precision"),
                ylab = "Gastric Cancer",
                xlab = "Broadenrich",
                ylim = c(0,100),
                col = c("salmon","tan","khaki","lavender"))
        abline(h=50, col="red")
      }
      
      ## Gastric Cancer and Seq2pathway
      if(input$d == 'gc' && input$t == 'sy')
      {
        db_sn <- sensitivity_data[,which(grepl("Seq2pathway_Gastric", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Seq2pathway_Gastric", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Seq2pathway_Gastric", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Seq2pathway_Gastric", colnames(precision_data)))]
        db <- cbind(db_sn,db_sp,db_pr,db_pn)
        boxplot(db,
                boxwex = 0.1,
                names = c("Sensitivity","Specificity","Prioritization","Precision"),
                ylab = "Gastric Cancer",
                xlab = "Seq2pathway",
                ylim = c(0,100),
                col = c("salmon","tan","khaki","lavender"))
        abline(h=50, col="red")
      }
      
      ## Gastric Cancer and Enrichr
      if(input$d == 'gc' && input$t == 'er')
      {
        db_sn <- sensitivity_data[,which(grepl("Enrichr_Gastric", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Enrichr_Gastric", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Enrichr_Gastric", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Enrichr_Gastric", colnames(precision_data)))]
        db <- cbind(db_sn,db_sp,db_pr,db_pn)
        boxplot(db,
                boxwex = 0.1,
                names = c("Sensitivity","Specificity","Prioritization","Precision"),
                ylab = "Gastric Cancer",
                xlab = "Enrichr",
                ylim = c(0,100),
                col = c("salmon","tan","khaki","lavender"))
        abline(h=50, col="red")
      }
      
      ## Gastric Cancer and GREAT
      if(input$d == 'gc' && input$t == 'gt')
      {
        db_sn <- sensitivity_data[,which(grepl("GREAT_Gastric", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("GREAT_Gastric", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("GREAT_Gastric", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("GREAT_Gastric", colnames(precision_data)))]
        db <- cbind(db_sn,db_sp,db_pr,db_pn)
        boxplot(db,
                boxwex = 0.1,
                names = c("Sensitivity","Specificity","Prioritization","Precision"),
                ylab = "Gastric Cancer",
                xlab = "GREAT",
                ylim = c(0,100),
                col = c("salmon","tan","khaki","lavender"))
        abline(h=50, col="red")
      }
      
      ## Prostate Cancer and Chipenrich
      if(input$d == 'pc' && input$t == 'ce')
      {
        db_sn <- sensitivity_data[,which(grepl("Chipenrich_Prostate", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Chipenrich_Prostate", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Chipenrich_Prostate", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Chipenrich_Prostate", colnames(precision_data)))]
        db <- cbind(db_sn,db_sp,db_pr,db_pn)
        boxplot(db,
                boxwex = 0.1,
                names = c("Sensitivity","Specificity","Prioritization","Precision"),
                ylab = "Prostate Cancer",
                xlab = "Chipenrich",
                ylim = c(0,100),
                col = c("salmon","tan","khaki","lavender"))
        abline(h=50, col="red")
      }
      
      ## Prostate Cancer and Broadenrich
      if(input$d == 'pc' && input$t == 'be')
      {
        db_sn <- sensitivity_data[,which(grepl("Broadenrich_Prostate", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Broadenrich_Prostate", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Broadenrich_Prostate", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Broadenrich_Prostate", colnames(precision_data)))]
        db <- cbind(db_sn,db_sp,db_pr,db_pn)
        boxplot(db,
                boxwex = 0.1,
                names = c("Sensitivity","Specificity","Prioritization","Precision"),
                ylab = "Prostate Cancer",
                xlab = "Broadenrich",
                ylim = c(0,100),
                col = c("salmon","tan","khaki","lavender"))
        abline(h=50, col="red")
      }
      
      ## Prostate Cancer and Seq2pathway
      if(input$d == 'pc' && input$t == 'sy')
      {
        db_sn <- sensitivity_data[,which(grepl("Seq2pathway_Prostate", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Seq2pathway_Prostate", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Seq2pathway_Prostate", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Seq2pathway_Prostate", colnames(precision_data)))]
        db <- cbind(db_sn,db_sp,db_pr,db_pn)
        boxplot(db,
                boxwex = 0.1,
                names = c("Sensitivity","Specificity","Prioritization","Precision"),
                ylab = "Prostate Cancer",
                xlab = "Seq2pathway",
                ylim = c(0,100),
                col = c("salmon","tan","khaki","lavender"))
        abline(h=50, col="red")
      }
      
      ## Prostate Cancer and Enrichr
      if(input$d == 'pc' && input$t == 'er')
      {
        db_sn <- sensitivity_data[,which(grepl("Enrichr_Prostate", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Enrichr_Prostate", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Enrichr_Prostate", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Enrichr_Prostate", colnames(precision_data)))]
        db <- cbind(db_sn,db_sp,db_pr,db_pn)
        boxplot(db,
                boxwex = 0.1,
                names = c("Sensitivity","Specificity","Prioritization","Precision"),
                ylab = "Prostate Cancer",
                xlab = "Enrichr",
                ylim = c(0,100),
                col = c("salmon","tan","khaki","lavender"))
        abline(h=50, col="red")
      }
      
      ## Prostate Cancer and GREAT
      if(input$d == 'pc' && input$t == 'gt')
      {
        db_sn <- sensitivity_data[,which(grepl("GREAT_Prostate", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("GREAT_Prostate", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("GREAT_Prostate", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("GREAT_Prostate", colnames(precision_data)))]
        db <- cbind(db_sn,db_sp,db_pr,db_pn)
        boxplot(db,
                boxwex = 0.1,
                names = c("Sensitivity","Specificity","Prioritization","Precision"),
                ylab = "Prostate Cancer",
                xlab = "GREAT",
                ylim = c(0,100),
                col = c("salmon","tan","khaki","lavender"))
        abline(h=50, col="red")
      }
      
      ## Alzheimer's Disease and Chipenrich
      if(input$d == 'ad' && input$t == 'ce')
      {
        db_sn <- sensitivity_data[,which(grepl("Chipenrich_Alzheimer", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Chipenrich_Alzheimer", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Chipenrich_Alzheimer", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Chipenrich_Alzheimer", colnames(precision_data)))]
        db <- cbind(db_sn,db_sp,db_pr,db_pn)
        boxplot(db,
                boxwex = 0.1,
                names = c("Sensitivity","Specificity","Prioritization","Precision"),
                ylab = "Alzheimer's Disease",
                xlab = "Chipenrich",
                ylim = c(0,100),
                col = c("salmon","tan","khaki","lavender"))
        abline(h=50, col="red")
      }
      
      ## Alzheimer's Disease and Broadenrich
      if(input$d == 'ad' && input$t == 'be')
      {
        db_sn <- sensitivity_data[,which(grepl("Broadenrich_Alzheimer", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Broadenrich_Alzheimer", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Broadenrich_Alzheimer", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Broadenrich_Alzheimer", colnames(precision_data)))]
        db <- cbind(db_sn,db_sp,db_pr,db_pn)
        boxplot(db,
                boxwex = 0.1,
                names = c("Sensitivity","Specificity","Prioritization","Precision"),
                ylab = "Alzheimer's Disease",
                xlab = "Broadenrich",
                ylim = c(0,100),
                col = c("salmon","tan","khaki","lavender"))
        abline(h=50, col="red")
      }
      
      ## Alzheimer's Disease and Seq2pathway
      if(input$d == 'ad' && input$t == 'sy')
      {
        db_sn <- sensitivity_data[,which(grepl("Seq2pathway_Alzheimer", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Seq2pathway_Alzheimer", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Seq2pathway_Alzheimer", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Seq2pathway_Alzheimer", colnames(precision_data)))]
        db <- cbind(db_sn,db_sp,db_pr,db_pn)
        boxplot(db,
                boxwex = 0.1,
                names = c("Sensitivity","Specificity","Prioritization","Precision"),
                ylab = "Alzheimer's Disease",
                xlab = "Seq2pathway",
                ylim = c(0,100),
                col = c("salmon","tan","khaki","lavender"))
        abline(h=50, col="red")
      }
      
      ## Alzheimer's Disease and Enrichr
      if(input$d == 'ad' && input$t == 'er')
      {
        db_sn <- sensitivity_data[,which(grepl("Enrichr_Alzheimer", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("Enrichr_Alzheimer", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("Enrichr_Alzheimer", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("Enrichr_Alzheimer", colnames(precision_data)))]
        db <- cbind(db_sn,db_sp,db_pr,db_pn)
        boxplot(db,
                boxwex = 0.1,
                names = c("Sensitivity","Specificity","Prioritization","Precision"),
                ylab = "Alzheimer's Disease",
                xlab = "Enrichr",
                ylim = c(0,100),
                col = c("salmon","tan","khaki","lavender"))
        abline(h=50, col="red")
      }
      
      ## Alzheimer's Disease and GREAT
      if(input$d == 'ad' && input$t == 'gt')
      {
        db_sn <- sensitivity_data[,which(grepl("GREAT_Alzheimer", colnames(sensitivity_data)))]
        db_sp <- specificity_data[,which(grepl("GREAT_Alzheimer", colnames(specificity_data)))]
        db_pr <- prioritization_data[,which(grepl("GREAT_Alzheimer", colnames(prioritization_data)))]
        db_pn <- precision_data[,which(grepl("GREAT_Alzheimer", colnames(precision_data)))]
        db <- cbind(db_sn,db_sp,db_pr,db_pn)
        boxplot(db,
                boxwex = 0.1,
                names = c("Sensitivity","Specificity","Prioritization","Precision"),
                ylab = "Alzheimer's Disease",
                xlab = "GREAT",
                ylim = c(0,100),
                col = c("salmon","tan","khaki","lavender"))
        abline(h=50, col="red")
      }
    }
  }
  
  
  ## Deploying function for plotting at the click of the action button.
  
  output$userplot <- renderPlot(
    {
      # Take a dependency on 'input$submit'. This will run once initially, because the value changes from NULL to 0.
      input$submit
      for_plot()
    }, height = 500, width = 500)
  
}

# Complete app with UI and server components
shinyApp(ui = ui, server = server)