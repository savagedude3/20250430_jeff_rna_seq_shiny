#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(shiny)


#<< ui
#define choices for input
directions <- c("downregulated", "upregulated", "both")
#plot_list <- c("dotplot", "heatplot", "upsetplot")
ont_list <- c("MF", "CC", "BP")

data_in <- read_csv("20250430_deseq2_entrez.csv")

ui <- fluidPage(
  fluidRow(
    column(6,
           numericInput("p_value", "p_value", value = 0.05),
           numericInput("log2fc", "log2fc", value = 1),
           selectInput("direction_of_change", "direction_of_change", choices = directions),
           #selectInput("plot_type", "plot_type", choices = plot_list),
           selectInput("ont_type", "ont_type", choices = ont_list)
           #actionButton("plot", "Generate Plot")
    )
  ),
  
  # fluidRow(
  #   column(12, tableOutput("sig_genes"))
  # ),
  
  downloadButton("sig_genes_download", "Download sig_genes"
  ),
  
  # fluidRow(
  #   column(12, tableOutput("go_table"))
  # ),
  
  downloadButton("go_ora_download", "Download go_ora"
  ),
  
  downloadButton("go_dotplot_download", "Download go_dotplot"
  ),
  
  fluidRow(
    column(12, plotOutput("go_dotplot"))
  )
  
)
#>>

#<< server
server <- function(input, output, session) {
  
  
  filtered_data <- reactive({
    data_in <- read_csv("20250430_deseq2_entrez.csv")
    
    #print(head(data_in))
    
    if(input$direction_of_change == "upregulated"){
      current_data <- data_in %>% 
        dplyr::filter(data_in, padj < input$p_value, log2FoldChange > input$log2fc) 
    }
    
    if(input$direction_of_change == "downregulated"){
      current_data <- data_in %>% 
        dplyr::filter(padj < input$p_value, log2FoldChange < (input$log2fc * -1)) 
    }
    
    if(input$direction_of_change == "both"){
      current_data <- data_in %>% 
        dplyr::filter(padj < input$p_value, abs(log2FoldChange) > input$log2fc) 
    }
    
    current_data
    
  })
  
  output$sig_genes <- renderTable(width = 800, height = 800,
    {
    
    filtered_data()

    }
  )
  
  output$sig_genes_download <- downloadHandler(
    filename = function() {
      paste("sig_genes-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(filtered_data(), file)
    }
  )
  
  go_ora <- reactive(
                               {
                                 data_in <- read_csv("20250430_deseq2_entrez.csv")
                                 
                                 sig_genes <- filtered_data() %>% pull(entrezgene_id)
                                 
                                 enrichGO(gene = as.character(sig_genes),
                                                    OrgDb = org.Mm.eg.db,
                                                    universe = as.character(data_in$entrezgene_id),
                                                    ont = input$ont_type,
                                                    readable = TRUE) # maps gene IDs to gene names
                                 
  
                               })
  
  output$go_ora_download <- downloadHandler(
    filename = function() {
      paste("go_ora-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(go_ora(), file)
    }
  )
  
  output$go_table <- renderTable(width = 800, height = 800,{
    go_ora()
  })
  
  plot1 <- reactive({
    go_ora() %>% 
      dotplot(showCategory = 15) +
      ggtitle("dotplot for ORA") 
  })
  
  output$go_dotplot <- renderPlot(width = 800, height = 800,
                           {
                            plot1()
                             })

  
  output$go_dotplot_download <- downloadHandler(
    filename = function() {
      paste("go_dotplot-", Sys.Date(), ".png", sep="")
    },
    content = function(file) {
      plot1()
      ggsave(file, width = 800, height = 800)
    }
  )
  
}
#>>

shinyApp(ui, server)
