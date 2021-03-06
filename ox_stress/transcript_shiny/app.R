# ------------------------------------------------------------------------------
# Emily's Embryo Transcriptomic data -- Shiny App
# TS O'Leary
# ------------------------------------------------------------------------------

# Use pacman package to install & load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(shiny, 
               DESeq2, 
               tidyverse, 
               kableExtra)

# Load data
ddslrtreg <- readRDS("ddslrtreg.rds")
region_res <- read_csv("tropvtempREGION.csv") # %>%
#   filter(!is.na(padj))
genes_df <- read_csv("emily_transcript_genes.csv") %>%
  filter(transcript_id %in% region_res$transcript_id)
res <- full_join(region_res, genes_df)
lt_df <- read_delim("lockwood_2018_LT50_data.txt", delim = "\t")

# Define UI for app that draws boxplots ----------------------------------------
ui <- fluidPage(
  
  # App title ----
  titlePanel("Emily's transcriptomics data: heat shock of 0-1 hr embryos"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(width = 3,
                 
                 # Input: Selct genes ----
                 selectInput("gene", 
                             "Pick a gene", 
                             selected = "Ucp4A",
                             sort(unique(genes_df$gene)), 
                             selectize = TRUE),
                 
                 # Input: Group plot by region or population ----
                 radioButtons("factor", "Group by",
                              c("Region" = "region",
                                "Population" = "pop"), 
                              selected = "region"),
                 
                 # Input: Filter significant genes ----
                 checkboxInput("sig_check", 
                               "Plot only significant", 
                               value = TRUE)
                 
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(width = 9,
              
              # Output: Tabset w/ plot and table ----
              tabsetPanel(type = "tabs",
                          tabPanel("Box Plot", 
                                   h2(htmlOutput(outputId = "plot_title")), 
                                   plotOutput(outputId = "boxPlot")),
                          tabPanel("Expression vs. LT50",
                                   h2(htmlOutput(outputId = "plot_title_lt")),
                                   plotOutput(outputId = "boxPlotLT50")),
                          tabPanel("DESeq Region Results", 
                                   h2(htmlOutput(outputId = "table_title")),
                                   tableOutput(outputId = "table")),
                          tabPanel("All DE Transcripts",
                                   h2("All differentially expressed transcripts between regions"),
                                   tableOutput(outputId = "table_sig")))
              
    )
  )
)

# Define server logic required to draw plot ------------------------------------
server <- function(input, output) {
  
  # Create reactive title for plots
  output$plot_title <- reactive({
    paste0("<i>", input$gene, "</i>", " expression")
  })
  
  # Create reactive title for table
  output$table_title <- reactive({
    paste0("DESeq2 results for <i>", input$gene, "</i>", " expression between regions")
  })
  
  # Create reactive title for lt50 expression
  output$plot_title_lt <- reactive({
    paste0("<i>", input$gene, "</i>", " expression vs embryonic LT50")
  })
  
  # Filter data for selected gene and/if significance ----
  data <- reactive({
    y <- res %>% 
      filter(gene == input$gene) %>%
      arrange(desc(baseMean))
    
    if (input$sig_check) {
      y <- y %>%
        filter(padj < 0.05)
    }
    
    validate(
      need(nrow(y) > 0,
           paste0("There are no transcripts for that are significantly different",
                  " between regions for this gene.\nUncheck 'Plot only significant'",
                  " or choose another gene."))
    )
    
    y 
    
  })
  
  
  # Define the plot height dynamically based on the number transcripts ---------
  
  # For Box plot
  plot_height <- reactive({
    y <- data()
    
    nrow(y)*300
  }
  )
  
  # For LT50 exp data
  plot_height_lt <- reactive({
    y <- data()
    
    nrow(y)*600
  }
  )
  
  # Create the box plots based on selected inputs ----
  output$boxPlotLT50 <- renderPlot({
    
    y <- data()
    
    plots <- vector(mode = "list", length = nrow(y))
    
    for (i in 1:nrow(y)){
      
      d <- plotCounts(ddslrtreg, 
                      gene = y$transcript_id[i], 
                      intgroup = c("pop", "temp"), 
                      returnData = TRUE) %>%
        left_join(lt_df, by = c("pop" = "Pop"))
      
      p <- ggplot(d, aes(x = Embryo_LT50, 
                         y = count)) +
        geom_smooth(aes(color = temp), 
                    method = "lm", 
                    se = TRUE) +
        geom_boxplot(aes(fill = pop), 
                     color = "grey50", 
                     alpha = 0.8) +
        scale_fill_manual(name = "Population",
                          values = c("coral", "coral1", 
                                     "coral2", "coral3", "coral4",
                                     "deepskyblue", "deepskyblue1", 
                                     "deepskyblue2", "deepskyblue3", 
                                     "deepskyblue4")) +
        labs(y = "Normalized\ntranscript abundance", 
             x = paste(expression("LT[50]"), "(°C)"),
             title = y$transcript_id[i]) +
        theme_classic(base_size = 14) +
        facet_wrap(~ temp) 
      
      plots[[i]] <- p
    }
    
    cowplot::plot_grid(plotlist = plots, nrow = length(plots))
    
  }, width = 800, height = plot_height_lt)
  
  # Create the box plots of expression vs LT50 for each temperature ----
  output$boxPlot <- renderPlot({
    
    y <- data()
    
    plots <- vector(mode = "list", length = nrow(y))
    
    for (i in 1:nrow(y)){
      
      d <- plotCounts(ddslrtreg, 
                      gene = y$transcript_id[i], 
                      intgroup = c(input$factor, "temp"), 
                      returnData = TRUE) 
      
      p <- ggplot(d, aes(x = temp, 
                         y = count, 
                         fill = d[, input$factor])) + 
        geom_boxplot() +
        labs(y = "Normalized\ntranscript abundance", 
             x = "Heat Shock (°C)",
             title = y$transcript_id[i]) +
        theme_classic(base_size = 14) 
      
      if (input$factor == "region") {
        p <- p + 
          scale_fill_manual(name = "Region",
                            labels = c("Tropical", "Temperate"),
                            values = c("#F8766D", "#00BFC4"))
        
      } else {
        p <- p +
          scale_fill_manual(name = "Population",
                            values = c("coral", "coral1", 
                                       "coral2", "coral3", "coral4",
                                       "deepskyblue", "deepskyblue1", 
                                       "deepskyblue2", "deepskyblue3", 
                                       "deepskyblue4"))
      }
      
      plots[[i]] <- p
    }
    
    cowplot::plot_grid(plotlist = plots, nrow = length(plots))
    
  }, width = 800, height = plot_height)
  
  # Create table with baseMean, log2FC, and padj ----
  output$table <- function() {
    
    res_filt <- res %>% 
      filter(gene == input$gene) %>%
      select(gene, transcript_id, baseMean, log2FoldChange, padj) %>%
      arrange(desc(baseMean))
    
    res_filt  %>%
      kable(format = "html", escape = F) %>%
      kable_styling("striped") %>%
      row_spec(which(res_filt$padj < 0.05), color = "red") 
    
  }
  
  # Create table for all significant genes
  # Create table with baseMean, log2FC, and padj ----
  output$table_sig <- function() {
    
    res %>% 
      filter(padj < 0.05) %>%
      select(gene, transcript_id, baseMean, log2FoldChange, padj) %>%
      arrange(padj) %>%
      kable(format = "html", escape = F) %>%
      kable_styling("striped")
    
  }
  
}

# Create Shiny app -------------------------------------------------------------
shinyApp(ui = ui, server = server)