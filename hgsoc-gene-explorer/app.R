library(shiny)
library(shinycssloaders)
library(dplyr)

# Specify the application port
options(shiny.host = "0.0.0.0")
options(shiny.port = 3838)

# --- Load helper functions ---
source("R/load_data.R")
source("R/plot_functions.R")
source("R/utils.R")

# --- Load data once at startup ---
datasets <- load_all_data()

ui <- fluidPage(
  
  titlePanel("HGSOC NanoString & RNAseq\nGene Explorer"),

  sidebarLayout(
    sidebarPanel(
      p(HTML("<strong>Click Run</strong> to check out how your favorite gene looks in a large meta-cohort of 
      HGSOC patients with paired samples pre and post NACT.")),
      p(HTML("If the gene is found in the NanoString panel, you can also see how well the gene aligns
      in expression according to samples sequenced with both NanoString and RNA-seq")), 
      br(),
      textInput("gene", "Enter Gene Name:", value = "RRM2"),
      actionButton("go", "Run"),
      br(),
      br(),
      p("Please cite the following:"),
      width = 3
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Pre & Post Expression\nPFS Binary",
                helpText("Normalized expression (RPKM for RNA-seq) before and after NACT treatment of paired-patients
                split by Progression-Free survival being >12 months or <=12 months. A linear regression is performed with the model and R^2 value.
                If NanoString Experiments are not included, it is because the gene is not in the NanoString panel."), 
                br(),
                 plotOutput("exprPlot") %>% withSpinner()   
        ),
        tabPanel("Expression & PFS\nCorrelation",
                helpText("Normalized expression (RPKM for RNA-seq) before and after NACT treatment and their
                correlation to PFS (in months). A linear regression is performed with the model and R^2 value. 
                The horizontal line pinpoints the 12 month cutoff. "), 
                br(),
                 plotOutput("exprCorrPlot") %>% withSpinner()
                
        ), 
        tabPanel("Log2FC\nPFS Binary",
                 helpText("Log2 Fold Change between Post and Pre NACT treatment across patients
                 split by Progression-Free survival being >12 months or <=12 months. 
                   T-tests are performed for each of the experiments and if combining all experiments, with p-values shown"), 
                br(),
                 plotOutput("l2fcPlot") %>% withSpinner()   
        ),
        tabPanel("Log2FC & PFS\nCorrelation",
                helpText("Log2 Fold Change between Post and Pre NACT treatment and their correlation to PFS (in months). 
                   A linear regression is performed with the model and R^2 value.  
                   The horizontal line pinpoints the 12 month cutoff. "), 
                br(),
                plotOutput("l2fcCorrPlot") %>% withSpinner()
                
        ), 
        tabPanel("NanoString/RNAseq Correlation",
          helpText(HTML("Scatterplots of normalized counts from NanoString (y-axis) and from RNA-seq using three 
                pre-processing approaches for RNA-seq. All RNA-seq was normalized using RPKM (TPM and MR showed no clear differences). 
                Data from FFPE tumor samples sequenced with both RNA-seq and NanoString.
                    <ul>
                      <li>Isoform_RPKM: Used counts over isoforms best matching the NanoString probe.</li>
                      <li>RNAIsoform_RPKM: Used counts over isoforms most highly expressed in the RNA-seq.</li>
                      <li>Exon_RPKM: Using counts over the exon(s) aligning with the NanoString probes.</li>
                    </ul>")), 
                br(),
                 plotOutput("Matched") %>% withSpinner()
                
        )
      )
    )
  )
)

server <- function(input, output, session) {

  ###################
  ## SHARED DATA ####
  selected_gene <- eventReactive(input$go, {
    req(input$gene)
    clean_gene <- sanitize_gene(input$gene)
    return(clean_gene)
  })
  # determine if nano gene or not ONLY when gene changes
  nano_gene <- eventReactive(selected_gene(),{
    gene <- selected_gene()
    gene %in% rownames(datasets$exp_data$Bitler)
  })

  # get the plot df ONLY when gene changes
  plot_df <- eventReactive(selected_gene(), {
    get_matched_plot_df(selected_gene(), datasets$exp_data, datasets$meta, nano_gene())
  })

  ################
  ## PLOTTING ####
  output$exprPlot <- renderPlot({
    # plot the expression with PFS binary
    plot_matched_gene(selected_gene(), plot_df())
  },
  height = 700)

  output$exprCorrPlot <- renderPlot({
    # plot the correlation with expression
    plot_exp_corr_PFS(selected_gene(), plot_df())
  },
  height = 800)

  output$l2fcPlot <- renderPlot({
    # plot the L2FC split by PFS binary
    plot_l2fc(selected_gene(), datasets$l2fc_data, nano_gene())
  },
  width = 800)

  output$l2fcCorrPlot <- renderPlot({
    # plot hte correlation with Log2FC
    plot_l2fc_corr_PFS(selected_gene(), datasets$l2fc_data, nano_gene())
  }, 
  height=1000)

  output$Matched <- renderPlot({
    # plot hte correlation with Log2FC
    plot_l2fc_corr_PFS(selected_gene(), datasets$l2fc_data, nano_gene())
  }, 
  height=1000)
  
  # plot the correlation between NanoString and RNA-seq
  output$Matched <- renderPlot({
    is_nano = nano_gene()
    gene = selected_gene()
    if (is_nano) {
      plot_gene_scatter_main(gene, datasets$matched_data, datasets$matched_data$Nano)
    } else {
      plot.new()
      text(0.5, 0.5, paste("Gene", gene, "not included in tested NanoString panel"))
      return()

    }
    
  }, height=400)
}


shinyApp(ui, server)