library(shiny)
library(shinycssloaders)
library(dplyr)

testing=TRUE
# Specify the application port
options(shiny.host = "0.0.0.0")
options(shiny.port = 3838)

#getwd()

# for some reason stupid VSCode can't realize where files are
if (testing) {
  # --- Load helper functions ---
source("hgsoc-gene-explorer/R/load_data.R")
source("hgsoc-gene-explorer/R/plot_functions.R")
source("hgsoc-gene-explorer/R/utils.R")
} else {
  # --- Load helper functions ---
source("R/load_data.R")
source("R/plot_functions.R")
source("R/utils.R")
}



# --- Load data once at startup ---
datasets <- load_all_data(testing)

ui <- fluidPage(
  
  titlePanel("HGSOC NanoString & RNAseq\nGene Explorer"),

  sidebarLayout(
    sidebarPanel(
      p(HTML("<strong>Click Run</strong> to check out how your favorite gene looks in a large meta-cohort of 
      HGSOC patients with paired transcriptomics samples pre and post NACT (more info in Cohort Descriptions).")),
      p(HTML("If the gene is found in the NanoString panel (nCounter PanCancer IO 360 panel), you can also see how well the gene aligns
      in expression according to samples sequenced with both NanoString and RNA-seq (NanoString/RNA Correlation tab)")), 
      br(),
      selectizeInput(
        "gene",
        "Enter Gene Name:",
        choices = NULL,
        options = list(
          placeholder = "RRM2",
          maxOptions = 20
        )
      ),
      actionButton("go", "Run"),
      br(),
      br(),
      p("Please cite the following:"),
      
      width = 3
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Cohort Descriptions", 
        p("This meta-cohort includes data from 7 studies with paired transcriptomics data pre and post NACT treatment. The studies include:"),
        p(HTML("Cohort Table:")),
        p(HTML("<strong>NanoString:</strong>
        <ul>
          <li>Manso: (GEO=GSE181597)</li>
          <li>James: (GEO=GSE201600)</li>
          <li>Bitler: 36 patients (GEO=)</li>
        </ul>")),
        p(HTML("<strong>RNA-seq:</strong>
        <ul>
          <li>Adzibolosu: (GEO=GSE227100, <a href='https://pubmed.ncbi.nlm.nih.gov/37435088/'>Adzibolosu et al<\a>)</li>
          <li>Javellana: 28 patients (RPKM-normalized counts at Supp Table 5 in <a href='https://pubmed.ncbi.nlm.nih.gov/34737212/'>Javellana et al</a>)</li>
          <li>EGA: The raw sequencing files are accessible via the European Genome Archive (<a href='https://ega-archive.org/datasets/EGAD00001006456'>EGAD00001006456</a>)</li>
        </ul>")),
        p(HTML("<strong>scRNA-seq:</strong>
        <ul>
          <li>EGA: (counts from GEO=GSE165897 with metadata available from <a href='https://www.science.org/doi/10.1126/sciadv.abm1831'>Zhang et al</a> and <a href='https://pubmed.ncbi.nlm.nih.gov/38383551/'>Perkiö et al</a>.)</li>
        </ul>")), 
        p(HTML("<strong>Microarray:</strong>
        <ul>
          <li><a href='https://pubmed.ncbi.nlm.nih.gov/32483290/'>Jiménez-Sánchez et al</a>: (counts from GEO=GSE146963)</li>
        </ul>"))
        ),
        tabPanel("Pre & Post Expression\nPFS Binary",
                helpText("Normalized expression (RPKM for RNA-seq) before and after NACT treatment of paired-patients
                split by Progression-Free survival being >12 months or <=12 months. A linear regression is performed with the model and R^2 value.
                If NanoString Experiments are not included, it is because the gene is not in the NanoString panel."), 
                br(),
                div(style="float:right;",
                  downloadButton("download_exp", label=NULL, icon=icon("download"))
                  ),        
                 plotOutput("exprPlot") %>% withSpinner()
        ),
        tabPanel("Expression & PFS\nCorrelation",
                helpText("Normalized expression (RPKM for RNA-seq) before and after NACT treatment and their
                correlation to PFS (in months). A linear regression is performed with the model and R^2 value. 
                The horizontal line pinpoints the 12 month cutoff. "), 
                br(),
                #div(style="float:right;",
                #  downloadButton("download_exp", label=NULL, icon=icon("download"))
                #  ),   
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
         tabPanel("Single-cell & Trends",
                helpText("Counts per million of genes within each cell type (normalized per cell type) at Pre-NACT, Post-NACT, and combined"), 
                br(),
                plotOutput("SC_trends") %>% withSpinner()
                
        ), 
        tabPanel("NanoString/RNAseq Correlation",
          helpText(HTML("Scatterplots of normalized counts from NanoString (y-axis) and from RNA-seq using three 
                pre-processing approaches for RNA-seq. All RNA-seq was normalized using RPKM (TPM and MR showed no clear differences). 
                Data from FFPE tumor samples sequenced with both RNA-seq and NanoString.
                    <ul>
                      <li>Isoform_RPKM: Used counts over isoforms best matching the NanoString probe.</li>
                      <li>RNAIsoform_RPKM: Used counts over isoforms most highly expressed in the RNA-seq.</li>
                      <li>Exon_RPKM: Using counts over the exon(s) aligning with the NanoString probes.</li>
                    </ul>
                Patients that show outlier expression values are labeled based on their patient id (e.g. P_9 and P_12).")), 
                br(),
                 plotOutput("Matched") %>% withSpinner()
                
        )
      )
    )
  )
)

server <- function(input, output, session) {

  ####################
  # Runs once on start up to provide gene options
  gene_list <-  union(rownames(datasets$exp_data$Adzib), rownames(datasets$exp_data$Jav))
  #cat(gene_list[1:2])
  updateSelectizeInput(
    session,
    "gene",
    choices = gene_list,
    selected = "RRM2",
    server = TRUE
  )

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
  
  ######################
  ### DATA DOWNLOAD ####
  # write out the expression data
  output$download_exp <- downloadHandler(
  filename = function() {
    paste0(selected_gene(), "_exp_data_", Sys.Date(), ".csv")
  },
  content = function(file) {
    df <- plot_df()
    validate(need(!is.null(df), "No data to download"))
    write.csv(df, file, row.names = FALSE)
  }
)

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
  width = 800, height=800)

  output$l2fcCorrPlot <- renderPlot({
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
      text(0.5, 0.5, paste("Gene", gene, "not included in NanoString panel (360)"))
      return()
    }
    
  }, height=400)

  # Plot the SC data
  output$SC_trends <- renderPlot({
    gene = selected_gene()
    graph_gene_prepost_heatmap(datasets$sc_data, gene)
  }, height=300)
  
}


shinyApp(ui, server)