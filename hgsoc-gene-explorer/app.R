library(shiny)
library(shinycssloaders)
library(dplyr)

testing=FALSE
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
  
  titlePanel("HGSOC Multi-Transcriptomics\nGene Explorer"),

  sidebarLayout(
    sidebarPanel(
      p(HTML("<strong>Click Run</strong> to check out how your favorite gene looks in a large meta-cohort of 
      HGSOC patients with paired transcriptomics samples pre and post NACT (more info in Cohort Descriptions).")),
      p(HTML("If the gene is found in the NanoString panel (nCounter PanCancer IO 360 panel), you can also see how well the gene aligns
      in expression according to 24 samples sequenced with both NanoString and RNA-seq (NanoString/RNA Correlation tab)")), 
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
      #downloadButton("download_exp", label="Download Expression Data", icon=icon("download")),
      br(),
      p(HTML("<strong>Download</strong> the full data <a href='https://github.com/Hope2925/NanoString-RNAseq-HGSOC/tree/main/hgsoc-gene-explorer/data'>here</a> with explanations of the data types found <a href='https://github.com/Hope2925/NanoString-RNAseq-HGSOC/blob/main/README.md#data'>here</a> or gene-specific data in the Tabs")),
      p(HTML("<strong>Important Notes:</strong> 
      <ul>
          <li>Microarray tends to show a lower dynamic range and therefore is graphed separately from NanoString/RNA-seq</li>
          <li>Some genes log2FC are not calculated due to hitting limits of detection and therefore being unreliable. Expression is still available.</li>
        </ul>"
      
      )),
      p("Please cite the following:"),
      
      width = 3
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Cohort Descriptions", 
        p("This meta-cohort includes data from 7 studies with paired transcriptomics data pre and post NACT treatment along with Platinum Free Survival (PFS). The studies include:"),
        p(HTML("Cohort Table:")),
        imageOutput("cohort_img"),
        #tags$img(src = "Cohort_Table_2.png", width = "70%"),
        #HTML('<img src="Cohort_Table.png"/>'),
        # height="100px" width="100px"
        #tags$img(src = "Cohort_Table.png", 
        #   height = "100px", 
        #   width = "300px", 
        #  alt = "Cohort Table"),
        p(HTML("<strong>NanoString:</strong>
        <ul>
          <li>Manso: 17 patients (GEO=GSE181597, <a href='https://www.nature.com/articles/s41698-021-00247-3'>Lodewijk et al</a>)</li>
          <li>James: 20 patients with paired pre/post (31 in original dataset when not requiring pairing) (GEO=GSE201600, <a href='https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.965331/full'>James et al</a>)</li>
          <li>Bitler: 35 patients (GEO=TBA, <a href='https://aacrjournals.org/clincancerres/article/26/23/6362/83054/The-Capacity-of-the-Ovarian-Cancer-Tumor'>Jordan et al</a> and LUCY PAPER)</li>
        </ul>")),
        p(HTML("<strong>RNA-seq:</strong>
        <ul>
          <li>Adzibolosu: 15 patients with PFS provided by original authors (GEO=GSE227100, <a href='https://pubmed.ncbi.nlm.nih.gov/37435088/'>Adzibolosu et al</a>)</li>
          <li>Javellana: 28 patients (RPKM-normalized counts at Supp Table 5 in <a href='https://pubmed.ncbi.nlm.nih.gov/34737212/'>Javellana et al</a>)</li>
          <li>EGA: A total of 22 patients had the PFS information and paired bulk sequencing data. The exact samples used can be found in Supplemental Table 1 of X. The raw sequencing files are accessible via the European Genome Archive (<a href='https://ega-archive.org/datasets/EGAD00001006456'>EGAD00001006456</a>)</li>
        </ul>")),
        p(HTML("<strong>scRNA-seq:</strong>
        <ul>
          <li>EGA: 11 patients with paired Pre and Post, (counts from GEO=GSE165897 with metadata available from <a href='https://www.science.org/doi/10.1126/sciadv.abm1831'>Zhang et al</a> and <a href='https://pubmed.ncbi.nlm.nih.gov/38383551/'>Perkio et al</a>.)</li>
        </ul>")), 
        p(HTML("<strong>Microarray:</strong>
        <ul>
          <li>JimSanchez: 28 patients from <a href='https://pubmed.ncbi.nlm.nih.gov/32483290/'>Jimenez-Sanchez et al</a>: (counts from GEO=GSE146963)</li>
        </ul>"))
        ),
        tabPanel("Pre & Post Expression\nPFS Binary",
                helpText("Normalized expression (RPKM for RNA-seq) before and after NACT treatment of paired-patients
                split by Progression-Free survival being >12 months or <=12 months. A linear regression is performed with the model and R^2 value.
                If NanoString Experiments are not included, it is because the gene is not in the NanoString panel."), 
                br(),
                #div(style="float:right;",
                downloadButton("download_exp", label="Download Data", icon=icon("download")),
                 # ),        
                 plotOutput("exprPlot") %>% withSpinner()
        ),
        tabPanel("Expression & PFS\nCorrelation",
                helpText("Normalized expression (RPKM for RNA-seq) before and after NACT treatment and their
                correlation to PFS (in months). A linear regression is performed with the model and R^2 value. 
                The horizontal line pinpoints the 12 month cutoff. "), 
                br(),
                helpText("Download data in Tab Pre & Post Expression PFS Binary"),
                 plotOutput("exprCorrPlot") %>% withSpinner()
                
        ), 
        tabPanel("Log2FC\nPFS Binary",
                 helpText("Log2 Fold Change between Post and Pre NACT treatment across patients
                 split by Progression-Free survival being >12 months or <=12 months. 
                   T-tests are performed for each of the experiments and if combining all experiments, with p-values shown"), 
                br(),
                downloadButton("download_l2fc", label="Download Data", icon=icon("download")),
                 plotOutput("l2fcPlot") %>% withSpinner()   
        ),
        tabPanel("Log2FC & PFS\nCorrelation",
                helpText("Log2 Fold Change between Post and Pre NACT treatment and their correlation to PFS (in months). 
                   A linear regression is performed with the model and R^2 value.  
                   The horizontal line pinpoints the 12 month cutoff. "), 
                br(),
                helpText("Download data in Tab Log2FC PFS Binary"),
                plotOutput("l2fcCorrPlot") %>% withSpinner()
                
        ), 
         tabPanel("Single-cell & Trends",
                helpText(HTML("Counts per million of genes within each cell type (normalized per cell type) at Pre-NACT, Post-NACT, and combined",
                "Cell Type Descriptions can be found below or <a href='https://www.science.org/doi/10.1126/sciadv.abm1831'>here</a><ul>
                  <li>EOC_C1 & EOC_C2: No clear overrepresented genes</li>
                  <li>EOC_C3: patient specific & mostly EMT (TGFB, focal adhesion) + Interferon signaling + RNA processing</li>
                  <li>EOC_C4: Differentiated (O-linked glycolyslation of mucins) + Interferon signaling (STAT2, OAS1)</li>
                  <li>EOC_C5: DNA repair & Cell cycle + Proteasomal degradation + TCA cycle</li>
                  <li>EOC_C6: No clear overrepresented genes</li>
                  <li>EOC_C7: Stress-associated (IL6, TNF)</li>
                  <li>EOC_C8: TCA cycle</li>
                  <li>EOC_C9: patient specific + Cytokine + Apoptosis</li>
                  <li>EOC_C10: patient specific + Antigen presentation + TCA Cycle</li>
                  <li>EOC_C11: RNA processing</li>
                  <li>EOC_C12: No clear overrepresented genes</li>
                  <li>CAF-1: Matrix metalloproteinases (MMPs)</li>
                  <li>CAF-2: Inflammatory CAF (iCAF) markers like IL6, CXCL12, LIF; enriched in stress-high tumors</li>
                  <li>CAF-3: Myofibroblast markers</li>

                </ul>")), 
                br(),
                downloadButton("download_SC", label="Download Data", icon=icon("download")),
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

  ###################
  # Image output since for some reason the imageOutput in the UI wasn't working with the tags$img approach
   output$cohort_img <- renderImage({
    if (testing) {
      file_path <- normalizePath("hgsoc-gene-explorer/www/Cohort_Table_2.png")  # absolute path
    } else {
      file_path <- normalizePath("www/Cohort_Table_2.png")  # absolute path
    }
    list(
      src = file_path,
      contentType = "image/png",
      width = 600,
      alt = "Cohort Table"
    )
  }, deleteFile = FALSE)

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

  # have reactive variables
  selected_gene <- reactiveVal(NULL)
  nano_gene     <- reactiveVal(NULL)
  plot_df    <- reactiveVal(NULL)
  plot_l2fc_df <- reactiveVal(NULL)
  plot_SC_df <- reactiveVal(NULL)

  # If a new input is provided, update the following
  observeEvent(input$go, {
    # require that the input have a gene
  req(input$gene)

  gene <- sanitize_gene(input$gene)

  is_nano <- gene %in% rownames(datasets$exp_data$Bitler)
  

  df <- get_plot_df(
    gene,
    datasets$exp_data,
    datasets$meta,
    is_nano
  )
  
  l2fc_df = get_l2fc_df(gene, datasets$l2fc_data, is_nano) 

  SC_df = get_SC_df(gene, datasets$sc_data)

  # assign the variables
  selected_gene(gene)
  nano_gene(is_nano)
  plot_df(df)
  plot_l2fc_df(l2fc_df)
  plot_SC_df(SC_df)
})


  #selected_gene <- eventReactive(input$go, {
  #  req(input$gene)
  #  clean_gene <- sanitize_gene(input$gene)
  #  return(clean_gene)
  #})
  
  # determine if nano gene or not ONLY when gene changes
  #nano_gene <- eventReactive(selected_gene(),{
  #  gene <- selected_gene()
  #  gene %in% rownames(datasets$exp_data$Bitler)
  #})

  # get the plot df ONLY when gene changes
  #plot_df <- eventReactive(selected_gene(), {
   # get_matched_plot_df(selected_gene(), datasets$exp_data, datasets$meta, nano_gene())
  #})
  
  ######################
  ### DATA DOWNLOAD ####
  # write out the expression data
  output$download_exp <- downloadHandler(
  filename = function() {
    req(selected_gene())
    paste0(selected_gene(), "_exp_data.csv")
  },
  content = function(file) {
    df <- plot_df()
    validate(
      need(!is.null(df) && nrow(df) > 0, "Click Run to generate data.")
    )
    write.csv(df, file, row.names = FALSE)
  }
  )

  # l2fc data
  output$download_l2fc <- downloadHandler(
  filename = function() {
    req(selected_gene())
    paste0(selected_gene(), "_l2fc_data.csv")
  },
  content = function(file) {
    req(plot_l2fc_df())
    #validate(need(!is.null(df), "No data to download"))
    write.csv(plot_l2fc_df(), file, row.names = FALSE)
  }
  )

  # SC data
  output$download_SC <- downloadHandler(
  filename = function() {
    req(selected_gene())
    paste0(selected_gene(), "_SCprepost_CPM.csv")
  },
  content = function(file) {
    req(plot_SC_df())
    #validate(need(!is.null(df), "No data to download"))
    write.csv(plot_SC_df(), file, row.names = FALSE)
  }
  )

  ################
  ## PLOTTING ####
  
  output$exprPlot <- renderPlot({
    validate(
      need(nrow(plot_df()) > 0, "Click Run to generate data.")
    )
    # plot the expression with PFS binary
    plot_matched_gene(selected_gene(), plot_df())
  },
  height = 900)

  output$exprCorrPlot <- renderPlot({
    validate(
      need(nrow(plot_df()) > 0, "Click Run to generate data.")
    )
    # plot the correlation with expression
    plot_exp_corr_PFS(selected_gene(), plot_df())
  },
  height = 800)

  output$l2fcPlot <- renderPlot({
    validate(
      need(nrow(plot_l2fc_df()) > 0, "Click Run to generate data.")
    )
    # plot the L2FC split by PFS binary
    plot_l2fc(selected_gene(), plot_l2fc_df())
  },
  width = 800, height=800)

  output$l2fcCorrPlot <- renderPlot({
    # plot hte correlation with Log2FC
    validate(
      need(nrow(plot_l2fc_df()) > 0, "Click Run to generate data.")
    )
    plot_l2fc_corr_PFS(selected_gene(), plot_l2fc_df())
  }, 
  height=1000)

  
  # plot the correlation between NanoString and RNA-seq
  output$Matched <- renderPlot({
    is_nano = nano_gene()
    gene = selected_gene()
    validate(
      need(!is.null(selected_gene()), "Click Run to generate data.")
    )
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
    validate(
      need(nrow(plot_SC_df()) > 0, "Click Run to generate data.")
    )
    graph_gene_prepost_heatmap(plot_SC_df(), selected_gene())
  }, height=300)
  
}


shinyApp(ui, server)