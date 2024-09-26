library(shiny)
library(shinydashboard)
library(shinyjs)

ui <- dashboardPage(
  dashboardHeader(title = "GenomeProt",
                  dropdownMenu(type = "messages",
                               tags$li(HTML('<li><a href="https://biomedicalsciences.unimelb.edu.au/sbs-research-groups/physiology/Parker-laboratory-Metabolic-Proteomics" target="_blank"><i class="fa fa-user"></i><h4>About us</h4><p>Parker Laboratory</p></a></li>')),
                               tags$li(HTML('<li><a href="https://biomedicalsciences.unimelb.edu.au/sbs-research-groups/anatomy-and-physiology-research/stem-cell-and-developmental-biology/clark-lab" target="_blank"><i class="fa fa-user"></i><h4>About us</h4><p>Clark Laboratory</p></a></li>')),
                               tags$li(HTML('<li><a href="mailto:genomeprot@outlook.com" target="_blank"><i class="fa fa-question"></i><h4>Support</h4><p>genomeprot@outlook.com</p></a></li>'))
                  )),
  # tabs
  dashboardSidebar(
    sidebarMenu(menuItem("Welcome", tabName = "welcome", icon = icon("house")),
                menuItem("Generate database", tabName = "db_generation", icon = icon("database")),
                # commenting out proteomics for now
                #menuItem("Analyse MS proteomics", tabName = "analyse_proteomics", icon = icon("gear")),
                menuItem("Integrate data", tabName = "integration", icon = icon("code-merge")),
                menuItem("Visualise results", tabName = "visualisation", icon = icon("chart-bar"))
    )
  ),
  # body
  dashboardBody(
    useShinyjs(),  # shinyjs
    tags$head(
      tags$style(HTML("
        .spinner {
          margin: 0 auto;
          width: 30px;
          height: 30px;
          border: 6px solid #ccc;
          border-top: 6px solid #333;
          border-radius: 50%;
          animation: spin 1s linear infinite;
        }
  
        @keyframes spin {
          0% { transform: rotate(0deg); }
          100% { transform: rotate(360deg); }
        }
  
        .loading-container {
          display: none;
          text-align: center;
          margin-top: 20px;
        }
        
        #downloadResults {
          background-color: #4CAF50; /* Green */
          border: none;
          color: white;
          padding: 15px 32px;
          text-align: center;
          text-decoration: none;
          display: inline-block;
          font-size: 12px;
        }
        
        #downloadResults:disabled {
          background-color: #d3d3d3; /* Gray */
          color: #a9a9a9; /* Dark gray */
        }
        
        .spacing {
          margin-top: 20px;
        }
      ")),
      tags$script(HTML("
        Shiny.addCustomMessageHandler('disableButton', function(params) {
          var button = document.getElementById(params.id);
          button.disabled = true;
          button.style.backgroundColor = 'grey';
          button.style.borderColor = 'grey';
          document.getElementById(params.spinnerId).style.display = 'block';
        });
  
        Shiny.addCustomMessageHandler('enableButton', function(params) {
          var button = document.getElementById(params.id);
          button.disabled = false;
          button.style.backgroundColor = '';
          button.style.borderColor = '';
          document.getElementById(params.spinnerId).style.display = 'none';
        });
      "))
    ),
    tabItems(
      tabItem(tabName = "welcome",
              fluidRow(
                column(12,
                       div(class = "box box-primary", style = "padding-right: 5%; padding-left: 5%; font-size:110%", 
                           div(class = "box-body", shiny::includeMarkdown("welcome-page-text.md")),
                           img(src = "images/workflow.png", width = "100%"),
                       )
                )
              )
      ),
      tabItem(tabName = "db_generation", 
              h2("Generate a custom proteogenomics database"),
              h5("Creates an amino acid FASTA of all ORFs in your data to use as input for FragPipe/MaxQuant etc."),
              fluidRow(
                column(6,
                       
                       # Choices
                       radioButtons("sequencing_type", h5(tags$b("Select sequencing type:")),
                                    choices = c("Long-read (ONT, PacBio)" = "long-read", 
                                                "Short-read" = "short-read")),
                       radioButtons("input_type", h5(tags$b("Select input type:")),
                                    choices = c("FASTQs" = "fastq_input",
                                                "BAMs" = "bam_input",
                                                "GTF (and/or transcript counts)" = "gtf_input")),
                       
                       # Constant options
                       selectInput("organism", label = "Organism:", 
                                   choices = list("Roundworm (C. elegans)" = "celegans", 
                                                  "Fruit fly (D. melanogaster)" = "drosophila", 
                                                  "Human (H. sapiens)" = "human", 
                                                  "Mouse (M. musculus)" = "mouse", 
                                                  "Rat (R. rattus)" = "rat", 
                                                  "Zebrafish (D. rerio)" = "zebrafish"), 
                                   selected = "human"),
                       selectInput("database_type", label = "ORFs to be included in proteomedb:", 
                                   choices = list("canonical", "all"),
                                   selected = "all"),
                       numericInput("min_orf_length", 
                                    label = "ORF length (amino acids):", 
                                    value = 30),
                       h5(tags$b("Find short (10 to 'ORF length' amino acids) ORFs in UTRs of reference transcripts:")),
                       checkboxInput("user_find_utr_5_orfs", label = "Upstream 5' ORFs",
                                     value = FALSE, width = NULL),
                       checkboxInput("user_find_utr_3_orfs", label = "Downstream 3' ORFs",
                                     value = FALSE, width = NULL),
                       numericInput("minimum_tx_count", 
                                    label = "Minimum expression threshold (sum per transcript):", 
                                    value = 5),
                       fileInput("user_reference_gtf", "Upload reference annotation GTF:", NULL, buttonLabel = "Browse...", multiple = FALSE),
                       
                       # Variable options
                       conditionalPanel(
                         condition = "input.input_type == 'fastq_input'",
                         numericInput("user_threads", label = "CPUs:", value = 4, min = 1, max = 46, step = 1),
                         #h5("Map FASTQs, identify (in long-reads) and quantify isoforms, and generate the database"),
                         fileInput("user_reference_genome", "Upload reference genome FASTA:", NULL, buttonLabel = "Browse...", multiple = FALSE),
                         conditionalPanel(condition = "input.sequencing_type == 'short-read'",
                                          fileInput("transcriptome_file", "Upload reference transcriptome FASTA:", NULL, buttonLabel = "Browse...", multiple = FALSE)
                         ),
                         fileInput("user_fastq_files", "Upload FASTQ file(s):", NULL, buttonLabel = "Browse...", multiple = TRUE)
                       ),
                       conditionalPanel(
                         condition = "input.input_type == 'bam_input'",
                         numericInput("user_threads", label = "CPUs:", value = 4, min = 1, max = 46, step = 1),
                         conditionalPanel(condition = "input.sequencing_type == 'short-read'",
                                          fileInput("user_reference_genome_bam", "Upload reference genome FASTA:", NULL, buttonLabel = "Browse...", multiple = FALSE)),
                         fileInput("user_bam_files", "Upload BAM file(s):", NULL, buttonLabel = "Browse...", multiple = TRUE)
                       ),
                       conditionalPanel(
                         condition = "input.input_type == 'gtf_input' & input.sequencing_type == 'long-read'",
                         fileInput("user_gtf_file", "Upload 'bambu_transcript_annotations.gtf':", NULL, buttonLabel = "Browse...", multiple = FALSE),
                         fileInput("user_tx_count_file", "Upload 'bambu_transcript_counts.txt' (optional):", NULL, buttonLabel = "Browse...", multiple = FALSE)
                       ),
                       conditionalPanel(
                         condition = "input.input_type == 'gtf_input' & input.sequencing_type == 'short-read'",
                         fileInput("user_tx_count_file", "Upload transcript counts:", NULL, buttonLabel = "Browse...", multiple = FALSE)
                       ),
                       actionButton("db_submit_button", "Submit", class = "btn btn-primary")
                ),
                column(6,
                       HTML("<h3>Download your results:</h3>"),
                       downloadButton("db_download_button", "Download results (zip)", disabled = TRUE, style = "width:70%;"), # initially disabled
                       div(id = "db-loading-container", class = "loading-container", div(class = "spinner"))
                )
              )
      ),
      # tabItem(tabName = "analyse_proteomics", 
      #         h2("Run MetaMorpheus with your custom proteogenomics database to analyse MS proteomics data"),
      #         h5("NOTE: this step requires significant computation and time (>8 CPUs and high memory requirements)"),
      #         fluidRow(
      #           column(4,
      #                  selectInput("protease", label = "Protease:", 
      #                              choices = list("trypsin" = "trypsin"), 
      #                              selected = "trypsin"),
      #                  numericInput("mm_cpu", 
      #                               label = "CPUs", 
      #                               value = 1),
      #                  fileInput("user_mm_fasta", "Upload 'proteome_database.fasta'", NULL, buttonLabel = "Browse...", multiple = FALSE),
      #                  fileInput("user_mm_data", "Upload mzML/raw file(s):", NULL, buttonLabel = "Browse...", multiple = TRUE),
      #                  actionButton("proteomics_submit_button", "Submit", class = "btn btn-primary")
      #           ),
      #           column(6,
      #                  HTML("<h3>Download your results:</h3>"),
      #                  downloadButton("proteomics_download_button", "Download results (zip)", disabled = TRUE, style = "width:70%;"), # initially disabled
      #                  div(id = "proteomics-loading-container", class = "loading-container", div(class = "spinner"))
      #           )
      #         )
      # ),
      tabItem(tabName = "integration", 
              h2("Integrate proteomics results with transcriptomics"),
              h5("Creates BED12s and GTFs of peptides, ORFs and transcripts for visualisation and produces summary data"),
              fluidRow(
                column(6,
                       fileInput("user_proteomics_file", "Upload proteomics results:", NULL, buttonLabel = "Browse...", multiple = FALSE),
                       fileInput("user_fasta_file", "Upload 'proteome_database.fasta':", NULL, buttonLabel = "Browse...", multiple = FALSE),
                       fileInput("user_metadata_file", "Upload 'proteome_database_metadata.txt':", NULL, buttonLabel = "Browse...", multiple = FALSE),
                       fileInput("user_post_gtf_file", "Upload 'proteome_database_transcripts.gtf':", NULL, buttonLabel = "Browse...", multiple = FALSE),
                       actionButton("integ_submit_button", "Submit", class = "btn btn-primary")
                ),
                column(6,
                       HTML("<h3>Download your results:</h3>"),
                       downloadButton("integ_download_button", "Download results (zip)", disabled = TRUE, style = "width:70%;"), # initially disabled
                       div(id = "integ-loading-container", class = "loading-container", div(class = "spinner"))
                )
              )
      ),
      tabItem(tabName = "visualisation", 
              h2("Visualise results"),
              h5("Plots your results using the GTFs created in the integration module."),
              fluidRow(
                column(4,
                       fileInput("user_vis_gtf_file", "Upload 'combined_annotations.gtf' file:", NULL, buttonLabel = "Browse...", multiple = FALSE),
                       fileInput("user_vis_tx_count_file", "Upload 'bambu_transcript_counts.txt' (optional):", NULL, buttonLabel = "Browse...", multiple = FALSE),
                       fileInput("user_pep_count_file", "Upload peptide intensities file (optional):", NULL, buttonLabel = "Browse...", multiple = FALSE),
                       actionButton("vis_submit_button", "Submit", class = "btn btn-primary")
                ),
                column(8,
                       selectInput("gene_selector", "Select Gene", choices = NULL),
                       checkboxInput("uniq_map_peptides", "Only display genes encoding ORFs with uniquely mapped peptides", value = FALSE),
                       div(id = "vis-loading-container", class = "loading-container", div(class = "spinner")),
                       plotOutput("plot"),
                       downloadButton("vis_download_button", "Download plot", disabled = TRUE, class = "spacing") # initially disabled
                )
              )
      )
    )
  ),  
  skin = "purple"
)
