library(shiny)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(InterMineR)
library(msaR)
library(ape)
library(ggtree)
library(limma)
library(scales)
library(shinycssloaders)
library(reshape2)
library(stringr)
library(NGLVieweR)
library(bio3d)
library(biomaRt)

# Environment variable associated with the container
# Set working directory to Docker path if it's set
docker <- Sys.getenv("REELGENE_DOCKER")
in_docker <- if (docker != "") TRUE else FALSE
if (in_docker) setwd("/reelgene")
addResourcePath("static", if (in_docker) "/reelgene/www" else "www")

# Sample data frame with GeneID, PanGeneID, Taxa, and Value columns
df <- read.csv('data/reelGene_allNAM_withConservation.csv')
#df <- reelGene
im <- initInterMine(mine=listMines()["MaizeMine"])
transcriptMatrix <- read.csv('data/inputMatrix/formatted_LSTM_metaTab_B73.txt')

# Define UI
ui <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "static/reelgene.css")
  ),
  titlePanel("Gene ID Lookup"),
  sidebarLayout(
    sidebarPanel(
      width = 8,
      textInput("geneID", "Enter Gene ID:", placeholder = "Gene ID"),
      actionButton("lookup", "Lookup"),
      # Add tableOutput for gene data
      tableOutput("geneData"),
      uiOutput("maizeGDB")
    ),
    mainPanel(
      width = 12,
      tabsetPanel(
        id = "tabs",
        tabPanel("Transcript Names and Scores",
                 h4("Transcript Names and Scores:"),
                 conditionalPanel(
                   condition = "input.lookup > 0",
                   withSpinner(tableOutput("transcripts"), color = "gold")
                 ),
                 conditionalPanel(
                   condition = "input.lookup > 0",
                   withSpinner(plotOutput("histogramPlot"), color = "gold")
                 ),
                 conditionalPanel(
                   condition = "input.lookup > 0",
                   withSpinner(plotOutput("scatterPlot"), color = "gold")
                 ),
                 conditionalPanel(
                   condition = "input.lookup > 0",
                   withSpinner(plotOutput("meanVarPlot"), color = "gold")
                 )
        ),
        tabPanel("Gene Model",
                 h4("Gene Model"),
                 conditionalPanel(
                   condition = "input.lookup > 0",
                   withSpinner(plotOutput("geneModelPlot"), color = "gold")
                 )
        ),
        tabPanel("PanGene Data",
                 h4("PanGene Data:"),
                 conditionalPanel(
                   condition = "input.lookup > 0",
                   withSpinner(plotOutput("boxplotPlot"), color = "gold")
                 ),
                 conditionalPanel(
                   condition = "input.lookup > 0",
                   withSpinner(tableOutput("panGeneData"), color = "gold")
                 )
        ),
        tabPanel("Expression Data",
                 h4("Expression Data:"),
                 conditionalPanel(
                   condition = "input.lookup > 0",
                   withSpinner(tableOutput("expressionData"), color = "gold")
                 )
        ),
        tabPanel("PanZea Variation (SNPs & QTL)",
                 h4("SNPs from PanZea"),
                 conditionalPanel(
                   condition = "input.lookup > 0",
                   withSpinner(tableOutput("panZeaSNPs"), color = "gold")
                 ),
                 h4("QTL from PanZea"),
                 conditionalPanel(
                   condition = "input.lookup > 0",
                   withSpinner(tableOutput("panZeaQTLs"), color = "gold")
                 )
        ),
        tabPanel("PanAnd Alignments",
                 fluidRow(
                   column(width = 12,
                          sidebarPanel(
                            selectInput(
                              "msa_file",
                              "Transcript:",
                              choices = "",
                              selected = "Unselected",
                              multiple = FALSE,
                              selectize = FALSE,
                              width = NULL,
                              size = NULL
                            ),
                            selectInput(
                              "filter",
                              "Filter:",
                              c(
                                "Unselected",
                                "Long Reads",
                                "Short Reads",
                                "Zea",
                                "Zea & Tripsacum"
                              ),
                              selected = "Unselected",
                              multiple = FALSE,
                              selectize = FALSE,
                              width = NULL,
                              size = NULL
                            )
                          )
                   )
                 ),
                 mainPanel(
                   width = 12,
                   h4("Multiple Sequence Alignment"),
                   conditionalPanel(
                     condition = "input.lookup > 0",
                     withSpinner(msaROutput("msa_viewer"), color = "gold")
                   ),
                   h4("Gene Tree"),
                   conditionalPanel(
                     condition = "input.lookup > 0",
                     withSpinner(plotOutput("tree_viewer"), color = "gold")
                   )
                 )
        ),
        tabPanel("NAM Alignments",
                 fluidRow(
                   column(width = 12,
                          sidebarPanel(
                            selectInput(
                              "msa_file_NAM",
                              "Transcript:",
                              choices = "",
                              selected = "Unselected",
                              multiple = FALSE,
                              selectize = FALSE,
                              width = NULL,
                              size = NULL
                            )
                          )
                   )
                 ),
                 mainPanel(
                   width = 12,
                   h4("Multiple Sequence Alignment"),
                   conditionalPanel(
                     condition = "input.lookup > 0",
                     withSpinner(msaROutput("msa_viewer_NAM"), color = "gold")
                   ),
                   h4("Gene Tree"),
                   conditionalPanel(
                     condition = "input.lookup > 0",
                     withSpinner(plotOutput("tree_viewer_NAM"), color = "gold")
                   )
                 )
        ),
        tabPanel(
          "Sequence",
          sidebarPanel(
            selectInput(
              "seq_choices",
              "Transcript:",
              choices = "",
              selected = "Unselected",
              multiple = FALSE,
              selectize = FALSE,
              width = NULL,
              size = NULL
            ),
            selectInput(
              "seqtype",
              "Sequence Type:",
              choices = c("Unselected", "cdna", "peptide", "3utr", "5utr", "gene_exon", "transcript_exon",
                          "transcript_exon_intron", "gene_exon_intron", "coding", "coding_transcript_flank",
                          "coding_gene_flank", "transcript_flank", "gene_flank"),
              selected = "Unselected",
              multiple = FALSE,
              selectize = FALSE,
              width = NULL,
              size = NULL
            )
          ),
          mainPanel(
            width = 12,
            h4("Sequence"),
            conditionalPanel(
              condition = "input.lookup > 0",
              withSpinner(
                tags$div(
                  style = "max-width: 100%; word-wrap: break-word;",
                  textOutput("sequence")
                ),
                color = "gold"
              )
            )
          )
        ),
        tabPanel("Protein Structure",
                 fluidRow(
                   column(width = 12,
                          sidebarPanel(
                            selectInput(
                              "protein_file",
                              "Isoform:",
                              choices = "",
                              selected = "Unselected",
                              multiple = FALSE,
                              selectize = FALSE,
                              width = NULL,
                              size = NULL
                            )
                          )
                   )
                 ),
                 mainPanel(
                   width = 12,
                   h4("Predicted Protein Structure"),
                   conditionalPanel(
                     condition = "input.lookup > 0",
                     withSpinner(NGLVieweROutput("protein_viewer"), color = "gold")
                   )
                 )
      )
    )
  )
),
tags$footer(paste0("Container: ", if (in_docker) docker else "?"))
)



# Define server
server <- function(input, output, session) {
  # Create a reactive value for available transcripts
  available_transcripts <- reactiveVal(character(0))
  
  observeEvent(input$lookup, {
    
    # OUTPUT GENE INFO
    gene_id <- isolate(input$geneID)
    queryInfo = getTemplateQuery(
      im = im,
      name = "gene_chromosomal_location"
    )
    queryInfo$where[[2]][["value"]] <- gene_id
    resInfo<- runQuery(im, queryInfo)
    output$geneData <- renderTable({
      if (nrow(resInfo) > 0) {
        colnames(resInfo) <- c("Primary Identifier", "Gene Name", "Biotype", "Source", "Length", "Chromosome", "End", "Start", "Strand", "Assembly")
        resInfo[,c(1:6,8,7,9:10)]
      } else {
        "No transcripts found for the given Gene ID."
      }
    })
    
    
    # OUTPUT TRANSCRIPT TABLE
    gene_id <- isolate(input$geneID)
    transcript_df <- df[df$Gene == gene_id, c("Transcript", "ProteinScore", "ExonScore", "average", "class", "conservation", "taxaShort", "Pan_gene_ID")]
    output$transcripts <- renderTable({
      if (nrow(transcript_df) > 0) {
        transcripts <- list(transcript_df$Transcript)
        available_transcripts(transcript_df$Transcript)
        transcript_df
      } else {
        "No transcripts found for the given Gene ID."
      }
    })
   
    
    # OUTPUT HISTOGRAM PLOT
     output$histogramPlot <- renderPlot({
      gene_id <- isolate(input$geneID)
      
      # Create the histogram plot
      p <- ggplot(df, aes(x = average)) +
        geom_histogram(fill = "lightblue", color = "black") +
        labs(x = "Average exon/protein model score", y = "Frequency", title = "Distribution of Scores of all transcripts")
      
      # Add vertical lines for the genes of interest
      if (gene_id %in% df$Gene) {
        gene_df <- df[df$Gene == gene_id, ]
        p <- p + geom_vline(aes(xintercept = average, color = Transcript, linetype = Transcript), 
                            data = gene_df, show.legend = FALSE) +
          geom_text(aes(x = average, y = 400000, label = Transcript), data = gene_df, vjust = -1.5, color = "black", size = 3)
      }
      p
    })
     
     
    # OUTPUT BOX PLOT
    output$boxplotPlot <- renderPlot({
      gene_id <- isolate(input$geneID)
      
      # Get the PanGeneID for the inputted gene
      pan_gene_id <- df$Pan_gene_ID[df$Gene == gene_id][1]
      
      # Filter the data for genes with the same PanGeneID
      plot_df <- df %>%
        filter(Pan_gene_ID == pan_gene_id)
      
      # Create the boxplot
      p <- ggplot(plot_df, aes(x = taxaShort, y = average)) +
        geom_boxplot() +
        theme(axis.text.x = element_text(angle = -90)) + ylim(0,1) + ylab("Average protein/exon score") + labs(title = "reelGene scores for pangene across all NAM") + xlab("NAM line")
      
      p
    })
    
    
    # OUTPUT PAN GENE DATA
    gene_id <- isolate(input$geneID)
    queryPanGene = getTemplateQuery(
      im = im,
      name = "Gene_to_all_syntelogs"
    )
    queryPanGene$where[[1]][["value"]] <- gene_id
    resPanGene <- runQuery(im, queryPanGene)
    output$panGeneData <- renderTable({
      if (nrow(resPanGene) > 0) {
        resPanGene
      } else {
        "No transcripts found for the given Gene ID."
      }
    })
    
    
    # OUTPUT SCATTER PLOT
    output$scatterPlot <- renderPlot({
      gene_id <- isolate(input$geneID)
      
      # Create the scatter plot with red points for the gene of interest
      gene_df <- df[df$Gene == gene_id, ]
      
      p <- ggplot(gene_df, aes(x = ExonScore, y = ProteinScore)) +
        geom_point(color = "red") +
        labs(x = "Exon Score", y = "Protein Score", title = "Exon Score vs Protein Score") + xlim(0,1) + ylim(0,1)
      
      p
    })
    
    
    # OUTPUT MEAN VARIANCE PLOT
    output$meanVarPlot <- renderPlot({
      gene_id <- isolate(input$geneID)
      
      # Get the PanGeneID for the inputted gene
      pan_gene_id <- df$Pan_gene_ID[df$Gene == gene_id][1]
      
      # Filter the data for genes with the same PanGeneID
      plot_df <- df %>%
        filter(Pan_gene_ID == pan_gene_id)
      
      # Calculate mean and variance for each PanGeneID
      mean_var_df <- df %>%
        group_by(Pan_gene_ID) %>%
        summarise(mean_score = mean(average),
                  var_score = var(average))
      
      # Create the scatter plot of mean and variance
      p <- ggplot(mean_var_df, aes(x = mean_score, y = var_score)) +
        geom_point(color = "gray", size = 3) +
        geom_point(data = mean_var_df[mean_var_df$Pan_gene_ID == pan_gene_id, ], 
                   aes(color = "red"), size = 3) +
        scale_color_manual(values = c("red" = "red"), guide = FALSE) +
        labs(x = "Mean Score", y = "Variance", title = "Mean Score vs Variance")
      
      p
    })
    
    
    # OUTPUT EXPRESSION DATA
    gene_id <- isolate(input$geneID)
    queryGeneExpress = getTemplateQuery(
      im = im,
      name = "Expression_info_reportPage"
    )
    queryGeneExpress$where[[1]][["value"]] <- gene_id
    resGeneExpress <- runQuery(im, queryGeneExpress)
    output$expressionData <- renderTable({
      if (nrow(resGeneExpress) > 0) {
        resGeneExpress
      } else {
        "No transcripts found for the given Gene ID."
      }
    })

    
    # OUTPUT SNP DATA
    gene_id <- isolate(input$geneID)
    querySNP = getTemplateQuery(
      im = im, 
      name = "SNV_alleles_and_consequences"
    )
    querySNP$where[[1]][["path"]] <- "SNV.consequences.transcript.gene.primaryIdentifier"
    querySNP$where[[1]][["value"]] <- gene_id
    querySNP$joins <- "SNV.chromosomeLocation"
    querySNP$select <- c("SNV.primaryIdentifier","SNV.chromosome.assembly","SNV.chromosome.primaryIdentifier","SNV.chromosomeLocation.start","SNV.chromosomeLocation.end","SNV.referenceAllele","SNV.alternateAllele","SNV.consequences.consequenceTypes.name","SNV.consequences.transcriptIdentifier","SNV.consequences.referenceResidue","SNV.consequences.alternateResidue","SNV.consequences.siftNumericalValue","SNV.consequences.siftQualitativePrediction")
    resSNP<- runQuery(im, querySNP)
    output$panZeaSNPs <- renderTable({
      if (nrow(resSNP) > 0) {
        resSNP
      } else {
        "No transcripts found for the given Gene ID."
      }
      })
    
    
    # OUTPUT GENE MODEL PLOT

    # Calculate the height dynamically based on the number of transcripts
    num_transcripts <- nrow(transcript_df)
    h <- 200 + (num_transcripts * 250)  # Adjust the factor (50) to control the height
    
    output$geneModelPlot <- renderPlot({
      
      # Create the scatter plot with red points for the gene of interest
      
      workingMatrix <- transcriptMatrix %>% filter(Gene == gene_id)
      meltedMatrix <- melt(workingMatrix[, c(1:6, 11:14)], id.vars = c(1:5))
      meltedMatrix$value[meltedMatrix$value == -1] <- NA
      meltedMatrix$scale <- ifelse(substr(meltedMatrix$feature, 1, 1) == "C", 14, ifelse(substr(meltedMatrix$feature, 1, 1) == "I", 2, 6))
      
      # Create a new dataframe for arrows
      arrow_df <- meltedMatrix %>%
        group_by(transcript, variable) %>%
        summarise(arrow_pos = (start + end) / 2, 
                  arrow_label = ifelse(strand == "+", ">", "<")) %>%
        ungroup()
      
      # Plot the gene models with arrows
      p <- ggplot(meltedMatrix) +
        geom_segment(aes(x = start, y = factor(variable), xend = end, yend = factor(variable), color = value, size = scale)) +
        geom_text(aes(x = (start + end) / 2, y = 0, label = feature), vjust = -1.5, size = 3, angle = 45) +
        scale_color_gradient(low = "blue", high = "red", na.value = "gray") +
        labs(x = "Genomic Position", y = "Transcript Scores") +
        theme_minimal() +
        scale_size_area(max_size = 10) +
        facet_wrap(~ transcript, ncol = 1) +
        geom_text(data = arrow_df, aes(x = arrow_pos, y = factor(variable), label = arrow_label),
                  size = 4, color = "black", fontface = "bold")

      p
    }, height = h) # dynamically based on number of transcipts
    
    # MSA AND TREE DATA
    
    # Generate a list of choices for the selectInput
    transcript_choices <- c("Unselected")
    if (nrow(transcript_df) > 0) {
      transcript_choices <- c("Unselected", setNames(paste0("data/msas/panand/", transcript_df$Transcript, "_msa_padded10N_rc.fa"), transcript_df$Transcript))
    }
    
    # Update the choices of the selectInput for "Transcript"
    updateSelectInput(
      session = session,
      inputId = "msa_file",
      label = "Transcript:",
      choices = transcript_choices,
      selected = "Unselected"
    )
    
    # Reactive function to load the MSA file
    msa_data <- reactive({
      req(input$msa_file)
      if (input$msa_file != "Unselected") {
        ape::read.FASTA(input$msa_file)
      }
    })
    
    # Define a reactive value to hold the filtered MSA data
    filtered_msa <- reactiveVal(NULL)
    
    # Function to filter the MSA data based on sequence names
    filter_msa <- function(sequence_names) {
      if (is.null(sequence_names) || sequence_names == "Unselected") {
        return(msa_data())
      } else if (sequence_names == "Long Reads") {
        pattern <- "^..-|^(B73)"
        list <- names(msa_data())[grepl(pattern, names(msa_data()), ignore.case = FALSE)]
        filtered_msa_data <- subset(msa_data(), names(msa_data()) %in% list)
        return(filtered_msa_data)
      } else if (sequence_names == "Short Reads") {
        pattern <- c("^..-", "^Z", "^Td")
        list <- names(msa_data())[!grepl(paste(pattern, collapse = "|"), names(msa_data()), ignore.case = FALSE)]
        filtered_msa_data <- subset(msa_data(), names(msa_data()) %in% list)
        return(filtered_msa_data)
      } else if (sequence_names == "Zea") {
        pattern <- "^Z|^(B73)"
        list <- names(msa_data())[grepl(pattern, names(msa_data()), ignore.case = FALSE)]
        filtered_msa_data <- subset(msa_data(), names(msa_data()) %in% list)
        return(filtered_msa_data)
      } else if (sequence_names == "Zea & Tripsacum") {
        pattern <- "^Td|^(B73)|^Z"
        list <- names(msa_data())[grepl(pattern, names(msa_data()), ignore.case = FALSE)]
        filtered_msa_data <- subset(msa_data(), names(msa_data()) %in% list)
        return(filtered_msa_data)
      }
    }
    
    # Automatically load the MSA file when a new option is chosen
    observeEvent(input$msa_file, {
      filtered_msa_data <- msa_data()
      filtered_msa(filtered_msa_data)
    }, ignoreInit = TRUE)
    
    # Apply the filter when the "Filter" dropdown changes
    observe({
      filtered_msa_data <- filter_msa(input$filter)
      filtered_msa(filtered_msa_data)
    })
    
    # Render the MSA viewer
    output$msa_viewer <- renderMsaR({
      msaR(menu = TRUE, filtered_msa(), seqlogo = FALSE, alignmentHeight = 550, conservation = TRUE)
    })
    
    # Render the tree viewer
    output$tree_viewer <- renderPlot({
      if (!is.null(msa_data())) {
        alnm <- as.matrix(msa_data())
        alnm.filtered <- alnm[!apply(alnm, 1, function(x) all(x == "04")), !alnm[2, ] == "04"]
        d <- dist.dna(alnm.filtered, model = "K80", pairwise.deletion = TRUE)
        tre <- njs(d)
        ggtree(tre,size =1.1,ladderize = T ) + 
          geom_tiplab(hjust = .1,as_ylab = T)+
          guides(shape = 19) +
          theme(legend.position = c(.1,.95),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 12)) 
      }
    })

    
    
    # MSA AND TREE DATA - NAM
    
    # Generate a list of choices for the selectInput
    transcript_choices_NAM <- c("Unselected")
    if (nrow(transcript_df) > 0) {
      transcript_choices_NAM <- c("Unselected", setNames(paste0("data/msas/nam/", transcript_df$Transcript, "_msa.fastA"), transcript_df$Transcript))
    }
    
    # Update the choices of the selectInput for "Transcript"
    updateSelectInput(
      session = session,
      inputId = "msa_file_NAM",
      label = "Transcript:",
      choices = transcript_choices_NAM,
      selected = "Unselected"
    )
    
    # Reactive function to load the MSA file
    msa_data_NAM <- reactive({
      req(input$msa_file_NAM)
      if (input$msa_file_NAM != "Unselected") {
        ape::read.FASTA(input$msa_file_NAM)
      }
    })
    
    
    # Automatically load the MSA file when a new option is chosen
    observeEvent(input$msa_file_NAM, {
      msa_data_NAM <- msa_data_NAM()
    }, ignoreInit = TRUE)
    
    # Render the MSA viewer
    output$msa_viewer_NAM <- renderMsaR({
      msaR(menu = TRUE, msa_data_NAM(), seqlogo = FALSE, alignmentHeight = 550, conservation = TRUE)
    })
    
    # Render the tree viewer
    output$tree_viewer_NAM <- renderPlot({
      if (!is.null(msa_data_NAM())) {
        alnm <- as.matrix(msa_data_NAM())
        alnm.filtered <- alnm[!apply(alnm, 1, function(x) all(x == "04")), !alnm[2, ] == "04"]
        d <- dist.dna(alnm.filtered, model = "K80", pairwise.deletion = TRUE)
        tre <- njs(d)
        ggtree(tre,size =1.1,ladderize = T ) + 
          geom_tiplab(hjust = .1,as_ylab = T)+
          guides(shape = 19) +
          theme(legend.position = c(.1,.95),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 12)) 
      }
    })

    # Generate a list of choices for the selectInput
    seq_choices <- c("Unselected")
    if (nrow(transcript_df) > 0) {
      seq_choices <- c("Unselected",
          setNames(transcript_df$Transcript, transcript_df$Transcript
          ))
    }
    
    # OUTPUT SEQUENCE DATA
    # Update the choices of the selectInput for "Sequence"
    updateSelectInput(
      session = session,
      inputId = "seq_choices",
      label = "Transcript:",
      choices = seq_choices,
      selected = "Unselected"
    )
    
    # render the sequence output
    output$sequence <- renderText({
      if (input$seq_choices != "Unselected" && input$seqtype != "Unselected") {
          maize_mart <- useMart(biomart = "plants_mart", host = "https://plants.ensembl.org", dataset = "zmays_eg_gene")
          maize_mart@biomart <- 'ENSEMBL_MART_ENSEMBL'
          seq <- getSequence(id = input$seq_choices, mart = maize_mart, seqType = input$seqtype, type = "ensembl_transcript_id")
          paste(seq)
      }
    })

    ######### Output protein plots
    # Generate a list of choices for the selectInput
    isoform_choices <- c("Unselected")
    if (nrow(transcript_df) > 0) {
      df_modified <- transcript_df %>%
        mutate(across('Transcript', str_replace, '_T', '_P'))
      isoform_choices <- c("Unselected", setNames(paste0("data/proteinStructure/", df_modified$Transcript, ".pdb"), df_modified$Transcript))
    }
    
    # Update the choices of the selectInput for "Transcript"
    updateSelectInput(
      session = session,
      inputId = "protein_file",
      label = "Isoform:",
      choices = isoform_choices,
      selected = "Unselected"
    )
    
    # Render the protein viewer
    output$protein_viewer <- renderNGLVieweR({
      protein_file <- input$protein_file
      if (input$protein_file != "Unselected") {
        viewer <- NGLVieweR(protein_file) %>%
          addRepresentation("cartoon")
        
        viewer
      }
    })

    ### MaizeGDB link
    gene_id <- isolate(input$geneID)
    output$maizeGDB <- renderUI({
      link <- paste0("https://www.maizegdb.org/gene_center/gene/", gene_id)
      tags$a(href = link, "MaizeGDB Link")
    })

  })
}

options(shiny.host = if (in_docker) '0.0.0.0' else '127.0.0.1')
options(shiny.port = 8000)

# Run the app
shinyApp(ui, server)
