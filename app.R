## Author: Jason Yeung
## jy357@bu.edu
## BF591 Final Project

library(shiny)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(comprehenr)
library(DT)
library(stringr)
library(igraph)
library(tidyverse)

ui <- fluidPage(
  titlePanel("BF591 Final Project"),
  h5("R Shiny App was built on se_soybean_cn_sub dataset provided by:"),
  h5("Brown AV, Hudson KA (2015) Developmental profiling of gene expression in soybean trifoliate leaves and cotyledons. BMC Plant Biol 15:169"),
  h5("Datasets and reference can be accessed at:"),
  tags$a(href="https://rdrr.io/bioc/bigPint/man/se_soybean_cn_sub.html", "https://rdrr.io/bioc/bigPint/man/se_soybean_cn_sub.html"), br(),
  br(),
  mainPanel(
    tabsetPanel(type = "tabs",
                tabPanel("Samples",
                         sidebarLayout(
                           sidebarPanel(br(),
                             fileInput("sample_info", "Sample Information (.csv)"),
                             ),
                           mainPanel(br(),
                             tabsetPanel(
                               tabPanel("Summary", br(), tableOutput("sample_summary")),
                               tabPanel("Table", br(), tableOutput("sample_table")),
                               tabPanel("Plot", br(), plotOutput("sample_plot"))
                             )
                           )
                         )
                         ),
                tabPanel("Counts",
                         sidebarLayout(
                           sidebarPanel(br(),
                             fileInput("counts_matrix", "Normalized Counts (.csv)"),
                             sliderInput("variance_slider", "Genes with at Least X Percentile Variance", min=0.00, max=1.00, value=0),
                             sliderInput("nonzero_slider", "Genes with at Least X Non-zero Samples", min=0, max=9, value=0, ticks = FALSE)
                           ),
                           mainPanel(br(),
                             tabsetPanel(
                               tabPanel("Summary", br(), tableOutput("counts_summary")),
                               tabPanel("Diagnostic Plots", br(), plotOutput("counts_vs_variance"), plotOutput("counts_vs_zeros")),
                               tabPanel("Heatmap", br(), plotOutput("counts_heatmap", width = "600px", height = "600px")),
                               tabPanel("PCA", br(), inputPanel(selectInput("xpc_select", label="X-Axis PC",
                                                                            choices = to_list(for(x in 1:9) paste0("PC", x)), 
                                                                            selected = "PC1"),
                                                                selectInput("ypc_select", label="Y-Axis PC",
                                                                            choices = to_list(for(x in 1:9) paste0("PC", x)), 
                                                                            selected = "PC2")
                                                                ), br(),
                                        plotOutput("counts_pca")
                                        )
                               )
                             )
                           )
                         ),
                tabPanel("DE",
                         sidebarLayout(
                           sidebarPanel(br(),
                             fileInput("de_matrix", "DE Analysis Results (.csv)"),
                           ),
                           mainPanel(br(),
                             tabsetPanel(
                               tabPanel("Table", br(), DT::dataTableOutput("de_table", width="800px")),
                               tabPanel("Volcano Plot", br(), inputPanel(selectInput("de_select", label="Select DE Comparison",
                                                                                     choices = list("S1_S2", "S1_S3", "S2_S3"),
                                                                                     selected = "S1_S2"), br(),
                                                                         sliderInput("FDR_slider", "FDR Threshold:",
                                                                                     min = -3,
                                                                                     max = 0,
                                                                                     value = -0.1,
                                                                                     step = 0.1,
                                                                                     ticks = FALSE)
                                                                         ), br(),
                                        plotOutput("de_volcano"), br(),
                                        DT::dataTableOutput("de_comp_table", width="800px"), br()
                                        )
                               )
                             )
                           )
                         ),
                tabPanel("Correlation",
                         sidebarLayout(
                           sidebarPanel(br(),
                             fileInput("counts_network", "Normalized Counts (.csv)"),
                             textAreaInput("gene_list", label = "Enter one gene name per line for correlation analysis:"),
                             inputPanel(HTML("<strong>Example:</strong><br/><br/>Glyma20g32770.1<br/>Glyma03g29150.1<br/>Glyma18g25845.1<br/>Glyma02g45690.2<br/>Glyma08g21390.2")), br(),
                             sliderInput("min_cor", "Minimum Correlation:", min = 0, max = 1, value = 0.5, step = 0.01, ticks = FALSE), br(),
                             actionButton("submit_genes", "Submit")
                           ),
                           mainPanel(br(),
                             tabsetPanel(
                               tabPanel("Heatmap", br(), plotOutput("network_heatmap", width = "600px", height = "600px")),
                               tabPanel("Network", br(), plotOutput("network_graph")),
                               tabPanel("Metrics", br(), 
                                        DT::dataTableOutput("network_table", width="800px"), br()
                                        )
                               )
                             )
                           )
                         )
                )
    )
)

server <- function(input, output, session) {
  # SAMPLE TAB
  load_info <- reactive({
    req(input$sample_info)
    sample_info <- input$sample_info
    if (is.null(sample_info)) {
      return(NULL)
    }
    df <- read.csv(sample_info$datapath, header=TRUE)
    return(df)
  })
  
  draw_sample_summary <- function(dataframe) {
    if (is.null(dataframe)) {
      return(NULL)
    }
    sample_summary <- data.frame(sapply(dataframe, typeof)) %>%
      cbind(Column = row.names(.), .) %>%
      rename_with(.cols = 2, ~"Type")
    return(sample_summary)
  }
  
  draw_sample_bar <- function(dataframe) {
    dataframe["stage"] %>%
      group_by(stage) %>%
      ggplot(aes(x=stage, fill=stage)) + geom_bar() + ylim(0, 5) + theme_bw() + theme(legend.position = "none") + labs(title="Number of Replicates") %>%
      return()
  }
    
  output$sample_summary <- renderTable(draw_sample_summary(dataframe=load_info()), rownames = FALSE, colnames = TRUE)
  output$sample_table <- renderTable(load_info(), rownames = FALSE, colnames = TRUE)
  output$sample_plot <- renderPlot(draw_sample_bar(load_info()))
  
  # COUNTS TAB
  load_counts <- reactive({
    req(input$counts_matrix)
    counts_matrix <- input$counts_matrix
    if (is.null(counts_matrix)) {
      return(NULL)
    }
    df <- read.csv(counts_matrix$datapath, header=TRUE)
    return(df)
  })
  
  calc_variance_genelist <- function(matrix, variance_slider) {
    variance_df <- data.frame(apply(matrix, 1, var)) %>%
      rename_with(.cols=1, ~"Variance") %>%
      filter(Variance >= quantile(Variance, variance_slider))
    genelist <- as.vector(rownames(variance_df))
    return(genelist) # returns list of genes above variance slider
  }
  
  calc_nonzero_genelist <- function(matrix, nonzero_slider) {
    nonzero_df <- matrix[rowSums(matrix == 0) <= nonzero_slider, ]
    genelist <- as.vector(rownames(nonzero_df))
    return(genelist) # returns list of genes above nonzero slider
  }
  
  filter_counts_matrix <- function(matrix, variance_slider, nonzero_slider) {
    filtered_matrix <- matrix %>%
      filter(rownames(.) %in% calc_variance_genelist(matrix, variance_slider)) %>%
      filter(rownames(.) %in% calc_nonzero_genelist(matrix, nonzero_slider)) %>%
      return()  # returns matrix filtered by variance and nonzero gene list
  }
  
  filter_summary <- function(matrix, variance_slider, nonzero_slider) {
    filtered_matrix <- filter_counts_matrix(matrix, variance_slider, nonzero_slider)
    unfiltered_matrix <- matrix
    statistic <- c("# of Samples", "Total # of Genes", "# of Genes Passed", "# of Genes Failed", "% Pass", "% Fail")
    
    num_samples <- ncol(matrix)
    total_num_genes <- nrow(matrix)
    num_pass <- nrow(filtered_matrix)
    num_fail <- nrow(matrix) - nrow(filtered_matrix)
    perc_pass <- nrow(filtered_matrix) / nrow(matrix)
    perc_fail <- (nrow(matrix) - nrow(filtered_matrix)) / nrow(matrix)
    
    value <- c(num_samples, total_num_genes, num_pass, num_fail, perc_pass, perc_fail)
    
    df <- data.frame(statistic, value)
    return(df)
  }
  
  counts_vs_variance_plot <- function(matrix, variance_slider) {
    variance_df <- data.frame(apply(matrix, 1, var)) %>%
      rename_with(.cols=1, ~"Variance")
    median_df <- data.frame(apply(matrix, 1, median)) %>%
      rename_with(.cols=1, ~"Median Count")
    plot <- merge(variance_df, median_df, by="row.names") %>%
      ggplot(aes(x=`Median Count`,
                 y=Variance,
                 color = Variance > quantile(Variance, variance_slider))) +
      geom_point() +
      scale_color_manual(values=c("grey60", "orangered4")) +
      theme_bw() +
      theme(legend.position="none") +
      labs(title="Median Counts vs Variance")
    return(plot)
  }
  
  counts_vs_zeros_plot <- function(matrix, nonzero_slider) {
    nonzero_df <- as.data.frame(matrix[rowSums(matrix == 0) <= 4, ]) 
    nonzero_df$nonzero_tally <- rowSums(nonzero_df != 0)
    median_df <- data.frame(apply(matrix, 1, median)) %>%
      rename_with(.cols=1, ~"Median Count")
    plot <- merge(median_df, nonzero_df, by="row.names") %>%
      ggplot(aes(x=`Median Count`,
                 y=nonzero_tally,
                 color = nonzero_tally >= nonzero_slider)) +
      geom_point() +
      scale_color_manual(values=c("orangered4", "grey60")) +
      theme_bw() +
      theme(legend.position="none") +
      labs(title="Median Counts vs # of Zeros",
           y="# of Zeros")
    return(plot)
  }
  
  counts_heatmap_plot <- function(matrix, variance_slider, nonzero_slider) {
    filtered_matrix <- filter_counts_matrix(matrix, variance_slider, nonzero_slider)
    plot <- heatmap(t(as.matrix(filtered_matrix)), Rowv = FALSE, scale = "row", col = colorRampPalette(brewer.pal(8, "Blues"))(5))
    legend(x="topleft", inset=c(-0.01, -0.01), bty = "n", legend=c("low", "", "", "", "high"), fill=colorRampPalette(brewer.pal(8, "Blues"))(5))
    return(plot)
  }
  
  counts_pca_plot <- function(matrix, variance_slider, nonzero_slider, xpc_select, ypc_select) {
    filtered_matrix <- filter_counts_matrix(matrix, variance_slider, nonzero_slider)
    pca <- prcomp(filtered_matrix, center=TRUE, scale. = TRUE)
    var_explained <- pca$sdev^2/sum(pca$sdev^2)
    plot <- pca$x %>%
      as.data.frame %>%
      ggplot(aes(x=!!sym(xpc_select), y=!!sym(ypc_select))) + geom_point(aes(color="orangered4"), size=2) +
      theme_bw() +
      labs(x=paste0(xpc_select, ": ", round(var_explained[as.integer(substr(xpc_select, 3, 3))]*100, 1), "%"),
           y=paste0(ypc_select, ": ", round(var_explained[as.integer(substr(ypc_select, 3, 3))]*100, 1), "%")) +
      theme(legend.position = "none")
    return(plot)
  }
    
  output$counts_summary <- renderTable(filter_summary(load_counts(), input$variance_slider, input$nonzero_slider))
  output$counts_vs_variance <- renderPlot(counts_vs_variance_plot(load_counts(), input$variance_slider))
  output$counts_vs_zeros <- renderPlot(counts_vs_zeros_plot(load_counts(), input$nonzero_slider))
  output$counts_heatmap <- renderPlot(counts_heatmap_plot(load_counts(), input$variance_slider, input$nonzero_slider))
  output$counts_pca <- renderPlot(counts_pca_plot(load_counts(), input$variance_slider, input$nonzero_slider, input$xpc_select, input$ypc_select))
  
  # DIFFERENTIAL EXPRESSION TAB
  load_de <- reactive({
    req(input$de_matrix)
    de_matrix <- input$de_matrix
    if (is.null(de_matrix)) {
      return(NULL)
    }
    df <- read.csv(de_matrix$datapath, header=TRUE)
    rownames(df)=df[,1]
    df=df[,-1]
    return(df)
  })
  
  filter_de_by_comp <- function(matrix, selection, FDR_threshold) {
    filtered_matrix <- matrix[, grepl(selection, names(matrix))]
    FDR_column <- paste0(selection, ".FDR")
    filtered_matrix %>%
      filter(!is.na(!!sym(FDR_column))) %>%
      filter(!!sym(FDR_column) < 1*10^FDR_threshold) %>%
      return()
  }
  
  de_volcano_plot <- function(matrix, selection, FDR_threshold) {
    FDR_column <- paste0(selection, ".FDR")
    logFC_column <- paste0(selection, ".logFC")
    
    plot <- matrix %>%
      filter(!is.na(!!sym(FDR_column))) %>%
      ggplot(aes(x=!!sym(logFC_column),
                 y=-log10(!!sym(FDR_column)),
                 color=!!sym(FDR_column) < 1*10^FDR_threshold)) +
      geom_point(size=2) +
      scale_color_manual(values=c("grey60", "orangered4")) +
      geom_hline(yintercept = -1*FDR_threshold, linetype='dotted', size=0.8, col='blue') +
      theme_bw() +
      theme(legend.position = "none")
    return(plot)
  }
  
  output$de_table <- DT::renderDataTable({DT::datatable(load_de(), options = list(scrollX = TRUE, pageLength = 10))})
  output$de_comp_table <- DT::renderDataTable({DT::datatable(filter_de_by_comp(load_de(), input$de_select, input$FDR_slider), options = list(pageLength = 10))})
  output$de_volcano <- renderPlot(de_volcano_plot(load_de(), input$de_select, input$FDR_slider))
  
  # CORRELATION TAB
  
  load_counts_network <- eventReactive(input$submit_genes, {
    req(input$counts_network)
    req(input$gene_list)
    
    matrix <- input$counts_network
    if (is.null(matrix)) {
      return(NULL)
    }
    df <- read.csv(matrix$datapath, header=TRUE)
    return(df)
  })
  
  load_gene_list <- eventReactive(input$submit_genes, {
    req(input$gene_list)
    
    gene_list <- input$gene_list
    gene_list <- unlist(strsplit(str_squish(gene_list), split = ' '))
    return(gene_list)
  })
  
  load_min_cor <- eventReactive(input$submit_genes, {
    min_cor <- input$min_cor
    return(min_cor)
  })
  
  filter_matrix_by_gene <- function(matrix, gene_list) {
    filtered_matrix <- matrix[gene_list, ]
    return(filtered_matrix)
  }
  
  network_heatmap_plot <- function(matrix, gene_list) {
    filtered_matrix <- filter_matrix_by_gene(matrix, gene_list)
    plot <- heatmap(t(as.matrix(filtered_matrix)), Rowv = FALSE, scale = "row", margins = c(16,4), col = colorRampPalette(brewer.pal(8, "Blues"))(5))
    legend(x="topleft", inset=c(-0.01, -0.01), bty = "n", legend=c("low", "", "", "", "high"), fill=colorRampPalette(brewer.pal(8, "Blues"))(5))
    return(plot)
  }
  
  generate_cormatrix <- function(matrix, gene_list) {
    filtered_matrix <- filter_matrix_by_gene(matrix, gene_list)
    cormatrix <- cor(t(filtered_matrix)) # transpose so we correlate genes not samples
    diag(cormatrix) <- 0 # set middle diag correlation to 0
    return(cormatrix)
  }
  
  network_graph <- function(cormatrix, min_cor) {
    graph <- graph.adjacency(cormatrix, weighted=TRUE, mode="undirected")
    graph <- delete.edges(graph, E(graph)[abs(weight) < min_cor])
    names <- vertex_attr(graph)$name
    vertex_attr(graph) <- list(name=names, color = rep("slategray1", gorder(graph)))
    return(graph)
  }
  
  plot_network <- function(cormatrix, graph) {
    plot <- plot.igraph(graph, vertex.label.cex=1, vertex.label.dist=2, layout=layout_with_graphopt,
                 edge.width=abs(E(graph)$weight)*3, 
                 edge.color=ifelse(cormatrix > 0, "blue", "red"))
    legend(x="topleft", bty = "n", legend=c("positive", "negative"), fill=c("blue", "red"))
    return(plot)
  }
  
  network_table <- function(graph) {
    old_weights <- E(graph)$weight
    # min-max normalization to handle negative values
    new_weights <- (old_weights - min(old_weights)) / (max(old_weights) - min(old_weights))
    # replace any 0s with 1*10^-9
    new_weights[new_weights == 0] <- 1*10^-9
    # calculate degree, closeness, betweenness
    gene <- V(graph)$name
    degree <- degree(graph)
    closeness <- closeness(graph, weights = new_weights)
    betweenness <- betweenness(graph, weights = new_weights)
    df <- data.frame(gene, degree, closeness, betweenness)
    df <- df %>% remove_rownames %>% column_to_rownames(var="gene") %>% as.matrix
    return(df)
  }
  
  output$network_heatmap <- renderPlot(network_heatmap_plot(load_counts_network(), load_gene_list()))
  output$network_graph <- renderPlot(plot_network(generate_cormatrix(load_counts_network(), load_gene_list()), 
                                                  network_graph(
                                                    generate_cormatrix(load_counts_network(), 
                                                                       load_gene_list()), 
                                                    load_min_cor())), 
                                     width=800, 
                                     height=800)
  output$network_table <- DT::renderDataTable({DT::datatable(network_table(network_graph(generate_cormatrix(load_counts_network(),
                                                                                                            load_gene_list()),
                                                                                         load_min_cor())))})
}


shinyApp(ui=ui, server=server)