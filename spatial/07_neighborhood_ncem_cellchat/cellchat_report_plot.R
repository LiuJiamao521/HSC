#' Generate comprehensive CellChat visualization report
#'
#' @param cellchat CellChat object
#' @param params List containing visualization parameters:
#'        - pathways.show: character vector of signaling pathways to visualize
#'        - sources.use: numeric vector of source cell indices (default NULL)
#'        - targets.use: numeric vector of target cell indices (defau plt NULL)
#'        - top: numeric threshold for filtering (default 0.005)
#' @param output_file Path to output PDF file (default "CellChat_Re./port.pdf")
#' @param fig.width Width of figures in inches (default 10)
#' @param fig.height Height of figures in inches (default 7)
#' @param title Report title (default "CellChat Visualization Report")
#'
#' @return None (generates PDF file)
generate_cellchat_report <- function(cellchat, params, 
                                   output_dir = "./plot",
                                   fig.width = 10, fig.height = 7) {
  
  # Parameter validation
  if (!"pathways.show" %in% names(params)) stop("pathways.show must be specified")
  
  # Create output directory
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Generate dynamic file name and title
  pathway_label <- gsub("[^[:alnum:]]", "_", paste(params$pathways.show, collapse = "+"))
  output_file <- file.path(output_dir, paste0("CellChat_", pathway_label, "_Report.pdf"))
  title <- paste("Signaling Analysis:", paste(params$pathways.show, collapse = " & "))

  # Set default values
  if (!"sources.use" %in% names(params)) params$sources.use <- NULL
  if (!"targets.use" %in% names(params)) params$targets.use <- NULL
  if (!"top" %in% names(params)) params$top <- 0.1
  
  # Load required libraries
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    message("Installing gridExtra package for better plot layout...")
    install.packages("gridExtra")
    library(gridExtra)
  } else {
    library(gridExtra)
  }
  
  # Create PDF device
  pdf(output_file, width = fig.width, height = fig.height, onefile = TRUE)
  
  # Add title page
  plot.new()
  text(0.5, 0.5, title, cex = 2, font = 2)
  text(0.5, 0.4, paste("Signaling Pathway:", paste(params$pathways.show, collapse = ", ")), 
       cex = 1.5, font = 2)
  text(0.5, 0.3, date(), cex = 1)
  
  # Use tryCatch so report generation continues even if one plot fails
  
  # 1. Circle plot
  tryCatch({
    # No plot.new() needed; draw directly
    netVisual_aggregate(cellchat, 
                       signaling = params$pathways.show, 
                       layout = "circle",  
                       top = params$top, 
                       remove.isolate = FALSE,
                       title = "Circle Plot of Cell-Cell Communication")
  }, error = function(e) {
    # If an error occurs, draw an error message panel
    plot.new()
    text(0.5, 0.5, "Error generating Circle Plot", col = "red")
    text(0.5, 0.4, as.character(e), col = "red", cex = 0.8)
  })

  # 2. Chord diagram
  tryCatch({
    # Chord plotting usually returns a plot object
    p <- netVisual_aggregate(cellchat, 
                           signaling = params$pathways.show, 
                           layout = "chord",
                           title = "Chord Diagram of Cell-Cell Communication")
    # Ensure plot object is rendered correctly
    if (inherits(p, "ggplot")) {
      print(p)
    }
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, "Error generating Chord Diagram", col = "red")
    text(0.5, 0.4, as.character(e), col = "red", cex = 0.8)
  })

  # 3. Ligand-receptor pair contribution
  tryCatch({
    p <- netAnalysis_contribution(cellchat, 
                                signaling = params$pathways.show,
                                title = "Ligand-Receptor Pair Contribution")
    # Ensure plot object is rendered correctly
    if (inherits(p, "ggplot")) {
      print(p)
    }
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, "Error generating Contribution Plot", col = "red")
    text(0.5, 0.4, as.character(e), col = "red", cex = 0.8)
  })

  # 4. Bubble plot (if source/target are specified)
  if (!is.null(params$sources.use) || !is.null(params$targets.use)) {
    tryCatch({
      p <- netVisual_bubble(cellchat, 
                          sources.use = params$sources.use,
                          targets.use = params$targets.use,
                          signaling = params$pathways.show, 
                          remove.isolate = FALSE,
                          title = "Bubble Plot of Specific Communication")
      # Ensure plot object is rendered correctly
      if (inherits(p, "ggplot")) {
        print(p)
      }
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, "Error generating Bubble Plot", col = "red")
      text(0.5, 0.4, as.character(e), col = "red", cex = 0.8)
    })
  }
  
  # 5. Signaling gene expression
  tryCatch({
    p <- plotGeneExpression(cellchat, 
                          signaling = params$pathways.show, 
                          enriched.only = TRUE, 
                          type = "violin")
    
    # Add title and render
    if (inherits(p, "ggplot")) {
      p <- p + ggtitle("Signaling Gene Expression")
      print(p)
    } else if (is.list(p) && length(p) > 0) {
      # If result is a plot list, display with grid.arrange
      title_grob <- grid::textGrob("Signaling Gene Expression", gp = grid::gpar(fontsize = 14, fontface = "bold"))
      gridExtra::grid.arrange(title_grob, do.call(gridExtra::arrangeGrob, c(p, ncol = 2)), 
                             heights = grid::unit.c(grid::unit(0.05, "npc"), grid::unit(0.95, "npc")))
    } else {
      # If neither ggplot nor list, use fallback display
      plot.new()
      text(0.5, 0.9, "Signaling Gene Expression", font = 2, cex = 1.2)
      text(0.5, 0.5, "No visualization available", col = "darkgray")
    }
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, "Error generating Gene Expression Plot", col = "red")
    text(0.5, 0.4, as.character(e), col = "red", cex = 0.8)
  })

  # 6. Network centrality scores
  tryCatch({
    plot.new()
    text(0.5, 0.97, "Network Centrality Scores", font = 2, cex = 1.5)
    
    # Draw network centrality score plot on the current page
    # Use smaller margins via par() to free more plotting space
    par(mar = c(4, 3, 0.5, 1) + 0.1)  # Reduce all margins
    
    # Draw directly without capturing return value
    # Increase width/height and enlarge font
    netAnalysis_signalingRole_network(cellchat, 
                                     signaling = params$pathways.show, 
                                     width = fig.width * 1.5,  # Increase width
                                     height = fig.height * 0.6, # Increase height
                                     font.size = 12)          # Larger font
    
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, "Error generating Network Centrality Plot", col = "red")
    text(0.5, 0.4, as.character(e), col = "red", cex = 0.8)
  })
  
  # Close PDF device
  dev.off()
  
  message("Report saved to: ", normalizePath(output_file))

  # Optionally save heatmap separately
  tryCatch({
    # Generate heatmap file name
    heatmap_file <- file.path(output_dir, paste0("Heatmap_", pathway_label, ".pdf"))
    
    # Generate heatmap
    p <- netVisual_heatmap(cellchat, signaling = params$pathways.show, color.heatmap = "Reds")
    
    # Save heatmap
    pdf(heatmap_file, width = 7, height = 7.5)
    print(p)
    dev.off()
    
    message("Heatmap also saved separately to: ", normalizePath(heatmap_file))
  }, error = function(e) {
    message("Failed to save separate heatmap: ", as.character(e))
  })
}