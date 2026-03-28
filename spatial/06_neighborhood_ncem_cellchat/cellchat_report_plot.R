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
  
  # 参数检查
  if (!"pathways.show" %in% names(params)) stop("pathways.show must be specified")
  
  # 创建输出目录
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # 生成动态文件名和标题
  pathway_label <- gsub("[^[:alnum:]]", "_", paste(params$pathways.show, collapse = "+"))
  output_file <- file.path(output_dir, paste0("CellChat_", pathway_label, "_Report.pdf"))
  title <- paste("Signaling Analysis:", paste(params$pathways.show, collapse = " & "))

  # 设置默认值
  if (!"sources.use" %in% names(params)) params$sources.use <- NULL
  if (!"targets.use" %in% names(params)) params$targets.use <- NULL
  if (!"top" %in% names(params)) params$top <- 0.1
  
  # 加载必要的库
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    message("Installing gridExtra package for better plot layout...")
    install.packages("gridExtra")
    library(gridExtra)
  } else {
    library(gridExtra)
  }
  
  # 创建PDF设备
  pdf(output_file, width = fig.width, height = fig.height, onefile = TRUE)
  
  # 添加标题页
  plot.new()
  text(0.5, 0.5, title, cex = 2, font = 2)
  text(0.5, 0.4, paste("Signaling Pathway:", paste(params$pathways.show, collapse = ", ")), 
       cex = 1.5, font = 2)
  text(0.5, 0.3, date(), cex = 1)
  
  # 使用tryCatch来确保即使某个图形失败，报告仍会继续生成
  
  # 1. 圆形图
  tryCatch({
    # 不需要plot.new()，直接绘制
    netVisual_aggregate(cellchat, 
                       signaling = params$pathways.show, 
                       layout = "circle",  
                       top = params$top, 
                       remove.isolate = FALSE,
                       title = "Circle Plot of Cell-Cell Communication")
  }, error = function(e) {
    # 如果出错，创建一个错误信息图
    plot.new()
    text(0.5, 0.5, "Error generating Circle Plot", col = "red")
    text(0.5, 0.4, as.character(e), col = "red", cex = 0.8)
  })

  # 2. 弦图
  tryCatch({
    # 弦图通常返回一个绘图对象
    p <- netVisual_aggregate(cellchat, 
                           signaling = params$pathways.show, 
                           layout = "chord",
                           title = "Chord Diagram of Cell-Cell Communication")
    # 确保绘图对象被正确显示
    if (inherits(p, "ggplot")) {
      print(p)
    }
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, "Error generating Chord Diagram", col = "red")
    text(0.5, 0.4, as.character(e), col = "red", cex = 0.8)
  })

  # 3. 配体-受体对贡献
  tryCatch({
    p <- netAnalysis_contribution(cellchat, 
                                signaling = params$pathways.show,
                                title = "Ligand-Receptor Pair Contribution")
    # 确保绘图对象被正确显示
    if (inherits(p, "ggplot")) {
      print(p)
    }
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, "Error generating Contribution Plot", col = "red")
    text(0.5, 0.4, as.character(e), col = "red", cex = 0.8)
  })

  # 4. 气泡图 (如果有指定来源/目标)
  if (!is.null(params$sources.use) || !is.null(params$targets.use)) {
    tryCatch({
      p <- netVisual_bubble(cellchat, 
                          sources.use = params$sources.use,
                          targets.use = params$targets.use,
                          signaling = params$pathways.show, 
                          remove.isolate = FALSE,
                          title = "Bubble Plot of Specific Communication")
      # 确保绘图对象被正确显示
      if (inherits(p, "ggplot")) {
        print(p)
      }
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, "Error generating Bubble Plot", col = "red")
      text(0.5, 0.4, as.character(e), col = "red", cex = 0.8)
    })
  }
  
  # 5. 信号基因表达
  tryCatch({
    p <- plotGeneExpression(cellchat, 
                          signaling = params$pathways.show, 
                          enriched.only = TRUE, 
                          type = "violin")
    
    # 添加标题并显示
    if (inherits(p, "ggplot")) {
      p <- p + ggtitle("Signaling Gene Expression")
      print(p)
    } else if (is.list(p) && length(p) > 0) {
      # 如果是图形列表，使用grid.arrange显示
      title_grob <- grid::textGrob("Signaling Gene Expression", gp = grid::gpar(fontsize = 14, fontface = "bold"))
      gridExtra::grid.arrange(title_grob, do.call(gridExtra::arrangeGrob, c(p, ncol = 2)), 
                             heights = grid::unit.c(grid::unit(0.05, "npc"), grid::unit(0.95, "npc")))
    } else {
      # 如果不是ggplot对象，也不是列表，尝试其他方式显示
      plot.new()
      text(0.5, 0.9, "Signaling Gene Expression", font = 2, cex = 1.2)
      text(0.5, 0.5, "No visualization available", col = "darkgray")
    }
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, "Error generating Gene Expression Plot", col = "red")
    text(0.5, 0.4, as.character(e), col = "red", cex = 0.8)
  })

  # 6. 网络中心性评分
  tryCatch({
    plot.new()
    text(0.5, 0.97, "Network Centrality Scores", font = 2, cex = 1.5)
    
    # 在当前页面上绘制网络中心性评分图
    # 使用par设置更小的边距，以便有更多空间显示图形
    par(mar = c(4, 3, 0.5, 1) + 0.1)  # 减小所有边距
    
    # 直接调用函数绘制图形，不捕获返回值
    # 增加宽度和高度，增大字体大小
    netAnalysis_signalingRole_network(cellchat, 
                                     signaling = params$pathways.show, 
                                     width = fig.width * 1.5,  # 增加宽度
                                     height = fig.height * 0.6, # 增加高度
                                     font.size = 12)          # 增大字体
    
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, "Error generating Network Centrality Plot", col = "red")
    text(0.5, 0.4, as.character(e), col = "red", cex = 0.8)
  })
  
  # 关闭PDF设备
  dev.off()
  
  message("Report saved to: ", normalizePath(output_file))

  # 额外单独保存热图（如果需要）
  tryCatch({
    # 生成热图文件名
    heatmap_file <- file.path(output_dir, paste0("Heatmap_", pathway_label, ".pdf"))
    
    # 生成热图
    p <- netVisual_heatmap(cellchat, signaling = params$pathways.show, color.heatmap = "Reds")
    
    # 保存热图
    pdf(heatmap_file, width = 7, height = 7.5)
    print(p)
    dev.off()
    
    message("Heatmap also saved separately to: ", normalizePath(heatmap_file))
  }, error = function(e) {
    message("Failed to save separate heatmap: ", as.character(e))
  })
}