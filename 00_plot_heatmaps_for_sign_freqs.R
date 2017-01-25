


plot_heatmaps_for_sign_freqs <- function(expr_all, md, FDR_cutoff, pval_name_list, adjpval_name_list, results2plot, breaks, legend_breaks, color_groups, outdir, prefix, prefix2, suffix){
  
  
  md <- md[order(md$condition), , drop = FALSE]
  rownames(md) <- md$shortname
  
  samples2plot <- md$shortname[md$shortname %in% colnames(expr_all)]
  
  md <- md[samples2plot, ]
  md$condition <- factor(md$condition)
  md$condition1 <- factor(md$condition1)
  md$condition2 <- factor(md$condition2)
  
  gaps_col <- cumsum(table(md$condition))
  
  
  # -----------------------------
  ### Plot one heatmap with all the conditions and all the p-values
  
  ### order the expression by adjpval
  for(i in length(pval_name_list):1){
    expr_all <- expr_all[order(expr_all[, pval_name_list[i]]), , drop = FALSE]
  }
  
  which_top_pvs <- rep(TRUE, nrow(expr_all))
  
  
  if(sum(which_top_pvs) > 0){
    print("Plot pheatmap0")
    
    expr_heat <- expr_all[which_top_pvs, , drop = FALSE]
    expr <- expr_heat[, samples2plot, drop = FALSE]
    
    gaps_row <- NULL
    labels_row <- paste0(expr_heat$label)
    labels_col <- colnames(expr)
    
    annotation_col <- data.frame(condition = factor(md[colnames(expr), "condition"]))
    rownames(annotation_col) <- colnames(expr)
    
    annotation_colors <- list(condition = color_groups[levels(annotation_col$condition)])
    
    pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_col = annotation_col, annotation_colors = annotation_colors, annotation_legend = FALSE, filename = file.path(outdir, paste0(prefix, "frequencies_", prefix2, "pheatmap0", suffix, ".pdf")))
    
    
    pvs_heat <- expr_heat[, adjpval_name_list, drop = FALSE]
    
    labels_col <- colnames(pvs_heat)
    
    pheatmap(pvs_heat, cellwidth = 60, cellheight = 24, color = c("grey50", "grey70", "grey90"), breaks = c(0, 0.05, 0.1, 1), legend_breaks = c(0, 0.05, 0.1, 1), legend = FALSE, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = NULL, gaps_row = gaps_row, display_numbers = TRUE, number_format = "%.02e", number_color = "black", fontsize_row = 14, fontsize_col = 14, fontsize = 12, filename = file.path(outdir, paste0(prefix, "frequencies_", prefix2, "pheatmap0pvs", suffix, ".pdf")))
    
    
  }
  
  
  # -----------------------------
  ### Plot specified heatmaps
  
  results_names <- names(results2plot)
  
  for(i in 1:length(results2plot)){
    # i = 1
    
    adjpval_name <- results2plot[[i]][["adjpval_name"]]
    pval_name <- results2plot[[i]][["pval_name"]]
    
    samples2plot <- md$shortname[md$condition1 %in% results2plot[[i]][["condition1"]] & md$condition2 %in% results2plot[[i]][["condition2"]]]
    
    
    expr_all <- expr_all[order(expr_all[, pval_name]), , drop = FALSE]
    
    which_top_pvs <- expr_all[, adjpval_name] < FDR_cutoff & !is.na(expr_all[, adjpval_name])
    which(which_top_pvs)
    
    
    if(sum(which_top_pvs) > 0) {
      
      expr_heat <- expr_all[which_top_pvs, , drop = FALSE]
      expr <- expr_heat[, samples2plot, drop = FALSE]
      
      gaps_row <- NULL
      gaps_col <- cumsum(table(factor(md[samples2plot, "condition"])))
      labels_row <- paste0(expr_heat$label, " (", sprintf( "%.02e", expr_heat[, adjpval_name]), ")")
      labels_col <- colnames(expr)
      
      annotation_col <- data.frame(condition = factor(md[colnames(expr), "condition"]))
      rownames(annotation_col) <- colnames(expr)
      
      annotation_colors <- list(condition = color_groups[levels(annotation_col$condition)])
      
      pheatmap(expr, cellwidth = 28, cellheight = 24, color = colorRampPalette(c("#87CEFA", "#56B4E9", "#0072B2", "#000000", "#D55E00", "#E69F00", "#FFD700"), space = "Lab")(100), breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, gaps_row = gaps_row, fontsize_row = 14, fontsize_col = 14, fontsize = 12, annotation_col = annotation_col, annotation_colors = annotation_colors, annotation_legend = FALSE, filename = file.path(outdir, paste0(prefix, "frequencies_", prefix2, "pheatmap", i, "_", results_names[i], suffix, ".pdf")))
      
      
    }else{
      
      pdf(file.path(outdir, paste0(prefix, "frequencies_", prefix2, "pheatmap", i, "_", results_names[i],  suffix, ".pdf")))
      plot(1, type="n", axes=F, xlab="", ylab="")
      dev.off()
      
    }
    
    
    
  }
  
  
  return(NULL)
  
}



