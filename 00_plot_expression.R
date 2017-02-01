


plot_expression <- function(ggdf, ggds, color_groups, outdir, prefix, prefix2){
  
  
  clusters <- levels(ggdf$cluster)
  
  # ------------------------------------
  ## plot each cluster as a separate page in the pdf file

  ggp <- list()

  for(i in 1:nlevels(ggdf$cluster)){
    # i = 1
    
    df <- ggdf[ggdf$cluster == clusters[i], , drop = FALSE]
    
    ggp[[i]] <- ggplot(df, aes(x = group, y = expr, color = group)) +
      geom_boxplot(outlier.colour = NA) +
      geom_point(size = 2, shape = 16, alpha = 0.8, position = position_jitterdodge(jitter.width = 1, jitter.height = 0, dodge.width = 1)) +
      facet_wrap(~ marker, scales = "free") +
      ggtitle(clusters[i]) +
      theme_bw() +
      ylab("Expression") +
      xlab("") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=12, face="bold"), 
        axis.title.y = element_text(size=12, face="bold"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line.x = element_line(size = 0.5, linetype = "solid", color = "black"), 
        axis.line.y = element_line(size = 0.5, linetype = "solid", color = "black"),
        legend.title = element_blank(), legend.position = "right", legend.key = element_blank()) +
      scale_color_manual(values = color_groups)
    
  }

  pdf(file.path(outdir, paste0(prefix, "expr_", prefix2, ".pdf")), w = 18, h = 12, onefile=TRUE)
  for(i in seq(length(ggp)))
    print(ggp[[i]])
  dev.off()

  
  return(NULL)
  
}



