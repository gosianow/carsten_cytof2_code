


if(identical(levels(md$condition2), c("base", "tx")) && identical(levels(md$condition1), c("NR", "R", "HD"))){
  ## create formulas
  formula_lm <- y ~ condition1 + condition2 + condition1:condition2
  formula_lmer <- y ~ condition1 + condition2 + condition1:condition2 + (1|patient_id)
  
  formula_glm_binomial <- cbind(y, total-y) ~ condition1 + condition2 + condition1:condition2
  formula_glm_beta <- y/total ~ condition1 + condition2 + condition1:condition2
  formula_glmer_binomial <- y/total ~ condition1 + condition2 + condition1:condition2 + (1|patient_id)
  
  ## create contrasts
  contrast_names <- c("NRvsR", "NRvsR_base", "NRvsR_tx", "NRvsR_basevstx")
  k0 <- c(0, 1, 0, 0, 1/2, 0) # NR vs R
  k1 <- c(0, 1, 0, 0, 0, 0) # NR vs R in base
  k2 <- c(0, 1, 0, 0, 1, 0) # NR vs R in tx
  k3 <- c(0, 0, 0, 0, 1, 0) # whether NR vs R is different in base and tx
  K <- matrix(c(k0, k1, k2, k3), nrow = 4, byrow = TRUE)
  rownames(K) <- contrast_names
  
  ### p-value for sorting the output
  pval_name1 <- "pval_NRvsR"
  ### p-value for plotting the pheatmap0
  adjpval_name_list <- c("adjp_NRvsR", "adjp_NRvsR_base", "adjp_NRvsR_tx", "adjp_NRvsR_basevstx")
  pval_name_list <- c("pval_NRvsR", "pval_NRvsR_base", "pval_NRvsR_tx", "pval_NRvsR_basevstx")
  
  ### p-value for plotting the other pheatmaps
  
  results2plot <- list(both = list(adjpval_name = "adjp_NRvsR", pval_name = "pval_NRvsR", condition1 = c("NR", "R"), condition2 = c("base", "tx")), base = list(adjpval_name = "adjp_NRvsR_base", pval_name = "pval_NRvsR_base", condition1 = c("NR", "R"), condition2 = c("base")), tx = list(adjpval_name = "adjp_NRvsR_tx", pval_name = "pval_NRvsR_tx", condition1 = c("NR", "R"), condition2 = c("tx")))
  
  plot_condition1_levels <- c("NR", "R")
  
}else{
  
  stop("Metadata does not fit to any the models that are specified !!!")  
  
}


normalize_expression <- function(expr2norm, md, th){
  
  md <- md[md$condition1 %in% plot_condition1_levels, , drop = FALSE]
  
  expr2norm <- expr2norm[, md$shortname]
  
  conditions2 <- levels(md$condition2)

  ### Normalized to mean = 0 and sd = 1 per day
  for(i in conditions2){
    # i = "base"
    expr2norm[, md$condition2 == i] <- t(apply(expr2norm[, md$condition2 == i, drop = FALSE], 1, function(x){ 
      
      if(sum(!is.na(x)) == 0)
        return(x)
      
      if(sum(!is.na(x)) < 2)
        return(x-mean(x, na.rm = TRUE))
      
      sdx <- sd(x, na.rm = TRUE)
      if(sdx == 0)
        x <- (x-mean(x, na.rm = TRUE))
      else 
        x <- (x-mean(x, na.rm = TRUE))/sdx
      
      x[x > th] <- th
      x[x < -th] <- -th
      
      return(x)}))
  }

  return(expr2norm)
  
}




































