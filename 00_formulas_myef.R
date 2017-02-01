


if(identical(levels(md$condition2), c("PMN", "MNC")) && 
  identical(levels(md$condition1), c("HD", "Patient"))){
  
  ## create formulas
  formula_lm <- y ~ condition1 + condition2 + condition1:condition2
  formula_lmer <- y ~ condition1 + condition2 + condition1:condition2 + (1|patient_id)
  
  formula_glm_binomial <- cbind(y, total-y) ~ condition1 + condition2 + condition1:condition2
  formula_glm_beta <- y/total ~ condition1 + condition2 + condition1:condition2
  formula_glmer_binomial <- y/total ~ condition1 + condition2 + condition1:condition2 + (1|patient_id)
  
  print(model.matrix(y ~ condition1 + condition2 + condition1:condition2, data = data.frame(y = 1, unique(md[, c("condition1", "condition2")]))))
  
  ## create contrasts
  contrast_names <- c("HDvsP_PMN", "HDvsP_MNC", "HD_PMNvsP_MNC")
  k1 <- c(0, 1, 0, 0)
  k2 <- c(0, 1, 0, 1)
  k3 <- c(0, 1, 1, 1)
  K <- matrix(c(k1, k2, k3), nrow = 3, byrow = TRUE)
  rownames(K) <- contrast_names
  
  ### p-value for sorting the output
  pval_name1 <- "pval_HDvsP_PMN"
  ### p-value for plotting the pheatmap0
  adjpval_name_list <- paste0("adjp_", contrast_names)
  pval_name_list <- paste0("pval_", contrast_names)
  
  ### p-value for plotting the other pheatmaps
  
  results2plot <- list(PMN = list(adjpval_name = "adjp_HDvsP_PMN", pval_name = "pval_HDvsP_PMN", condition1 = c("HD", "Patient"), condition2 = c("PMN")), 
    MNC = list(adjpval_name = "adjp_HDvsP_MNC", pval_name = "pval_HDvsP_MNC", condition1 = c("HD", "Patient"), condition2 = c("MNC")), 
    PMN_MNC = list(adjpval_name = "adjp_HD_PMNvsP_MNC", pval_name = "pval_HD_PMNvsP_MNC", condition1 = c("HD", "Patient"), condition2 = c("PMN", "MNC")))
  
}else{
  
  stop("Metadata does not fit to any of the models that are specified !!!")  
  
}


normalize_expression <- function(expr2norm, md, th){
  
  expr2norm <- expr2norm[, md$shortname]
  
  conditions2 <- levels(md$condition2)

  ### Normalized to mean = 0 and sd = 1 per condition2
  for(i in conditions2){

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





































   
































