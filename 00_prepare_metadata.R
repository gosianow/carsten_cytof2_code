

library(limma)

prepare_metadata <- function(md){
  
  # add more info about samples
  cond_split <- strsplit2(md$condition, "_")
  colnames(cond_split) <- c("day", "response")
  
  md[, c("day", "response")] <- cond_split
  md$response <- factor(md$response, levels = c("NR", "R", "HD"))
  md$response <- factor(md$response)
  md$day <- factor(md$day, levels = c("base", "tx"))
  md$day <- factor(md$day)
  md$patient_id <- factor(md$patient_id)
  
  
  md$condition1 <- md$response
  md$condition2 <- md$day
  
  md$condition <- interaction(md$condition2, md$condition1, lex.order = FALSE, sep = "_")

  md <- md[order(md$condition), , drop = FALSE]


  return(md)
  
}

























