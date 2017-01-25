

prepare_metadata <- function(md){
  
  md$condition1 <- factor(md$condition1)
  md$condition2 <- factor(md$condition2)
  
  md$condition <- interaction(md$condition2, md$condition1, lex.order = TRUE, sep = "_")

  md <- md[order(md$condition), , drop = FALSE]
  
  return(md)
  
}

























