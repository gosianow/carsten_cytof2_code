##############################################################################
## <<04_expression.R>>

# BioC 3.3
# Created 22 Aug 2016
# Updated 13 Oct 2016

##############################################################################
Sys.time()
##############################################################################

# Load packages
library(gdata)
library(Repitools)
library(gplots)
library(ggplot2)
library(plyr)
library(reshape2)
library(RColorBrewer)
library(limma)
library(pheatmap)

##############################################################################
# Test arguments
##############################################################################

rwd='/home/Shared/data/cytof/carsten_cytof/CK_2016-06-23_01'
expr_prefix='23_01_pca1_mergingNEW2_raw_'
expr_outdir='080_expression'
path_data='010_data/23_01_expr_raw.rds'
path_metadata='/home/Shared/data/cytof/carsten_cytof/CK_metadata/metadata_23_01.rds'
path_clustering_observables='030_heatmaps/23_01_pca1_clustering_observables.xls'
path_clustering='030_heatmaps/23_01_pca1_mergingNEW2_clustering.xls'
path_clustering_labels='030_heatmaps/23_01_pca1_mergingNEW2_clustering_labels.xls'
path_fun_models='/home/gosia/R/carsten_cytof_code/00_models.R'
path_fun_formulas='/home/gosia/R/carsten_cytof_code/00_formulas_1dataset_3responses.R'
path_fun_plot_heatmaps <- "/home/gosia/R/carsten_cytof_code/00_plot_heatmaps_for_sign_expr.R"
path_fun_plot_expression <- "/home/gosia/R/carsten_cytof_code/00_plot_expression.R"
analysis_type='all'


##############################################################################
# Read in the arguments
##############################################################################

rm(list = ls())

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

cat(paste0(args, collapse = "\n"), fill = TRUE)


##############################################################################

setwd(rwd)

prefix <- expr_prefix
outdir <- expr_outdir
suffix <- ""

if( !file.exists(outdir) ) 
  dir.create(outdir, recursive = TRUE)

if(!analysis_type %in% c("clust", "all"))
  stop("analysis_type must be 'all' or 'clust'!")

out_name <- analysis_type

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

# read expression data
if(grepl(".txt", path_data)){
  expr <- read.table(path_data, header = TRUE, sep = "\t", as.is = TRUE)
}
if(grepl(".rds", path_data)){
  expr <- readRDS(path_data)
}


fcs_colnames <- colnames(expr)[!grepl("cell_id|sample_id", colnames(expr))]
e <- expr[, fcs_colnames]

# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

md <- readRDS(path_metadata)

### Colors 

colors <- unique(md[, c("condition", "color")])

color_groups <- colors$color
names(color_groups) <- colors$condition

color_groupsb <- adjustcolor(color_groups, alpha = 0.3)
names(color_groupsb) <- colors$condition


# ------------------------------------------------------------
# Load clustering data
# ------------------------------------------------------------

# clustering observables
clustering_observables <- read.table(path_clustering_observables, header = TRUE, sep = "\t", as.is = TRUE)
rownames(clustering_observables) <- clustering_observables$mass

clust_observ <- clustering_observables[clustering_observables$clustering_observable, "mass"]


# clustering labels
labels <- read.table(path_clustering_labels, header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))
rownames(labels) <- labels$cluster

# clustering
clustering <- read.table(path_clustering, header = TRUE, sep = "\t", as.is = TRUE)

## drop the "drop" cluster
if("drop" %in% labels$label){
  
  clust2drop <- labels$cluster[labels$label == "drop"]
  cells2drop <- clustering$cluster != clust2drop
  clustering <- clustering[cells2drop, , drop = FALSE]
  labels <- labels[labels$label != "drop", ,drop = FALSE]
  labels$label <- factor(labels$label)
  e <- e[cells2drop, , drop = FALSE]
  
}


clust <- clustering[, "cluster"]

samp <- clustering[, "sample_id"]


# ------------------------------------------------------------
# get the isotope and antigen for fcs markers

m <- match(fcs_colnames, clustering_observables$mass)

fcs_panel <- data.frame(colnames = fcs_colnames, Isotope = clustering_observables$mass[m], Antigen = clustering_observables$marker[m], stringsAsFactors = FALSE)


# ------------------------------------------------------------

# Indeces of observables used for clustering 
scols <- which(fcs_colnames %in% clust_observ)

# Indeces of other observables
xcols <- which(!fcs_colnames %in% clust_observ)


# ordered by decreasing pca score
if("avg_score" %in% colnames(clustering_observables)){
  scols <- scols[order(clustering_observables[fcs_colnames[scols], "avg_score"], decreasing = TRUE)]
  xcols <- xcols[order(clustering_observables[fcs_colnames[xcols], "avg_score"], decreasing = TRUE)]
}


# ------------------------------------------------------------
# Get the median expression per cluster, and if sample has not enough cells, set expression to NA
# ------------------------------------------------------------


if(analysis_type == "clust"){
  
  min_cells <- 20
  
  table_samp <- aggregate(e[, 1, drop = FALSE], by = list(clust, samp), FUN = length, drop = FALSE)
  
  keep_samps <- table_samp[, 3] > min_cells
  
  colnames(e) <- fcs_panel$Antigen
  
  a <- aggregate(e, by = list(clust, samp), FUN = median, drop = FALSE)
  
  mlab <- match(a$Group.1, labels$cluster)
  a$label <- labels$label[mlab]
  
  colnames(a)[1:2] <- c("cluster", "sample")
  
  a <- a[, c("cluster", "label", "sample", fcs_panel$Antigen[c(scols, xcols)])]
  a$label <- factor(a$label, levels = labels$label)
  
  a[!keep_samps, fcs_panel$Antigen[c(scols, xcols)]] <- NA
  
  ### Save the median expression per cluster and sample
  write.table(a, file.path(outdir, paste0(prefix, "expr_clust.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  ### Skipp samples that have NA for all the clusters
  ### Clusters with at least one NA will be skipped in the expression analysis (see 00_models.R)
  
  table_sample <- table(a[complete.cases(a), "sample"])
  
  keep_samps <- names(table_sample)[table_sample > 0]
  
  a <- a[a$sample %in% keep_samps, , drop = FALSE]
  
  md <- md[md$shortname %in% keep_samps, , drop = FALSE]
  
  ## drop unused levels
  md$condition1 <- factor(md$condition1)
  md$condition2 <- factor(md$condition2)
  md$patient_id <- factor(md$patient_id)
  
}

if(analysis_type == "all"){
  
  min_cells <- 50
  
  table_samp <- table(samp)
  
  keep_samps <- names(table_samp)[which(table_samp > min_cells)]
  
  
  colnames(e) <- fcs_panel$Antigen
  
  a <- aggregate(e, by = list(samp), FUN = median)
  
  colnames(a)[1] <- c("sample")
  
  a$cluster <- -1
  a$label <- "all"
  
  a <- a[, c("cluster", "label", "sample", fcs_panel$Antigen[c(scols, xcols)])]
  a$label <- factor(a$label)
  
  a[!a$sample %in% keep_samps, fcs_panel$Antigen[c(scols, xcols)]] <- NA
  
  ### Save the median expression per cluster and sample
  write.table(a, file.path(outdir, paste0(prefix, "expr_all.xls")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  ### Keep only those samples that have enough cells 
  
  a <- a[a$sample %in% keep_samps, , drop = FALSE]
  
  md <- md[md$shortname %in% keep_samps, , drop = FALSE]
  
  ## drop unused levels
  md$condition1 <- factor(md$condition1)
  md$condition2 <- factor(md$condition2)
  md$patient_id <- factor(md$patient_id)
  
  
}



# -----------------------------------------------------------------------------
# Plot expression per cluster
# -----------------------------------------------------------------------------

source(path_fun_plot_expression)

ggdf <- melt(a, id.vars = c("cluster", "label", "sample"), value.name = "expr", variable.name = "marker")

## use labels as clusters
ggdf$cluster <- ggdf$label

## add group info
mm <- match(ggdf$samp, md$shortname)
ggdf$group <- factor(md$condition[mm])

## order markers as for heatmaps
ggdf$marker <- factor(ggdf$marker, levels = fcs_panel$Antigen[c(scols, xcols)])


plot_expression(ggdf = ggdf, color_groups = color_groups, outdir = outdir, prefix = prefix, prefix2 = out_name)



# -----------------------------------------------------------------------------
# Test for marker expression differences between groups overall and per cluster 
# -----------------------------------------------------------------------------
## The model functions do not anlyse a cluster with NAs; 
## For merged data it means such cluster was not present in all the datasets
## For expression data clusters with no cells are skipped


### Load functions fitting models
source(path_fun_models)
### Load formulas that are fit in the models - this function may change the md object!!!
source(path_fun_formulas)

source(path_fun_plot_heatmaps)


levels(md$day)
levels(md$response)


# -----------------------------------------------------------------------------
### normalize the expression
# -----------------------------------------------------------------------------

### Prepare the matrix with data (rows - markers X clusters; columns - samples)

expr <- a
exprm <- melt(expr, id.vars = c("cluster", "label", "sample"), value.name = "expr", variable.name = "marker")
exprc <- dcast(exprm, cluster + label + marker ~ sample, value.var = "expr")


expr2norm <- exprc[md$shortname]
th <- 2.5


expr2norm_out <- normalize_expression(expr2norm = expr2norm, md = md, th = th)


expr_norm <- data.frame(exprc[, c("cluster", "label", "marker")], expr2norm_out)

breaks = seq(from = -th, to = th, length.out = 101)
legend_breaks = seq(from = -round(th), to = round(th), by = 1)




# models2fit <- c("lm_interglht", "lmer_interglht", "rlm_interglht")
models2fit <- c("lmer_interglht")


for(k in models2fit){
  # k = "lm_interglht"
  print(k)
  
  
  switch(k,
    lm_interglht = {
      
      # Fit a LM with interactions + test contrasts with multcomp pckg
      fit_out <- fit_lm_interglht(data = exprc, md, method = "lm", formula = formula_lm, K = K)
      
    }, 
    lmer_interglht = {
      
      # Fit a lmer with interactions + test contrasts with multcomp pckg
      fit_out <- fit_lmer_interglht(data = exprc, md, formula = formula_lmer, K = K)
      
    },
    lmrob_interglht = {
      ## Problems with running lmrob!!!
      fit_out <- fit_lm_interglht(data = exprc, md, method = "lmrob", formula = formula_lm, K = K)
    },
    rlm_interglht = {
      fit_out <- fit_lm_interglht(data = exprc, md, method = "rlm", formula = formula_lm, K = K)
    }
    
  )
  
  
  # ----------------------------------------
  # Extract p-values and coeffs
  # ----------------------------------------
  
  pvs <- data.frame(exprc[, c("cluster", "label", "marker")], fit_out[["pvals"]])
  coeffs <- data.frame(exprc[, c("cluster", "label", "marker")], fit_out[["coeffs"]])
  
  oo <- order(pvs[, pval_name1], decreasing = FALSE)
  pvs <- pvs[oo, , drop = FALSE]
  coeffs <- coeffs[oo, , drop = FALSE]
  
  ## save the results
  write.table(pvs, file = file.path(outdir, paste0(prefix, "expr_", out_name, "_pvs_", k, suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  write.table(coeffs, file = file.path(outdir, paste0(prefix, "expr_", out_name, "_coeffs_", k, suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  
  
  # ----------------------------------------
  # Plot a heatmap of significant cases
  # ----------------------------------------
  
  ### add p-value info
  expr_all <- merge(pvs, expr_norm, by = c("cluster", "label", "marker"), all.x = TRUE, sort = FALSE)
  
  prefix2 <- paste0(out_name, "_", k, "_")
  FDR_cutoff <- 0.05
  
  plot_heatmaps_for_sign_expr(expr_all = expr_all, md = md, FDR_cutoff = FDR_cutoff, pval_name_list = pval_name_list, adjpval_name_list = adjpval_name_list, results2plot = results2plot, breaks = breaks, legend_breaks = legend_breaks, color_groups = color_groups, outdir = outdir, prefix = prefix, prefix2 = prefix2, suffix = suffix)
  
  
}







sessionInfo()














################################
### 04_expression done!
################################