##############################################################################
## <<04_frequencies.R>>

# BioC 3.3
# Created 27 July 2016
# Updated 1 Feb 2017

##############################################################################
Sys.time()
##############################################################################

# Load packages
library(gdata)
library(gplots)
library(ggplot2)
library(plyr)
library(reshape2)
library(limma) # for strsplit2
library(RColorBrewer)
library(pheatmap)
library(gtools) # for logit


##############################################################################
# Test arguments
##############################################################################

rwd='/home/Shared/data/cytof/carsten_cytof2/MyeEUNITERfinal_neutrophils_merging'
freq_prefix='myefNEUTROP_myef_pca0_cl10_'
freq_outdir='050_frequencies_auto//'
path_metadata='/home/Shared/data/cytof/carsten_cytof2/MyeEUNITERfinal_metadata/metadata_MyeEUNITERfinal.rds'
path_clustering='030_heatmaps/myefNEUTROP_myef_pca0_cl10_clustering.xls'
path_clustering_labels='030_heatmaps/myefNEUTROP_myef_pca0_cl10_clustering_labels.xls'
path_fun_models='/home/gosia/R/carsten_cytof2_code/00_models.R'
path_fun_formulas='/home/gosia/R/carsten_cytof2_code/00_formulas_myef.R'
path_fun_plot_heatmaps='/home/gosia/R/carsten_cytof2_code/00_plot_heatmaps_for_sign_freqs.R'
path_fun_plot_frequencies='/home/gosia/R/carsten_cytof2_code/00_plot_frequencies.R'

### Optional arguments
pdf_hight=4

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

prefix <- freq_prefix
suffix <- ""
outdir <- freq_outdir

if(!file.exists(outdir)) 
  dir.create(outdir, recursive = TRUE)



if(!any(grepl("pdf_hight=", args))){
  pdf_hight=4
}


# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

md <- readRDS(path_metadata)
md


### Colors 

colors <- unique(md[, c("condition", "color")])

color_groups <- colors$color
names(color_groups) <- colors$condition

color_groupsb <- adjustcolor(color_groups, alpha = 0.3)
names(color_groupsb) <- colors$condition

# ------------------------------------------------------------
# Load cluster data
# ------------------------------------------------------------

## clustering labels
labels <- read.table(path_clustering_labels, header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))

## clustering
clustering <- read.table(path_clustering, header = TRUE, sep = "\t", as.is = TRUE)

## drop the "drop" cluster
if("drop" %in% labels$label){
  
  clust2drop <- labels$cluster[labels$label == "drop"]
  cells2drop <- clustering$cluster != clust2drop
  clustering <- clustering[cells2drop, , drop = FALSE]
  labels <- labels[labels$label != "drop", ,drop = FALSE]
  labels$label <- factor(labels$label)
  
}

clust <- clustering[, "cluster"]

# ---------------------------------------
# Calculate the cluster frequencies per sample
# ---------------------------------------

samp <- clustering[, "sample_id"]

# calculate frequencies
freq <- table(clust, samp)

prop <- t(t(freq) / colSums(freq)) * 100 # proportion of clusters in samples

# use labels as names of clusters
mlab <- match(rownames(freq), labels$cluster)


### Save the frequencies and proportions
prop_out <- data.frame(labels[mlab, c("cluster", "label")], as.data.frame.matrix(prop))

freq_out <- data.frame(labels[mlab, c("cluster", "label")], as.data.frame.matrix(freq))

write.table(prop_out, file=file.path(outdir, paste0(prefix, "frequencies.xls")), row.names=FALSE, quote=FALSE, sep="\t")
write.table(freq_out, file=file.path(outdir, paste0(prefix, "counts.xls")), row.names=FALSE, quote=FALSE, sep="\t")



# ---------------------------------------
# Keep only those samples that have enough cells 
# ---------------------------------------

min_cells <- 0

table_samp <- colSums(freq_out[md$shortname], na.rm = TRUE)
names(table_samp) <- md$shortname

keep_samps <- names(table_samp)[which(table_samp > min_cells)]

prop_out <- prop_out[, colnames(prop_out) %in% c("cluster", "label", keep_samps), drop = FALSE]
freq_out <- freq_out[, colnames(freq_out) %in% c("cluster", "label", keep_samps), drop = FALSE]

md <- md[md$shortname %in% keep_samps, , drop = FALSE]



## drop unused levels
md$condition1 <- factor(md$condition1)
md$condition2 <- factor(md$condition2)
md$patient_id <- factor(md$patient_id)


# ------------------------------------------------------------
### Plot frequencies
# ------------------------------------------------------------

source(path_fun_plot_frequencies)

ggdf <- melt(prop_out, id.vars = c("cluster", "label"), value.name = "prop", variable.name = "samp")

## use labels as clusters
ggdf$cluster <- factor(ggdf$label, levels = labels$label)
ggdf <- ggdf[, c("cluster", "samp", "prop")]

## add group info
mm <- match(ggdf$samp, md$shortname)
ggdf$group <- factor(md$condition[mm])
ggdf$condition2 <- factor(md$condition2[mm])


plot_frequencies(ggdf = ggdf, color_groups = color_groups, color_groupsb = color_groupsb, outdir = outdir, prefix = prefix, pdf_hight = pdf_hight)



# ------------------------------------------------------------
# Test for frequency differences between groups
# ------------------------------------------------------------
## The model functions do not anlyse a cluster with NAs; 
## For merged data it means such cluster was not present in all the datasets
## For expression data clusters with no cells are skipped


### Load functions fitting models
source(path_fun_models)

### Load formulas that are fit in the models - this function may change the md object!!!
source(path_fun_formulas)

source(path_fun_plot_heatmaps)


levels(md$condition2)
levels(md$condition1)



# ------------------------------------------------------------
### normalize the expression
# ------------------------------------------------------------

ass_freq_out <- freq_out
ass_freq_out[md$shortname] <- asin(sqrt((t(t(freq_out[md$shortname]) / colSums(freq_out[md$shortname])))))

expr2norm <- ass_freq_out[md$shortname]
th <- 2.5


expr2norm_out <- normalize_expression(expr2norm = expr2norm, md = md, th = th)

expr_norm <- data.frame(ass_freq_out[, c("cluster", "label")], expr2norm_out)


breaks = seq(from = -th, to = th, length.out = 101)
legend_breaks = seq(from = -round(th), to = round(th), by = 1)




# models2fit <- c("glm_binomial_interglht", "glm_quasibinomial_interglht", "glmer_binomial_interglht", "lmer_arcsinesqrt_interglht", "lmer_logit_interglht", "lm_arcsinesqrt_interglht", "lm_logit_interglht")
models2fit <- c("glmer_binomial_interglht")

for(k in models2fit){
  # k = "glmer_binomial_interglht"
  print(k)
  
  switch(k,
    glm_binomial_interglht = {
      # Fit a GLM binomial with interactions + test contrasts with multcomp pckg
      fit_out <- fit_glm_interglht(data = freq_out, md, family = "binomial", formula = formula_glm_binomial, K = K)
      
    }, 
    glm_quasibinomial_interglht = {
      # Fit a GLM quasibinomial with interactions + test contrasts with multcomp pckg
      fit_out <- fit_glm_interglht(data = freq_out, md, family = "quasibinomial", formula = formula_glm_binomial, K = K)
      
    },
    glmer_binomial_interglht = {
      # Fit a GLMM binomial with interactions + test contrasts with multcomp pckg
      fit_out <- fit_glmer_interglht(data = freq_out, md, family = "binomial", formula = formula_glmer_binomial, K = K)
      
    },
    lmer_logit_interglht = {
      
      logit_freq_out <- freq_out
      logit_freq_out[md$shortname] <- logit(t(t(freq_out[md$shortname]) / colSums(freq_out[md$shortname])))
      ## Be carefull about Inf and -Inf for prop = 0, 1
      fit_out <- fit_lmer_interglht(data = logit_freq_out, md, formula = formula_lmer, K = K)
      
    },
    lm_logit_interglht = {
      
      logit_freq_out <- freq_out
      logit_freq_out[md$shortname] <- logit(t(t(freq_out[md$shortname]) / colSums(freq_out[md$shortname])))
      ## Be carefull about Inf and -Inf for prop = 0, 1
      fit_out <- fit_lm_interglht(data = logit_freq_out, md, formula = formula_lm, K = K)
      
    },
    lmer_arcsinesqrt_interglht = {
      
      ass_freq_out <- freq_out
      ass_freq_out[md$shortname] <- asin(sqrt((t(t(freq_out[md$shortname]) / colSums(freq_out[md$shortname])))))
      
      fit_out <- fit_lmer_interglht(data = ass_freq_out, md, formula = formula_lmer, K = K)
      
    },
    lm_arcsinesqrt_interglht = {
      
      ass_freq_out <- freq_out
      ass_freq_out[md$shortname] <- asin(sqrt((t(t(freq_out[md$shortname]) / colSums(freq_out[md$shortname])))))
      
      fit_out <- fit_lm_interglht(data = ass_freq_out, md, formula = formula_lm, K = K)
      
    },
    glmmadmb_fixed_beta_interglht = {
      fit_out <- fit_glm_interglht(data = freq_out, md, family = "beta", formula = formula_glm_beta, K = K)
    },
    glmmadmb_fixed_betabinomial_interglht = {
      fit_out <- fit_glm_interglht(data = freq_out, md, family = "betabinomial", formula = formula_glm_binomial, K = K)
    }
    
  )
  
  # ----------------------------------------
  # Extract p-values and coeffs
  # ----------------------------------------
  
  pvs <- data.frame(freq_out[, c("cluster", "label")], fit_out[["pvals"]])
  coeffs <- data.frame(freq_out[, c("cluster", "label")], fit_out[["coeffs"]])
  
  oo <- order(pvs[, pval_name1], decreasing = FALSE)
  pvs <- pvs[oo, , drop = FALSE]
  coeffs <- coeffs[oo, , drop = FALSE]
  
  ## save the results
  write.table(pvs, file=file.path(outdir, paste0(prefix, "frequencies_pvs_", k, suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  write.table(coeffs, file=file.path(outdir, paste0(prefix, "frequencies_coeffs_", k, suffix, ".xls")), row.names=FALSE, quote=FALSE, sep="\t")
  
  
  # ----------------------------------------
  # Plot a heatmap of significant cases - transform proportions with arcsin-sqrt so the dispersion is the same for low and high props.
  # ----------------------------------------

  ### add p-value info
  expr_all <- merge(pvs, expr_norm, by = c("cluster", "label"), all.x = TRUE, sort = FALSE)
  
  prefix2 <- paste0(k, "_")
  FDR_cutoff <- 0.05
  
  
  plot_heatmaps_for_sign_freqs(expr_all = expr_all, md = md, FDR_cutoff = FDR_cutoff, pval_name_list = pval_name_list, adjpval_name_list = adjpval_name_list, results2plot = results2plot, breaks = breaks, legend_breaks = legend_breaks, color_groups = color_groups, outdir = outdir, prefix = prefix, prefix2 = prefix2, suffix = suffix)
  
}






sessionInfo()













################################
### 04_frequencies done!
################################