##############################################################################
## <<03_plottsne.R>>

# BioC 3.3
# Created 28 July 2016
# Updated 

##############################################################################
Sys.time()
##############################################################################

# Load packages
library(gdata)
library(Repitools)
library(gplots)
library(ggplot2)
library(plyr)
library(reshape2) # melt
library(RColorBrewer)
library(Rtsne)
library(coop) # cosine
library(limma) 


##############################################################################
# Test arguments
##############################################################################

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
tsnep_prefix='23_01_pca1_cl20_raw_'
tsnep_outdir='040_tsnemaps'
path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_01.rds'
path_rtsne_out='040_tsnemaps/23_01_pca1_raw_rtsne_out.rda'
path_rtsne_data='040_tsnemaps/23_01_pca1_raw_rtsne_data.xls'
path_clustering='030_heatmaps/23_01_pca1_cl20_clustering.xls'
path_clustering_labels='030_heatmaps/23_01_pca1_cl20_clustering_labels.xls'

### Optional arguments
tsne_cmin=1000 ### minimal size of a cluster to plot
tsne_distse=1
tsne_quantse=0.9

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

prefix <- tsnep_prefix
suffix <- ""
outdir <- tsnep_outdir

if(!file.exists(outdir)) 
  dir.create(outdir, recursive = TRUE)


if(!file.exists(path_rtsne_data))
  stop(paste0("File ", path_rtsne_data, " does not exists!!!"))


# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

md <- readRDS(path_metadata)


# ------------------------------------------------------------
# Load cluster data
# ------------------------------------------------------------

## clustering
clustering <- read.table(path_clustering, header = TRUE, sep = "\t", as.is = TRUE)

clust <- clustering[, "cluster"]

## clustering labels
labels <- read.table(path_clustering_labels, header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))
labels

# ------------------------------------------------------------
### Colors for tSNE maps
# ------------------------------------------------------------

# ------------------------------ 
# palette 1

# ggplot palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=60 , c=100)[1:n]
}

# ------------------------------ 
# color blind palette

colors_muted <- c("#DC050C", "#E8601C", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6", "#F1932D", "#F6C141", "#F7EE55", "#4EB265", "#90C987", "#CAEDAB")
color_ramp <- c(colors_muted, gg_color_hue(max(1, nlevels(labels$label) - length(colors_muted))))

colors_tsne <- color_ramp[1:nlevels(labels$label)]
names(colors_tsne) <- levels(labels$label)


# ------------------------------------------------------------
# Load tSNE data
# ------------------------------------------------------------

load(path_rtsne_out)

rtsne_data <- read.table(path_rtsne_data, header = TRUE, sep = "\t", as.is = TRUE)

## keep the tSNE cells that were used in clustering
cells2keep_rtsne <- rtsne_data$cell_index %in% clustering[, "cell_id"]


# ------------------------------------------------------------
# tSNE plots
# ------------------------------------------------------------

# get clustering for cells that were used in tSNE
names(clust) <- clustering[, "cell_id"]
clust_tsne <- clust[as.character(rtsne_data$cell_index[cells2keep_rtsne])]

ggdf <- data.frame(tSNE1 = rtsne_out$Y[cells2keep_rtsne, 1], tSNE2 = rtsne_out$Y[cells2keep_rtsne, 2], cluster = clust_tsne, sample = rtsne_data$sample_name[cells2keep_rtsne])

# add condition1 info
mm <- match(ggdf$sample, md$shortname)
ggdf$condition1 <- md$condition1[mm]
ggdf$condition2 <- md$condition2[mm]


# use cluster labels instead of numbers
mm <- match(ggdf$cluster, labels$cluster)
ggdf$cluster <- labels$label[mm]
ggdf$cluster <- factor(ggdf$cluster, levels = levels(labels$label))


# skipp the "drop" cluster

skipp_drop <- ggdf$cluster != "drop"

ggdf <- ggdf[skipp_drop, ]
ggdf$cluster <- factor(ggdf$cluster)

clust_tsne <- clust_tsne[skipp_drop]

rtsne_data <- rtsne_data[cells2keep_rtsne, ]
rtsne_data <- rtsne_data[skipp_drop, ]



# -----------------------------------
### Plot of tsne - all cells, all clusters

## one plot 
ggp <- ggplot(ggdf,  aes(x = tSNE1, y = tSNE2, color = cluster)) +
  geom_point(size = 1) +
  labs(x = "tSNE 1", y="tSNE 2")+ 
  theme_bw() +
  theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold"), legend.key = element_blank()) +
  scale_color_manual(values = colors_tsne[levels(ggdf$cluster)]) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

pdf(file.path(outdir, paste0(prefix, "tSNEone", suffix, ".pdf")), width = 9, height = 7)                 
print(ggp)
dev.off()


## facet per sample
ggp <- ggplot(ggdf,  aes(x = tSNE1, y = tSNE2, color = cluster)) +
  geom_point(size=1) +
  facet_wrap(~ sample) +
  labs(x = "tSNE 1", y="tSNE 2")+
  theme_bw() +
  theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold"), legend.key = element_blank()) +
  scale_color_manual(values = colors_tsne[levels(ggdf$cluster)]) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

pdf(file.path(outdir, paste0(prefix, "tSNEsample", suffix, ".pdf")), width = 10 + 2, height = 10)
print(ggp)
dev.off()


## facet per condition1
ggp <- ggplot(ggdf,  aes(x = tSNE1, y = tSNE2, color = cluster)) +
  geom_point(size=1) +
  facet_grid(. ~ condition1) +
  labs(x = "tSNE 1", y="tSNE 2")+
  theme_bw() +
  theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold"), legend.key = element_blank()) +
  scale_color_manual(values = colors_tsne[levels(ggdf$cluster)]) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

pdf(file.path(outdir, paste0(prefix, "tSNEcondition1", suffix, ".pdf")), width = 5 * nlevels(ggdf$condition1) + 2, height = 5)
print(ggp)
dev.off()


## facet per condition1 and condition2
ggp <- ggplot(ggdf,  aes(x = tSNE1, y = tSNE2, color = cluster)) +
  geom_point(size=1) +
  facet_grid(condition2 ~ condition1) +
  labs(x = "tSNE 1", y="tSNE 2")+
  theme_bw() +
  theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold"), legend.key = element_blank()) +
  scale_color_manual(values = colors_tsne[levels(ggdf$cluster)]) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

pdf(file.path(outdir, paste0(prefix, "tSNEcondition12", suffix, ".pdf")), width = 5 * nlevels(ggdf$condition1) + 2, height = 5 * nlevels(ggdf$condition2))
print(ggp)
dev.off()






### Save numbers of cells in each cluster

cell_count <- table(ggdf$cluster)

cell_count <- data.frame(cluster = names(cell_count), cell_count = as.numeric(cell_count))

write.table(cell_count, file = file.path(outdir, paste0(prefix, "tSNEone", suffix, ".xls")), quote = FALSE, sep = "\t", row.names = FALSE)


cell_count <- table(ggdf$cluster, ggdf$condition1)

cell_count <- data.frame(cluster = rownames(cell_count), as.data.frame.matrix(cell_count))

write.table(cell_count, file = file.path(outdir, paste0(prefix, "tSNEcondition1", suffix, ".xls")), quote = FALSE, sep = "\t", row.names = FALSE)






# -----------------------------------
### Plot of tsne - large clusters
if(any(grepl("tsne_cmin=", args))){
  
  tc <- table(ggdf$cluster)
  tc
  
  pdf(file.path(outdir, paste0(prefix, "cmin", suffix, ".pdf")), width = 10, height = 5)        
  barplot(tc)
  abline(h = tsne_cmin, col = "blue")
  dev.off()
  
  kk <- names(tc[tc > tsne_cmin])
  
  ggdf_sub <- ggdf[ggdf$cluster %in% kk,]
  ggdf_sub$cluster <- factor(ggdf_sub$cluster)
  
  
  ## one plot 
  ggp <- ggplot(ggdf_sub,  aes(x = tSNE1, y = tSNE2, color = cluster )) +
    geom_point(size=1) +
    labs(x = "tSNE 1", y="tSNE 2")+
    theme_bw() +
    theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold"), legend.key = element_blank()) +
    scale_color_manual(values = colors_tsne[levels(ggdf_sub$cluster)]) +
    guides(colour = guide_legend(override.aes = list(size = 5)))
  
  pdf(file.path(outdir, paste0(prefix, "tSNEone_subset", suffix, ".pdf")), width = 9, height = 7)
  print(ggp)
  dev.off()
  
  # facet per condition1
  ggp <- ggplot(ggdf_sub,  aes(x = tSNE1, y = tSNE2, color = cluster )) +
    geom_point(size=1) +
    facet_grid(. ~ condition1) +
    labs(x = "tSNE 1", y="tSNE 2")+
    theme_bw() +
    theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold"), legend.key = element_blank()) +
    scale_color_manual(values = colors_tsne[levels(ggdf_sub$cluster)]) +
    guides(colour = guide_legend(override.aes = list(size = 5)))
  
  pdf(file.path(outdir, paste0(prefix, "tSNEcondition1_subset", suffix, ".pdf")), width = 5 * nlevels(ggdf_sub$condition1) + 2, height = 5)
  print(ggp)
  dev.off()
  
  
  # facet per condition1 and condition2
  ggp <- ggplot(ggdf_sub,  aes(x = tSNE1, y = tSNE2, color = cluster )) +
    geom_point(size=1) +
    facet_grid(condition2 ~ condition1)+
    labs(x = "tSNE 1", y="tSNE 2")+
    theme_bw() +
    theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold"), legend.key = element_blank()) +
    scale_color_manual(values = colors_tsne[levels(ggdf_sub$cluster)]) +
    guides(colour = guide_legend(override.aes = list(size = 5)))
  
  pdf(file.path(outdir, paste0(prefix, "tSNEcondition12_subset", suffix, ".pdf")), width = 5 * nlevels(ggdf_sub$condition1) + 2, height = 5 * nlevels(ggdf_sub$condition2))
  print(ggp)
  dev.off()
  
}




# -----------------------------------
# Plot of tsne - cells that are close to the centers of clusters (no outliers), all clusters
# Approach 1 - removing all the cells that are further than 1 from its cluster center
if(any(grepl("tsne_distse=", args))){
  
  el <- rtsne_data[, -grep("cell_index|sample_name", colnames(rtsne_data)), drop = FALSE]
  
  a <- aggregate(el, by = list(clust_tsne), FUN=median)
  
  dists <- distse <- rep(NA, nrow(el))
  
  clust_lev <- a$Group.1
  
  for(i in 1:length(clust_lev)) {
    # i = 1
    
    data <- el[clust_tsne == clust_lev[i], , drop=FALSE]
    cent <- a[i,-1,drop=FALSE]
    
    ### cosine distance of cells from the center - per cluster
    # d <- 1 - cosine( t(rbind(cent,data)) )[1,-1]
    # dists[clust_tsne == clust_lev[i]] <- d
    
    ### Euclidean distance
    cent <- matrix(rep( as.numeric(cent), nrow(data)), byrow=TRUE, ncol=ncol(cent))
    # distse[clust_tsne == clust_lev[i]] <- sqrt(rowSums((cent-data)^2))
    distse[clust_tsne == clust_lev[i]] <- rowMeans(abs(cent-data))
  }
  
  
  ### Plot the distances per cluster
  ggdf$distse <- distse
  
  ggp <- ggplot(ggdf, aes(x = distse)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = tsne_distse, color = "blue") +
    facet_wrap(~ cluster, scales = "free_y")
  
  pdf(file.path(outdir, paste0(prefix, "distse", suffix, ".pdf")), width = 10, height = 10)        
  print(ggp)
  dev.off()
  
  
  ggdf_sub <- ggdf[distse < tsne_distse, ]
  ggdf_sub$cluster <- factor(ggdf_sub$cluster)
  
  
  ## one plot
  ggp <- ggplot(ggdf_sub,  aes(x = tSNE1, y = tSNE2, color = cluster)) +
    geom_point(size=1) + 
    labs(x = "tSNE 1", y="tSNE 2")+ 
    theme_bw() +
    theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold"), legend.key = element_blank())+
    scale_color_manual(values = colors_tsne[levels(ggdf_sub$cluster)]) +
    guides(colour = guide_legend(override.aes = list(size = 5)))
  
  pdf(file.path(outdir, paste0(prefix, "tSNEone_filtered", suffix, ".pdf")), width = 9, height = 7)                 
  print(ggp)
  dev.off()
  
  
  ## facet per condition1
  ggp <- ggplot(ggdf_sub,  aes(x = tSNE1, y = tSNE2, color = cluster)) +
    geom_point(size=1) +
    facet_grid(. ~ condition1) +
    labs(x = "tSNE 1", y="tSNE 2")+
    theme_bw() +
    theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold"), legend.key = element_blank())+
    scale_color_manual(values = colors_tsne[levels(ggdf_sub$cluster)]) +
    guides(colour = guide_legend(override.aes = list(size = 5)))
  
  pdf(file.path(outdir, paste0(prefix, "tSNEcondition1_filtered", suffix, ".pdf")), width = 5 * nlevels(ggdf_sub$condition1) + 2, height = 5)
  print(ggp)
  dev.off()
  
  
  ## facet per condition1 and condition2
  ggp <- ggplot(ggdf_sub,  aes(x = tSNE1, y = tSNE2, color = cluster)) +
    geom_point(size=1) +
    facet_grid(condition2 ~ condition1) +
    labs(x = "tSNE 1", y="tSNE 2")+
    theme_bw() +
    theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold"), legend.key = element_blank())+
    scale_color_manual(values = colors_tsne[levels(ggdf_sub$cluster)]) +
    guides(colour = guide_legend(override.aes = list(size = 5)))
  
  pdf(file.path(outdir, paste0(prefix, "tSNEcondition12_filtered", suffix, ".pdf")), width = 5 * nlevels(ggdf_sub$condition1) + 2, height = 5 * nlevels(ggdf_sub$condition2))
  print(ggp)
  dev.off()
  
}



# -----------------------------------
# Plot of tsne - cells that are close to the centers of clusters (no outliers), all clusters
# Approach 2 - removing 5% of most distant cells from the center in each cluster
if(any(grepl("tsne_quantse=", args))){
  
  el <- rtsne_data[, -grep("cell_index|sample_name", colnames(rtsne_data))]
  a <- aggregate(el, by = list(clust_tsne), FUN=median)
  
  clust_lev <- a$Group.1
  
  distse <- dist_keep <- rep(NA, nrow(el))
  
  for(i in 1:length(clust_lev)){
    # i = 1
    
    data <- el[clust_tsne == clust_lev[i], , drop=FALSE]
    cent <- a[i,-1,drop=FALSE]
    
    ### Euclidean distance
    cent <- matrix(rep( as.numeric(cent), nrow(data)), byrow=TRUE, ncol=ncol(cent))
    
    # d <- sqrt(rowSums((cent-data)^2))
    d <- rowMeans(abs(cent-data))
    distse[clust_tsne == clust_lev[i]] <- d
    
    dist_quant <- quantile(d, probs = tsne_quantse)
    
    dist_keep[clust_tsne == clust_lev[i]] <- d < dist_quant
    
  }
  
  ### Plot the distances per cluster
  ggdf$distse <- distse
  
  ggp <- ggplot(ggdf, aes(x = distse)) +
    geom_histogram(bins = 50) +
    facet_wrap(~ cluster, scales = "free_y")
  
  pdf(file.path(outdir, paste0(prefix, "quantse", suffix, ".pdf")), width = 10, height = 10)        
  print(ggp)
  dev.off()
  
  
  ggdf_sub <- ggdf[dist_keep, ]
  ggdf_sub$cluster <- factor(ggdf_sub$cluster)
  
  ## one plot
  ggp <- ggplot(ggdf_sub,  aes(x = tSNE1, y = tSNE2, color = cluster)) +
    geom_point(size=1) + 
    labs(x = "tSNE 1", y="tSNE 2")+ 
    theme_bw() +
    theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold"), legend.key = element_blank())+
    scale_color_manual(values = colors_tsne[levels(ggdf_sub$cluster)]) +
    guides(colour = guide_legend(override.aes = list(size = 5)))
  
  pdf(file.path(outdir, paste0(prefix, "tSNEone_filtered2", suffix, ".pdf")), width = 9, height = 7)                 
  print(ggp)
  dev.off()
  
  
  ## facet per condition1
  ggp <- ggplot(ggdf_sub,  aes(x = tSNE1, y = tSNE2, color = cluster)) +
    geom_point(size=1) +
    facet_grid(. ~ condition1) +
    labs(x = "tSNE 1", y="tSNE 2") +
    theme_bw() +
    theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold"), legend.key = element_blank())+
    scale_color_manual(values = colors_tsne[levels(ggdf_sub$cluster)]) +
    guides(colour = guide_legend(override.aes = list(size = 5)))
  
  pdf(file.path(outdir, paste0(prefix, "tSNEcondition1_filtered2", suffix, ".pdf")), width = 5 * nlevels(ggdf_sub$condition1) + 2, height = 5)
  print(ggp)
  dev.off()
  
  
  ## facet per condition1 and condition2
  ggp <- ggplot(ggdf_sub,  aes(x = tSNE1, y = tSNE2, color = cluster)) +
    geom_point(size=1) +
    facet_grid(condition2 ~ condition1) +
    labs(x = "tSNE 1", y="tSNE 2") +
    theme_bw() +
    theme(strip.text = element_text(size=15, face="bold"), axis.title  = element_text(size=15, face="bold"), legend.key = element_blank())+
    scale_color_manual(values = colors_tsne[levels(ggdf_sub$cluster)]) +
    guides(colour = guide_legend(override.aes = list(size = 5)))
  
  pdf(file.path(outdir, paste0(prefix, "tSNEcondition12_filtered2", suffix, ".pdf")), width = 5 * nlevels(ggdf_sub$condition1) + 2, height = 5 * nlevels(ggdf_sub$condition2))
  print(ggp)
  dev.off()
  
}







sessionInfo()















################################
### 03_plottsne done!
################################