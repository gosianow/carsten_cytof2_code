##############################################################################
## <<02_cluster_extracting.R>>

# BioC 3.3
# Created 28 July 2016
# Updated 25 Aug 2016

##############################################################################
Sys.time()
##############################################################################

library(flowCore)
library(gdata)


##############################################################################
# Test arguments
##############################################################################

# rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01'
# extract_outdir='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_2016-06-23_01_CD4/010_cleanfcs'
# path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata/metadata_23_01.xlsx'
# path_clustering='030_heatmaps/23_01_pca1_mergingNEW2_clustering.xls'
# path_clustering_labels='030_heatmaps/23_01_pca1_mergingNEW2_clustering_labels.xls'
# extract_cluster='CD4'

rwd='/home/Shared/data/cytof/carsten_cytof/MyeEUNITERfinal'
extract_outdir='/home/Shared/data/cytof/carsten_cytof/MyeEUNITERfinal_neutrophils_merging/010_cleanfcs'
path_metadata='/home/Shared/data/cytof/carsten_cytof/MyeEUNITERfinal_metadata/metadata_MyeEUNITERfinal.xlsx'
path_clustering='030_heatmaps/myef_myef_pca0_merging_clustering.xls'
path_clustering_labels='030_heatmaps/myef_myef_pca0_merging_clustering_labels.xls'
extract_cluster=c('neutrophils','monocytes')


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

outdir <- extract_outdir

if(!file.exists(outdir)) 
  dir.create(outdir, recursive = TRUE)

fcsDir <- "010_cleanfcs"

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

# read metadata
md <- read.xls(path_metadata, stringsAsFactors=FALSE)


# define FCS file names
f <- file.path(fcsDir, md$filename)
names(f) <- md$shortname


# read raw FCS files in
fcs <- lapply(f, read.FCS)



# ------------------------------------------------------------
# Load clustering data
# ------------------------------------------------------------

clustering <- read.table(path_clustering, header = TRUE, sep = "\t", as.is = TRUE)
clust <- clustering[, "cluster"]

labels <- read.table(path_clustering_labels, header = TRUE, sep = "\t", as.is = TRUE)
labels <- labels[order(labels$cluster, decreasing = FALSE), ]
labels$label <- factor(labels$label, levels = unique(labels$label))


# ------------------------------------------------------------
# Save fcs files for a cluster to export
# ------------------------------------------------------------

## use labels as cluster names
mlab <- match(clust, labels$cluster)
clust <- labels$label[mlab]

clustList <- split(clust, factor(clustering[, "sample_id"], levels = names(fcs)))


if(!all(names(fcs) == names(clustList)))
  stop("Wrong order of samples!")


lapply(1:length(fcs), function(i){
  # i = 2
  
  fn <- file.path(outdir, basename(f[i]))
  
  cells2keep <- clustList[[i]] %in% extract_cluster
  
  out <- fcs[[i]][cells2keep, ]
  
  write.FCS(out, fn)
  
  return(table(cells2keep))
  
})






sessionInfo()












################################
### 02_cluster_extracting done!
################################
