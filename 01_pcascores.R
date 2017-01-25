##############################################################################
## <<01_pcascores.R>>

# BioC 3.3
# Created 27 July 2016
# Updated 31 Aug 2016

##############################################################################
Sys.time()
##############################################################################

# Load packages
library(gdata)
library(gplots)
library(ggplot2)
library(plyr)
library(reshape2)
library(coop)
library(RColorBrewer)

##############################################################################
# Test arguments
##############################################################################

rwd='/home/Shared/data/cytof/carsten_cytof/MyeEUNITERfinal'
pcas_prefix='myef_myef_'
pcas_outdir='020_pcascores'
path_data='010_data/myef_myef_expr_raw.rds'
path_metadata='/home/Shared/data/cytof/carsten_cytof/MyeEUNITERfinal_metadata/metadata_MyeEUNITERfinal.rds'
path_panel='/home/Shared/data/cytof/carsten_cytof/MyeEUNITERfinal_panels/panel_MyeEUNITERfinal.xlsx'


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

prefix <- pcas_prefix

outdir <- pcas_outdir

if( !file.exists(outdir) ) 
  dir.create(outdir, recursive = TRUE)

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

samp <- expr[, "sample_id"]
fcs_colnames <- colnames(expr)[!grepl("cell_id|sample_id", colnames(expr))]
e <- expr[, fcs_colnames]


# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

md <- readRDS(path_metadata)

color_samples <- md$color
names(color_samples) <- md$shortname

# ------------------------------------------------------------
# Load panel
# ------------------------------------------------------------

# read panel, pick which columns to use
panel <- read.xls(path_panel, stringsAsFactors=FALSE)

if(!all(c("fcs_colname", "Isotope", "Antigen", "transform") %in% colnames(panel)))
  stop("Wrong columns in panel!!!")


# get the isotope and antigen for fcs markers
m <- match(fcs_colnames, panel$fcs_colname)

fcs_panel <- data.frame(fcs_colname = fcs_colnames, Isotope = panel$Isotope[m], Antigen = panel$Antigen[m], stringsAsFactors = FALSE)


# --------------------------------------------------------------------------
# run Levine et al. 2015 marker scoring,
# plot marginal distributions
# --------------------------------------------------------------------------

min_samp <- nrow(md)


doPRINCOMP <- function(z, ncomp=3) {
  # z = es[[1]]; ncomp=3
  
  if(nrow(z) < min_samp)
    return(rep(NA, ncol(z)))
  
  ## Score by Levine
  pr <- prcomp(z, center = TRUE, scale. = FALSE)  # default is to center but not scale
  
  score <- rowSums (outer( rep(1, ncol(z)), pr$sdev[1:ncomp]^2 ) * abs(pr$rotation[,1:ncomp]) )
  
  return(score)
  
}


# split expression per sample
es <- split(e, samp)

prs <- sapply(es, doPRINCOMP)
rmprs <- rowMeans(prs, na.rm = TRUE)

prs <- data.frame(mass = fcs_panel$fcs_colname, marker = fcs_panel$Antigen, avg_score = rmprs, round(prs, 4))

### Save ordered PCA scores
o <- order(rmprs, decreasing=TRUE)
prs <- prs[o,]

write.table(prs, file = file.path(outdir, paste0(prefix, "princompscore_by_sample.xls")), sep="\t", row.names=FALSE, quote=FALSE)


## plot PCA scores

ggp1 <- ggplot(prs, aes(x = marker, y = avg_score)) +
  geom_bar(stat="identity") +
  theme_bw() +
  ylab("Average PCA score") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=12))


## plot PCA scores for ordered markers

prs$marker <- factor(prs$marker, levels = prs$marker)

ggp2 <- ggplot(prs, aes(x = marker, y = avg_score)) +
  geom_bar(stat="identity") +
  theme_bw() +
  ylab("Average PCA score") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=12))


pdf(file.path(outdir, paste0(prefix, "princompscore_average.pdf")), width = 10, height = 7)
print(ggp1)
print(ggp2)
dev.off()


### Plot chanel distributions


df <- data.frame(samp = samp, e, check.names = FALSE)
dfm <- melt(df, id.var = "samp")

mm <- match(dfm$variable, prs$mass)

dfm$marker <- prs$marker[mm]

dfm$variable_marker <- interaction(dfm$variable, dfm$marker, sep = "/")

dfm$variable_marker <- factor(dfm$variable_marker, levels = paste0(prs$mass, "/", prs$marker))

variable_markers <- levels(dfm$variable_marker)



pdf(file.path(outdir, paste0(prefix, "channel_distributions.pdf")))

for(i in 1:length(variable_markers)){
  
  ggdf <- dfm[dfm$variable_marker == variable_markers[i], , drop = FALSE]
  
  ggp <- ggplot(ggdf, aes(x = value, color = samp)) +
    geom_density(adjust = 1) +
    ggtitle(variable_markers[i]) +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_color_manual(values = color_samples)
  
  print(ggp)
  
}

dev.off()








sessionInfo()













################################
### 01_pcascores done!
################################