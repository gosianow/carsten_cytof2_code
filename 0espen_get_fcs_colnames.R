##############################################################################
## <<01_get_fcs_colnames_espen.R>>

# BioC 3.3
# Created 29 Dec 2016


##############################################################################
Sys.time()
##############################################################################

# Load packages
library(flowCore)
library(gdata)

##############################################################################
# Test arguments
##############################################################################


rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/Espen'
data_prefix='espen_'
data_outdir='/Users/gosia/Dropbox/UZH/carsten_cytof/Espen_panels'
path_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof/Espen_metadata/espen_metadata_2016_12_22.xlsx'
path_panel='/Users/gosia/Dropbox/UZH/carsten_cytof/Espen_panels/espen_panel_2016_12_22.xlsx'


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

prefix <- data_prefix

# ------------------------------------------------------------
# define directories
# ------------------------------------------------------------

fcsDir <- "010_cleanfcs"
outdir <- data_outdir

if(!file.exists(outdir)) 
  dir.create(outdir, recursive = TRUE)


# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

# read metadata
md <- read.xls(path_metadata, stringsAsFactors=FALSE)


# ------------------------------------------------------------
# Load fcs files
# ------------------------------------------------------------

# define FCS file names
f <- file.path(fcsDir, md$filename)
names(f) <- md$shortname

# read raw FCS files in
fcs <- lapply(f, read.FCS)

fcs_colnames <- colnames(fcs[[1]])
fcs_colnames


# ------------------------------------------------------------
# Load panel
# ------------------------------------------------------------

# read panel, pick which columns to use
panel <- read.xls(path_panel, stringsAsFactors=FALSE)

if(!all(c("Isotope", "Antigen", "transform") %in% colnames(panel)))
  stop("Wrong columns in panel!!!")


### Get the fcs_colname


fcs_isotope <- as.numeric(gsub("[^0-9]", "", fcs_colnames))

m <- match(panel$Isotope, fcs_isotope)

panel$fcs_colname <- fcs_colnames[m]


write.table(panel, file.path(outdir, paste0(prefix, "panel.xls")), quote = FALSE, sep = "\t", row.names = FALSE)







sessionInfo()













################################
### done!
################################