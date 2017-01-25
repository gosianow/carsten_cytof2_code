##############################################################################

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


rwd='/home/Shared/data/cytof/carsten_cytof/MyeEUNITERfinal'
data_prefix='myef_'
data_outdir='/home/Shared/data/cytof/carsten_cytof/MyeEUNITERfinal_panels'
path_metadata='/home/Shared/data/cytof/carsten_cytof/MyeEUNITERfinal_metadata/metadata_MyeEUNITERfinal.xlsx'
path_panel='/home/Shared/data/cytof/carsten_cytof/MyeEUNITERfinal_panels/CK_panel_2016-12-22.xlsx'


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