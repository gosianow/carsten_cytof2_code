##############################################################################
## <<01_prepare_metadata.R>>

# BioC 3.3
# Created 29 Dec 2016


##############################################################################
Sys.time()
##############################################################################

# Load packages
library(gdata)
library(tools)

##############################################################################
# Test arguments
##############################################################################

rwd='/Users/gosia/Dropbox/UZH/carsten_cytof/CK_metadata'
file_metadata='metadata_23_01.xlsx'
path_fun_prepare_metadata='/Users/gosia/Dropbox/UZH/carsten_cytof_code/00_prepare_metadata.R'

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


# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

# read metadata
md <- read.xls(file_metadata, stringsAsFactors=FALSE)


source(path_fun_prepare_metadata)


md_new <- prepare_metadata(md)


saveRDS(md_new, file.path(paste0(basename(file_path_sans_ext(file_metadata)), ".rds")))

write.table(md_new, file.path(paste0(basename(file_path_sans_ext(file_metadata)), ".xls")), quote = FALSE, sep = "\t", row.names = FALSE)










sessionInfo()










################################
### Done!
################################