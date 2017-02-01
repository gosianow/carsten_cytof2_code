#!/bin/bash

shopt -s expand_aliases
source ~/.bash_aliases

###############################################################################################################
## Define paths to software and reference files

RWD_MAIN=/home/Shared/data/cytof/carsten_cytof2
RCODE=/home/gosia/R/carsten_cytof2_code
METADATA=$RWD_MAIN/Espen_metadata
PANELS=$RWD_MAIN/Espen_panels

## Define which analysis to re-run
prepare_metadata=false
data_normalization=false
pcascores=false
select_observables=false
flowsom=false
flowsom_validation=false
heatmaps=false
runtsne=false
plottsne=true
plottsne_expr=false
cluster_merging=false

## global parameters
tsne_pmin="Inf" # In the CK analysis, I use 1500 per sample.

path_fun_prepare_metadata="0espen_prepare_metadata.R"
path_fun_formulas="00_formulas_espen.R"
outdir_fun_formulas="/"

###############################################################################################################
# Analysis of Espen data
# Use Analysis block 1
###############################################################################################################

data_dir="Espen"

file_panel="panel_espen_2016_12_22.xlsx"
file_metadata="metadata_espen_2016_12_22"

pca_score_cutoff=0 # We keep all the markers!
rand_seed_consensus=1234
nmetaclusts=20

prefix_data="espen_"
prefix_panel="espen_"
prefix_pca="pca0_"
prefix_clust="cl20_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --prepare_metadata ${prepare_metadata} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --flowsom_validation false --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --plottsne_expr ${plottsne_expr} --frequencies false --expression false --METADATA ${METADATA} --path_fun_prepare_metadata ${path_fun_prepare_metadata} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin} --path_fun_formulas ${path_fun_formulas} --outdir_fun_formulas ${outdir_fun_formulas}


for i in 6 7
do
  nmetaclusts=$i
  prefix_clust="cl${i}_"

  ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --prepare_metadata false --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --flowsom_validation false --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --plottsne_expr false --frequencies false --expression false --METADATA ${METADATA} --path_fun_prepare_metadata ${path_fun_prepare_metadata} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin} --path_fun_formulas ${path_fun_formulas} --outdir_fun_formulas ${outdir_fun_formulas}

done


















#
