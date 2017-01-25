#!/bin/bash

shopt -s expand_aliases
source ~/.bash_aliases

###############################################################################################################
## Define paths to software and reference files

RWD_MAIN=/home/Shared/data/cytof/carsten_cytof
RCODE=/home/gosia/R/carsten_cytof_code
METADATA=$RWD_MAIN/MyeEUNITERfinal_metadata
PANELS=$RWD_MAIN/MyeEUNITERfinal_panels

## Define which analysis to re-run
prepare_metadata=false
data_normalization=false
pcascores=false
select_observables=false
flowsom=false
flowsom_validation=false
heatmaps=false
runtsne=false
plottsne=false
plottsne_expr=false
frequencies=false
expression=false
cluster_merging=false
cluster_extracting=false

## global parameters
tsne_pmin=2000 # In the CK analysis, I use 1500 per sample.

path_fun_prepare_metadata="0myef_prepare_metadata.R"
path_fun_formulas="00_formulas_1dataset_3responses.R"
outdir_fun_formulas="3responses"

###############################################################################################################
# Analysis of MyeEUNITERfinal data
# Use Analysis block 1
###############################################################################################################

data_dir="MyeEUNITERfinal"

file_panel="panel_MyeEUNITERfinal.xlsx"
file_metadata="metadata_MyeEUNITERfinal"

pca_score_cutoff=0 # We keep all the markers!
rand_seed_consensus=1234
nmetaclusts=20

prefix_data="myef_"
prefix_panel="myef_"
prefix_pca="pca0_"
prefix_clust="cl20_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --prepare_metadata ${prepare_metadata} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --flowsom_validation false --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --plottsne_expr ${plottsne_expr} --frequencies false --expression false --METADATA ${METADATA} --path_fun_prepare_metadata ${path_fun_prepare_metadata} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin} --path_fun_formulas ${path_fun_formulas} --outdir_fun_formulas ${outdir_fun_formulas}


for i in 11
do
  nmetaclusts=$i
  prefix_clust="cl${i}_"

  ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --prepare_metadata false --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --flowsom_validation false --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --plottsne_expr false --frequencies false --expression false --METADATA ${METADATA} --path_fun_prepare_metadata ${path_fun_prepare_metadata} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin} --path_fun_formulas ${path_fun_formulas} --outdir_fun_formulas ${outdir_fun_formulas}

done


# --------------------------------------------------
# Analysis of MyeEUNITERfinal cluster_merging
# Use Analysis block 2
# --------------------------------------------------

prefix_clust="cl20_"

file_merging="010_helpfiles/${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_merging.xlsx"
prefix_merging="merging_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging} --path_fun_formulas ${path_fun_formulas} --outdir_fun_formulas ${outdir_fun_formulas}


# --------------------------------------------------
# MyeEUNITERfinal - neutrophils and monocytes cluster extracting from merging
# Use Analysis block 3
# --------------------------------------------------

prefix_merging="merging_"

extract_cluster="c('neutrophils','monocytes')"
extract_dir="MyeEUNITERfinal_neutrophils_merging"


./Analysis_block_3_cluster_extracting.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_extracting ${cluster_extracting} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging} --extract_cluster ${extract_cluster} --extract_dir ${extract_dir}


# --------------------------------------------------
# Analysis of MyeEUNITERfinal_neutrophils_merging
# Use Analysis block 1
# --------------------------------------------------

file_metadata="metadata_MyeEUNITERfinal"

rand_seed_consensus=1234
nmetaclusts=20

prefix_pca="pca0_"
prefix_clust="cl20_"

data_dir='MyeEUNITERfinal_neutrophils_merging'
prefix_data='myefNEUTROP_'
file_panel='panel_MyeEUNITERfinal.xlsx'
prefix_panel='myef_'
pca_score_cutoff=0


./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --prepare_metadata false --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --flowsom_validation false --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --plottsne_expr false --frequencies ${frequencies} --expression false --METADATA ${METADATA} --path_fun_prepare_metadata ${path_fun_prepare_metadata} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin} --path_fun_formulas ${path_fun_formulas} --outdir_fun_formulas ${outdir_fun_formulas}

















#
