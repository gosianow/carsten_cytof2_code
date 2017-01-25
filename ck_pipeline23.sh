#!/bin/bash

shopt -s expand_aliases
source ~/.bash_aliases

###############################################################################################################
## Define paths to software and reference files

RWD_MAIN=/home/Shared/data/cytof/carsten_cytof
RCODE=/home/gosia/R/carsten_cytof_code
METADATA=$RWD_MAIN/CK_metadata
PANELS=$RWD_MAIN/CK_panels


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
fcs_saving=false
cytokines_bimatrix=false
cytokines_bimatrix_main=false
pd1_bimatrix=false
pd1_bimatrix_main=false
cytokines_expression=false
pd1_expression=false
cd69_bimatrix=false
cd69_bimatrix_main=false

## global parameters
tsne_pmin=1500

path_fun_prepare_metadata="00_prepare_metadata.R"
path_fun_formulas="00_formulas_1dataset_3responses.R"
outdir_fun_formulas="3responses"

###############################################################################################################
# Analysis of CK_2016-06-23_01 data
# Use Analysis block 1
###############################################################################################################

DATA=23
PANEL=1
data_dir="CK_2016-06-23_01"

file_panel="panel1.xlsx"
file_metadata="metadata_23_01"

pca_score_cutoff=3
rand_seed_consensus=123
nmetaclusts=20

prefix_data="23_"
prefix_panel="01_"
prefix_pca="pca1_"
prefix_clust="cl20_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --prepare_metadata ${prepare_metadata} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --flowsom_validation ${flowsom_validation} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --plottsne_expr ${plottsne_expr} --frequencies ${frequencies} --expression false --METADATA ${METADATA} --path_fun_prepare_metadata ${path_fun_prepare_metadata} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin} --path_fun_formulas ${path_fun_formulas} --outdir_fun_formulas ${outdir_fun_formulas}


for i in 7 8
do
  nmetaclusts=$i
  prefix_clust="cl${i}_"

  ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --prepare_metadata false --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --flowsom_validation false --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --plottsne_expr false --frequencies ${frequencies} --expression false --METADATA ${METADATA} --path_fun_prepare_metadata ${path_fun_prepare_metadata} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin} --path_fun_formulas ${path_fun_formulas} --outdir_fun_formulas ${outdir_fun_formulas}

done


# --------------------------------------------------
# Analysis of CK_2016-06-23_01 cluster_merging
# Use Analysis block 2
# --------------------------------------------------

prefix_clust="cl20_"

file_merging="010_helpfiles/${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_mergingNEW2.xlsx"
prefix_merging="mergingNEW2_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging} --path_fun_formulas ${path_fun_formulas} --outdir_fun_formulas ${outdir_fun_formulas}


file_merging="010_helpfiles/${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_merging6.xlsx"
prefix_merging="merging6_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging} --path_fun_formulas ${path_fun_formulas} --outdir_fun_formulas ${outdir_fun_formulas}


# --------------------------------------------------
# CK_2016-06-23_01 - CD4 and CD8 cluster extracting from mergingNEW2
# Use Analysis block 3
# --------------------------------------------------

prefix_merging="mergingNEW2_"

extract_cluster=("'CD4'" "'CD8'")
extract_dir=('CK_2016-06-23_01_CD4_mergingNEW2' 'CK_2016-06-23_01_CD8_mergingNEW2')

for i in 0 1
do
  ./Analysis_block_3_cluster_extracting.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_extracting ${cluster_extracting} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging} --extract_cluster ${extract_cluster[$i]} --extract_dir ${extract_dir[$i]}
done


# --------------------------------------------------
# Analysis of CK_2016-06-23_01_CD4_mergingNEW2 using panel1CD4.xlsx
# and CK_2016-06-23_01_CD8_mergingNEW2 using panel1CD8.xlsx
# Use Analysis block 1
# --------------------------------------------------

DATA=23
PANEL=1

file_metadata="metadata_23_01"

rand_seed_consensus=1234
nmetaclusts=20

prefix_pca="pca1_"
prefix_clust="cl20_"

data_dir=('CK_2016-06-23_01_CD4_mergingNEW2' 'CK_2016-06-23_01_CD8_mergingNEW2')
prefix_data=('23CD4_' '23CD8_')
file_panel=('panel1CD4.xlsx' 'panel1CD8.xlsx')
prefix_panel=('01CD4_' '01CD8_')
pca_score_cutoff=(2 2)

for i in 0 1
do
  ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --prepare_metadata false --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --flowsom_validation ${flowsom_validation} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --plottsne_expr ${plottsne_expr} --frequencies ${frequencies} --expression false --METADATA ${METADATA} --path_fun_prepare_metadata ${path_fun_prepare_metadata} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff[$i]} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin} --path_fun_formulas ${path_fun_formulas} --outdir_fun_formulas ${outdir_fun_formulas}

done


rand_seed_consensus=1234
nmetaclusts=(5 8)
prefix_clust=("cl5_" "cl8_")

for i in 0 1
do
  ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --prepare_metadata false --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --flowsom_validation false --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --plottsne_expr false --frequencies ${frequencies} --expression false --METADATA ${METADATA} --path_fun_prepare_metadata ${path_fun_prepare_metadata} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust[$i]} --pca_score_cutoff ${pca_score_cutoff[$i]} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts[$i]} --tsne_pmin ${tsne_pmin} --path_fun_formulas ${path_fun_formulas} --outdir_fun_formulas ${outdir_fun_formulas}

done



# --------------------------------------------------
# Analysis of CK_2016-06-23_01_CD4_mergingNEW2 using panel1CD4.xlsx for cluster_merging
# and CK_2016-06-23_01_CD8_mergingNEW2 using panel1CD8.xlsx for cluster_merging
# Use Analysis block 2
# --------------------------------------------------


prefix_clust=('cl5_' 'cl8_')

file_merging=("010_helpfiles/${prefix_data[0]}${prefix_panel[0]}${prefix_pca}${prefix_clust[0]}cluster_merging5.xlsx" "010_helpfiles/${prefix_data[1]}${prefix_panel[1]}${prefix_pca}${prefix_clust[1]}cluster_merging5.xlsx")
prefix_merging=('merging5_' 'merging5_')

for i in  0 1
do
  ./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust[$i]} --prefix_merging ${prefix_merging[$i]} --file_merging ${file_merging[$i]} --path_fun_formulas ${path_fun_formulas} --outdir_fun_formulas ${outdir_fun_formulas}
done


###############################################################################################################
# Analysis of CK_2016-06-23_02 data
# Use Analysis block 1
###############################################################################################################

DATA=23
PANEL=2
data_dir="CK_2016-06-23_02"

file_panel="panel2.xlsx"
file_metadata="metadata_23_02"

pca_score_cutoff=1
rand_seed_consensus=1234
nmetaclusts=20

prefix_data="23_"
prefix_panel="02_"
prefix_pca="pca1_"
prefix_clust="cl20_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --prepare_metadata ${prepare_metadata} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --flowsom_validation ${flowsom_validation} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --plottsne_expr ${plottsne_expr} --frequencies ${frequencies} --expression false --METADATA ${METADATA} --path_fun_prepare_metadata ${path_fun_prepare_metadata} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin} --path_fun_formulas ${path_fun_formulas} --outdir_fun_formulas ${outdir_fun_formulas}

rand_seed_consensus=1234

for i in 5 8
do

  nmetaclusts=$i
  prefix_clust="cl${i}_"

  ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --prepare_metadata false --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --flowsom_validation false --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --plottsne_expr false --frequencies ${frequencies} --expression false --METADATA ${METADATA} --path_fun_prepare_metadata ${path_fun_prepare_metadata} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin} --path_fun_formulas ${path_fun_formulas} --outdir_fun_formulas ${outdir_fun_formulas}

done


# --------------------------------------------------
# Analysis of CK_2016-06-23_02 merging2
# Use Analysis block 2
# --------------------------------------------------

prefix_clust="cl20_"

file_merging="010_helpfiles/${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_merging2.xlsx"
prefix_merging="merging2_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging} --path_fun_formulas ${path_fun_formulas} --outdir_fun_formulas ${outdir_fun_formulas}



# --------------------------------------------------
# CK_2016-06-23_02 - CD4 and CD8 cluster extracting from merging2
# Use Analysis block 3
# --------------------------------------------------

prefix_merging="merging2_"

extract_cluster=("'CD4'" "'CD8'")
extract_dir=('CK_2016-06-23_02_CD4_merging2' 'CK_2016-06-23_02_CD8_merging2')

for i in 0 1
do
  ./Analysis_block_3_cluster_extracting.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_extracting ${cluster_extracting} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging} --extract_cluster ${extract_cluster[$i]} --extract_dir ${extract_dir[$i]}
done

# --------------------------------------------------
# Analysis of CK_2016-06-23_02_CD4_merging2 using panel2CD4.xlsx
# and CK_2016-06-23_02_CD8_merging2 using panel2CD8.xlsx
# Use Analysis block 1
# --------------------------------------------------

DATA=23
PANEL=2
file_metadata="metadata_23_02"

rand_seed_consensus=1234
nmetaclusts=20

prefix_pca="pca1_"
prefix_clust="cl20_"

data_dir=('CK_2016-06-23_02_CD4_merging2' 'CK_2016-06-23_02_CD8_merging2')
prefix_data=('23CD4_' '23CD8_')
file_panel=('panel2CD4.xlsx' 'panel2CD8.xlsx')
prefix_panel=('02CD4_' '02CD8_')
pca_score_cutoff=(0.96 1)

for i in 0 1
do
  ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --prepare_metadata false --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --flowsom_validation ${flowsom_validation} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --plottsne_expr ${plottsne_expr} --frequencies ${frequencies} --expression false --METADATA ${METADATA} --path_fun_prepare_metadata ${path_fun_prepare_metadata} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff[$i]} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin} --path_fun_formulas ${path_fun_formulas} --outdir_fun_formulas ${outdir_fun_formulas}
done


rand_seed_consensus=1234
nmetaclusts=(6 7)
prefix_clust=("cl6_" "cl7_")

for i in 0 1
do
  ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --prepare_metadata false --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --flowsom_validation false --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --plottsne_expr false --frequencies ${frequencies} --expression false --METADATA ${METADATA} --path_fun_prepare_metadata ${path_fun_prepare_metadata} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust[$i]} --pca_score_cutoff ${pca_score_cutoff[$i]} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts[$i]} --tsne_pmin ${tsne_pmin} --path_fun_formulas ${path_fun_formulas} --outdir_fun_formulas ${outdir_fun_formulas}
done

# --------------------------------------------------
# Analysis of CK_2016-06-23_02_CD4_merging2 using panel2CD4.xlsx for cluster_merging
# and CK_2016-06-23_02_CD8_merging2 using panel2CD8.xlsx for cluster_merging
# Use Analysis block 2
# --------------------------------------------------


prefix_clust=('cl6_' 'cl7_')

file_merging=("010_helpfiles/${prefix_data[0]}${prefix_panel[0]}${prefix_pca}${prefix_clust[0]}cluster_merging2.xlsx" "010_helpfiles/${prefix_data[1]}${prefix_panel[1]}${prefix_pca}${prefix_clust[1]}cluster_merging2.xlsx")
prefix_merging=('merging2_' 'merging2_')

for i in 0 1
do
  ./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust[$i]} --prefix_merging ${prefix_merging[$i]} --file_merging ${file_merging[$i]} --path_fun_formulas ${path_fun_formulas} --outdir_fun_formulas ${outdir_fun_formulas}
done


# --------------------------------------------------
# FCS (transformed with arcsineh and arcsineh01) saving for CK_2016-06-29_02_CD4 and CK_2016-06-29_02_CD8
# Then one can define the thresholds for positive cytokines
# --------------------------------------------------


### FCS saving
for i in 0 1
do
  ./Analysis_block_3_fcs_saving.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --fcs_saving ${fcs_saving} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]}
done



# --------------------------------------------------
# Analysis of cytokines based on bimatrix
# in CK_2016-06-23_02_CD4_merging2 and CK_2016-06-23_02_CD8_merging2
# Use Analysis block 4
# --------------------------------------------------
# Clustering of bimatrix is based on SOM only when som_dim^2 = nmetaclusts, otherwise cluster consesnsus is applied additionally

# -----------------------------
### Tmem cluster
# -----------------------------

DATA=23
PANEL=2
file_metadata="metadata_23_02"

prefix_pca="pca1_"
prefix_clust="cl20_"

data_dir=('CK_2016-06-23_02_CD4_merging2' 'CK_2016-06-23_02_CD8_merging2')
prefix_data=('23CD4_' '23CD8_')
file_panel=('panel2CD4.xlsx' 'panel2CD8.xlsx')
prefix_panel=('02CD4_' '02CD8_')

prefix_merging=('merging2_' 'merging2_')
clsubset=("c('CM','EM','TE','TM')" "c('CM','EM','TE','TM')")
prefix_clsubset=('Tmem_' 'Tmem_')

outdir="060_cytokines_bimatrix"


### for CD8, new cytokine cutoffs are defined - they include all the cytokines that were also used for CD4
file_cytokines_cutoffs=('panel2CD4_23_cytokines_CM.xlsx' 'panel2CD8_23_cytokines_CM_new.xlsx')
prefix_cytokines_cutoffs=('cytCM_raw2_' 'cytCMnew_raw2_')

som_dim=(5 5)
nmetaclusts=(25 25)

for i in 0 1
do
  ./Analysis_block_4_cytokines_bimatrix_0_bimatrix.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cytokines_bimatrix ${cytokines_bimatrix} --PANELS ${PANELS} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging[$i]} --prefix_clsubset ${prefix_clsubset[$i]} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs[$i]} --file_cytokines_cutoffs ${file_cytokines_cutoffs[$i]} --clsubset ${clsubset[$i]} --outdir ${outdir}

  ./Analysis_block_6_cytokines_bimatrix_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cytokines_bimatrix_main ${cytokines_bimatrix_main} --PANELS ${PANELS} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging[$i]} --prefix_clsubset ${prefix_clsubset[$i]} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs[$i]}  --som_dim ${som_dim[$i]} --nmetaclusts ${nmetaclusts[$i]} --outdir ${outdir}
done


# --------------------------------------------------
# Analysis of cytokines based on expression
# in CK_2016-06-23_02_CD4_merging2 and CK_2016-06-23_02_CD8_merging2
# Use Analysis block 4
# --------------------------------------------------

DATA=23
PANEL=2
file_metadata="metadata_23_02"

prefix_pca="pca1_"

data_dir=('CK_2016-06-23_02_CD4_merging2' 'CK_2016-06-23_02_CD8_merging2')
prefix_data=('23CD4_' '23CD8_')
prefix_panel=('02CD4_' '02CD8_')

prefix_merging=('merging2_' 'merging2_') # name of merging from which the Tmem clusters are extracted
clsubset=("c('CM','EM','TE','TM')" "c('CM','EM','TE','TM')")
prefix_clsubset=('Tmem_' 'Tmem_')
file_cytokines_cutoffs=('panel2CD4_23_cytokines_CM.xlsx' 'panel2CD8_23_cytokines_CM.xlsx')
prefix_cytokines_cutoffs=('cytCM_raw2_' 'cytCM_raw2_')


for i in 0 1
do
  ./Analysis_block_4_cytokines_expression.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cytokines ${cytokines_expression} --PANELS ${PANELS} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging[$i]} --prefix_clsubset ${prefix_clsubset[$i]} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs[$i]} --file_cytokines_cutoffs ${file_cytokines_cutoffs[$i]} --clsubset ${clsubset[$i]}
done



# ----------------------------------------------------------------------------------------------------
# Analysis of PD1+ and PD1- cells based on bimatrix
# in CK_2016-06-23_02_CD4_merging2 and CK_2016-06-23_02_CD8_merging2
# Use Analysis block 5
# ----------------------------------------------------------------------------------------------------

### Analysis of cytokines for PD1+

DATA=23
PANEL=2
file_metadata="metadata_23_02"

prefix_pca="pca1_"
prefix_clust="cl20_"

data_dir=('CK_2016-06-23_02_CD4_merging2' 'CK_2016-06-23_02_CD8_merging2')
prefix_data=('23CD4_' '23CD8_')
file_panel=('panel2CD4.xlsx' 'panel2CD8.xlsx')
prefix_panel=('02CD4_' '02CD8_')

prefix_merging=('merging2_' 'merging2_')
clsubset=("c('CM','EM','TE','TM')" "c('CM','EM','TE','TM')")
prefix_clsubset=('Tmem_' 'Tmem_')
file_cytokines_cutoffs=('panel2CD4_23_cytokines_CM.xlsx' 'panel2CD8_23_cytokines_CM_new.xlsx')
marker="PD-1"
prefix_cytokines_cutoffs=('cytCM_raw2_pd1_' 'cytCMnew_raw2_pd1_')

outdir="070_pd1_bimatrix"

som_dim=(5 5)
nmetaclusts=(25 25)

for i in 0 1
do
  ./Analysis_block_5_pd1_bimatrix_0_bimatrix.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --pd1_bimatrix ${pd1_bimatrix} --PANELS ${PANELS} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging[$i]} --prefix_clsubset ${prefix_clsubset[$i]} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs[$i]} --file_cytokines_cutoffs ${file_cytokines_cutoffs[$i]} --clsubset ${clsubset[$i]} --marker ${marker} --outdir ${outdir}

  ./Analysis_block_6_cytokines_bimatrix_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cytokines_bimatrix_main ${pd1_bimatrix_main} --PANELS ${PANELS} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging[$i]} --prefix_clsubset ${prefix_clsubset[$i]} --prefix_cytokines_cutoffs "${prefix_cytokines_cutoffs[$i]}positive_"  --som_dim ${som_dim[$i]} --nmetaclusts ${nmetaclusts[$i]} --outdir ${outdir}

done


# --------------------------------------------------
# Analysis of PD1+ and PD1- cells based on expression
# in CK_2016-06-23_02_CD4_merging2 and CK_2016-06-23_02_CD8_merging2
# Use Analysis block 5
# --------------------------------------------------

### Analysis of cytokines for PD1+

DATA=23
PANEL=2
file_metadata="metadata_23_02"

prefix_pca="pca1_"

data_dir=('CK_2016-06-23_02_CD4_merging2' 'CK_2016-06-23_02_CD8_merging2')
prefix_data=('23CD4_' '23CD8_')
prefix_panel=('02CD4_' '02CD8_')

prefix_merging=('merging2_' 'merging2_') # name of merging from which the Tmem clusters are extracted
clsubset=("c('CM','EM','TE','TM')" "c('CM','EM','TE','TM')")
prefix_clsubset=('Tmem_' 'Tmem_')
file_cytokines_cutoffs=('panel2CD4_23_cytokines_CM.xlsx' 'panel2CD8_23_cytokines_CM_new.xlsx')
prefix_cytokines_cutoffs=('cytCM_raw2_pd1_' 'cytCMnew_raw2_pd1_')


for i in 0 1
do
  ./Analysis_block_5_pd1_expression.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --pd1 ${pd1_expression} --PANELS ${PANELS} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging[$i]} --prefix_clsubset ${prefix_clsubset[$i]} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs[$i]} --file_cytokines_cutoffs ${file_cytokines_cutoffs[$i]} --clsubset ${clsubset[$i]}
done


# --------------------------------------------------
# Analysis of CD69+ and CD69- cells based on bimatrix
# in CK_2016-06-23_02_CD4_merging2 and CK_2016-06-23_02_CD8_merging2
# Use Analysis block 5
# --------------------------------------------------

### Analysis of cytokines for CD69+

DATA=23
PANEL=2
file_metadata="metadata_23_02"

prefix_pca="pca1_"
prefix_clust="cl20_"

data_dir=('CK_2016-06-23_02_CD4_merging2' 'CK_2016-06-23_02_CD8_merging2')
prefix_data=('23CD4_' '23CD8_')
file_panel=('panel2CD4.xlsx' 'panel2CD8.xlsx')
prefix_panel=('02CD4_' '02CD8_')

prefix_merging=('merging2_' 'merging2_')
clsubset=("c('CM','EM','TE','TM')" "c('CM','EM','TE','TM')")
prefix_clsubset=('Tmem_' 'Tmem_')
file_cytokines_cutoffs=('panel2CD4_23_cytokines_CM_CD69.xlsx' 'panel2CD8_23_cytokines_CM_CD69.xlsx')
marker="CD69"
prefix_cytokines_cutoffs=('cytCM_raw2_cd69_' 'cytCM_raw2_cd69_')

outdir="070_cd69_bimatrix"

som_dim=(5 5)
nmetaclusts=(25 25)

for i in 0 1
do
  ./Analysis_block_5_pd1_bimatrix_0_bimatrix.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --pd1_bimatrix ${cd69_bimatrix} --PANELS ${PANELS} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging[$i]} --prefix_clsubset ${prefix_clsubset[$i]} --prefix_cytokines_cutoffs ${prefix_cytokines_cutoffs[$i]} --file_cytokines_cutoffs ${file_cytokines_cutoffs[$i]} --clsubset ${clsubset[$i]} --marker ${marker} --outdir ${outdir}

  ./Analysis_block_6_cytokines_bimatrix_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cytokines_bimatrix_main ${cd69_bimatrix_main} --PANELS ${PANELS} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging[$i]} --prefix_clsubset ${prefix_clsubset[$i]} --prefix_cytokines_cutoffs "${prefix_cytokines_cutoffs[$i]}positive_"  --som_dim ${som_dim[$i]} --nmetaclusts ${nmetaclusts[$i]} --outdir ${outdir}

done



###############################################################################################################
# Analysis of CK_2016-06-23_03 data
# Use Analysis block 1
###############################################################################################################

DATA=23
PANEL=3
data_dir="CK_2016-06-23_03"

file_panel="panel3.xlsx"
file_metadata="metadata_23_03"

pca_score_cutoff=0.9
rand_seed_consensus=1234
nmetaclusts=20

prefix_data="23_"
prefix_panel="03_"
prefix_pca="pca1_"
prefix_clust="cl20_"

./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --prepare_metadata ${prepare_metadata} --data_normalization ${data_normalization} --pcascores ${pcascores} --select_observables ${select_observables} --flowsom ${flowsom} --flowsom_validation ${flowsom_validation} --heatmaps ${heatmaps} --runtsne ${runtsne} --plottsne ${plottsne} --plottsne_expr ${plottsne_expr} --frequencies ${frequencies} --expression false --METADATA ${METADATA} --path_fun_prepare_metadata ${path_fun_prepare_metadata} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin} --path_fun_formulas ${path_fun_formulas} --outdir_fun_formulas ${outdir_fun_formulas}


for i in 6
do
  nmetaclusts=$i
  prefix_clust="cl${i}_"

  ./Analysis_block_1_main.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --prepare_metadata false --data_normalization false --pcascores false --select_observables false --flowsom ${flowsom} --flowsom_validation false --heatmaps ${heatmaps} --runtsne false --plottsne ${plottsne} --plottsne_expr false --frequencies ${frequencies} --expression false --METADATA ${METADATA} --path_fun_prepare_metadata ${path_fun_prepare_metadata} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --pca_score_cutoff ${pca_score_cutoff} --rand_seed_consensus ${rand_seed_consensus} --nmetaclusts ${nmetaclusts} --tsne_pmin ${tsne_pmin} --path_fun_formulas ${path_fun_formulas} --outdir_fun_formulas ${outdir_fun_formulas}

done


# --------------------------------------------------
# Analysis of CK_2016-06-23_03 merging2
# Use Analysis block 2
# --------------------------------------------------

prefix_clust="cl20_"

file_merging="010_helpfiles/${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_merging4.xlsx"
prefix_merging="merging4_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging} --path_fun_formulas ${path_fun_formulas} --outdir_fun_formulas ${outdir_fun_formulas}







#
