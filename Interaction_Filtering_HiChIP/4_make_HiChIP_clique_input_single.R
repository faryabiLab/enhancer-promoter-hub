#### ---------------------------------------------------------------------------------------------//
#### This is the fourth of 4 scripts used to create a list of valid HiChIP EE/PP/EP interactions as input for hub calling. 
#### Overview: This script filters the EE/PP/EP interactions for the given condition (i.e. MB157 DMSO) to create a final 
# list of valid spatial interactions between putatitve regulatory elements for hub calling.
#
#### Dependencies: see library below
#
#### Required input: Dataframe containing EE/PP/EP interactions for the given condition
# 
#### Outputs: One .csv file of valid, filtered EE/PP/EP spatial interactions, which 
# can be directly used for subsequent command-line calling of hubs via clustering
# with format "enh/pro1", "enh/pro2", "interaction_weight/score", "annot1", "annot2"
#### ---------------------------------------------------------------------------------------------//

#### Load library and numeric preferences
rm(list = ls())
library("dplyr")
options(scipen=999) # Convert sci notation to decimal if necessary

#### Set workdir
setwd("test_HiChIP")
cmd <- paste("mkdir example_hub_calling")
system(cmd)
setwd("example_hub_calling")

#### Input: Load file containing EE/PP/EP spatial contacts (output of script #3)
load("test_HiChIP/MB157_SMC1_regulatory_int.rda")

#### Input: Define intermediate output file name/path
out1 <- "MB157_DMSO_PPEEEP_0.05Q.csv"

#### Clean up PPEEEP interaction dataframe
DMSO_Clique_PPEEEP = subset(DMSO.filt.clique, select=c(q_id, p_id, MB157_loop_dmso_norm, q_annot, p_annot))
colnames(DMSO_Clique_PPEEEP) <- c("item_id1", "item_id2", "loop_DMSO_norm", "item_annot1", "item_annot2")
DMSO_Clique_PPEEEP <- unique(DMSO_Clique_PPEEEP)

#### Weight all interactions equally for clustering given that they all surpass FitHiChIP significance threshold
DMSO_Clique_PPEEEP <- DMSO_Clique_PPEEEP %>% mutate(loop_DMSO_norm_2 = 1)
DMSO_Clique_PPEEEP_final <- DMSO_Clique_PPEEEP[,c(1,2,6,4,5)]

#### [Optional alternative, not illustrated] Score interactions by read count / significance for further filtering or clustering purposes

#### Save datatable to file for input into final helper bash shell script 
write.table(DMSO_Clique_PPEEEP_final, file = out1, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

#### Call helper bash shell script to reorganize .csv files to create final dataframe of EE/PP/EP interactions for input into hub-calling pipeline 
# in format (chr_A_start1_stop1,chrA_start2_stop2,score/weight) as well as enhancer/promoter annotation file for anchors 
cmd <- paste("Interaction_Filtering_HiChIP/generate_ct_input_single.sh", out1)
system(cmd)

