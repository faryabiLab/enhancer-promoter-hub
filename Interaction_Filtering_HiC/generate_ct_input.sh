#!/bin/bash
#### ---------------------------------------------------------------------------------------------//
#### This is a helper script to "make_HiC_Clique_RandomEP_input.R".
#### Overview: This script organizes the 2 dataframes of EE/PP/EP filtered interactions for condition 1 and condition 2 into 
# 2 .csv formatted files for input into hierarchical spectral clustering hub-calling pipeline. It also generates 2 .csv files
# containing annotations for each anchor as enhancer or promoter as well as [optional] 2 .bedpe files of interactions for each condition 
# for visualization with IGV or another platform. 
#
#### Args (passed from command line call): 1) .csv file of condition 1 filtered EE/PP/EP contacts (with header), 
# 2) .csv file of condition 2 filtered EE/PP/EP contacts (with header)
#
#### Outputs: 2 .csv files containing filtered EE/PP/EP contacts and scores for each condition for direct input into hub-calling pipeline 
# in format (chr_A_start1_stop1,chrA_start2_stop2,contactFreq), 2 .csv files containing anchor coordinates and annotations for
# each condition in format (chr_A_start1_stop1,"enhancer/promoter"), and 2 .bedpe formatted files with filtered EE/PP/EP contacts for each condition
#### ---------------------------------------------------------------------------------------------//

#### Input file locations and names
CONTROL_CSV=$1
CASE_CSV=$2

# Generate cluster-tree input file for each condition
awk -F, 'NR>1 {print $1","$2","$3}' $1 > Control_ct_filt.csv
awk -F, 'NR>1 {print $1","$2","$3}' $2 > Case_ct_filt.csv

# Generate label file for each condition 
awk -F, 'BEGIN{print "item,cluster"} NR>1 {print $1","$4"\n"$2","$5}' $1 > Control_label.csv
awk -F, 'BEGIN{print "item,cluster"} NR>1 {print $1","$4"\n"$2","$5}' $2 > Case_label.csv

# [Optional] Generate .bedpe files for each condition for manual input into IGV to visualize contact arches
awk -F, 'NR>1 {print $1"\t"$2}' $1 > Control_arches.int
awk -F_ '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' Control_arches.int > Control_arches.bedpe

awk -F, 'NR>1 {print $1"\t"$2}' $2 > Case_arches.int
awk -F_ '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' Case_arches.int > Case_arches.bedpe
