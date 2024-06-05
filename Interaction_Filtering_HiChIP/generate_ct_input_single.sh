#!/bin/bash
#### ---------------------------------------------------------------------------------------------//
#### This is a helper script to "make_HiChIP_clique_input_single.R".
#### Overview: This script organizes a dataframe of EE/PP/EP filtered interactions for condition 1 
# into a .csv formatted file for input into the hierarchical spectral clustering hub-calling pipeline. It also generates a .csv file
# containing annotations for each anchor as enhancer or promoter as well as an [optional] .bedpe file of interactions 
# for visualization with IGV or another platform. 
#
#### Args (passed from command line call): 1) .csv file of filtered EE/PP/EP contacts (with header)
#
#### Outputs: 1 .csv file containing filtered EE/PP/EP interactions with weights for direct input into hub-calling pipeline 
# in format (chr_A_start1_stop1,chrA_start2_stop2,score/weight), 1 .csv file containing anchor coordinates and annotations 
# in format (chr_A_start1_stop1,"enhancer/promoter"), and 1 .bedpe formatted file with filtered EE/PP/EP interactions
#### ---------------------------------------------------------------------------------------------//

#### Input file location and name
CONTROL_CSV=$1

#### Generate cluster-tree input file
awk -F, 'NR>1 {print $1","$2","$3}' $1 > Control_ct_filt.csv

#### Generate label file
awk -F, 'BEGIN{print "item,cluster"} NR>1 {print $1","$4"\n"$2","$5}' $1 > Control_label.csv

#### [Optional] Generate .bedpe file for manual input into IGV to visualize contact arches
awk -F, 'NR>1 {print $1"\t"$2}' $1 > Control_arches.int
awk -F_ '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' Control_arches.int > Control_arches.bedpe

