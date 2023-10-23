#!/bin/bash
#### --------------------------------------------------------------------------------------------------------------------------//
#### Overview: This script contains the complete pipeline for separately identifying hubs within 2 conditions AND for 
# creating a list of union hubs between the 2 conditions (i.e. control & case) for use as an in silico differential hub screen.                                                                                                                                 
#
#### Inputs: Filtered .csv datatables of PP/EE/EP interactions in the control and case conditions (chrom_start1_end1,chrom_start2_end2,contactFrequency), 
# .csv files of annotated enhancer/promoters in each condition (chrom_start_end,"enhancer/promoter"), BED4 file of all genes for given genome build (e.g. hg19) 
# in format (chrom \t start \t end \t gene_symbol \n), BED4 file of expressed genes in either case and/or control conditions, 
# and paths/names of relevant input/output dirs/files
#
#### Dependencies: cluster-tree (https://github.com/GregorySchwartz/hierarchical-spectral-clustering), bedtools suite,
# cluster_analysis_GW.py (plots.py, igraph, altair), plots.py (pandas, altair, mygene, igraph, functools), make_contiguous_hub_intervals.py,
# make_union_hublist_comprehensive_node.py, hub_interval_to_gene_node.py
#
#### Intermediate Outputs: Multiple lists of hubs. 
#
#### Final Outputs: 1) 2 lists of hubs in case & control conditions (respectively) with interaction and enhancer/promoter count attributes
# with columns "chrom /t start /t end /t hub_id /t interaction_count /t enhancer_promoter_count" 
# 2) 1 list of the union of hubs in case & control conditions with interaction and enhancer/promoter count attributes and gene annotations
# with columns "chrom /t start /t end /t control_hub_id /t control_interaction_count /t control_enhancer_promoter_count /t case_hub_id /t case_interaction_count /t case_enhancer_promoter_count /t log2_FC_interaction_count /t gene_annotations /t de_novo_gain_lost_status /t union_status"
#### -------------------------------------------------------------------------------------------------------------------------//

#### Inputs: Paths to case and control .csv tables of filtered PP/EE/EP interactions and .csv tables of labeled anchors (enh/pro),
# path to input directory, path of final desired hub unionlist output file name, and path to .bed file(s) containing
# the chromosomal intervals of (expressed) genes within either case or control conditions (optional)
INPUTDIR=/mnt/data0/brent/analysis/idea32_230824_final_hub_scripts_submission/230825_Finalized_Hub_Scripts/Hub_Pipeline
DATADIR=/mnt/data0/brent/analysis/idea32_230824_final_hub_scripts_submission/230825_Finalized_Hub_Scripts/Hub_Pipeline/example_input_data
CONTROL_INPUT=${DATADIR}/Control_ct_filt.csv
CASE_INPUT=${DATADIR}/Case_ct_filt.csv
CONTROL_LABEL=${DATADIR}/Control_label.csv
CASE_LABEL=${DATADIR}/Case_label.csv
OUT_NAME=${DATADIR}/REC1_hub_IBR_RES_unionlist.txt
ALL_GENES=${DATADIR}/hg19_ucsc_gene_name_interval_unique.bed
EXPRESSED_GENES_CONTROL=${DATADIR}/REC1_DMSO_24h_exp_genes.bed
EXPRESSED_GENES_CASE=${DATADIR}/REC1_RES_24h_exp_genes.bed
EXPRESSED_GENES_UNION=${DATADIR}/REC1_DMSO_RES_24h_exp_genes_uniq2.bed

#### Make folder to store intermediate output as working directory
WORKDIR=${DATADIR}/intermediate_files
mkdir ${WORKDIR}
cd ${WORKDIR}
echo "intermediate_files created as working directory, beginning clustering..."

#### [Optional] Move to conda environment with python that contains all dependencies
CONDA_PREFIX=/mnt/data0/apps/anaconda/anaconda2
source ${CONDA_PREFIX}/etc/profile.d/conda.sh
conda activate brent_py3

#### Call cluster-tree function to perform hierarchical spectral clustering on control & case EE/PP/EP interactions separately
cat $CONTROL_INPUT | ~/.local/bin/cluster-tree -c Dense -n 2 -s -o Control_hub_output_trees/tree.json > Control_hub_ct_output.csv
cat $CASE_INPUT | ~/.local/bin/cluster-tree -c Dense -n 2 -s -o Case_hub_output_trees/tree.json > Case_hub_ct_output.csv
echo "cluster-tree successfully called on Control and Case input files to generate lists of hubs"

#### [Optional/Adjustable] Modify output to eliminate subtrees to visualize consolidated hubs (i.e. merge all subtrees into parent tree)
# This step can be removed/altered if a more stringent definition of hubs is desired
awk -F/ 'BEGIN{print "item,cluster"} NR>1 {print $1"/1"}' ${WORKDIR}/Control_hub_ct_output.csv > ${WORKDIR}/Control_hsc.csv
awk -F/ 'BEGIN{print "item,cluster"} NR>1 {print $1"/1"}' ${WORKDIR}/Case_hub_ct_output.csv > ${WORKDIR}/Case_hsc.csv

#### Call python helper program to calculate the number of interactions, enhancers/promoters, etc. in each hub [longest step]
GW_FUNC=${INPUTDIR}/cluster_analysis_GW.py
python $GW_FUNC ${WORKDIR}/Control_hsc.csv $CONTROL_INPUT $CONTROL_LABEL ${WORKDIR}/Hub_properties_control/
python $GW_FUNC ${WORKDIR}/Case_hsc.csv $CASE_INPUT $CASE_LABEL ${WORKDIR}/Hub_properties_case/
echo "Control and Case hubs are successfully characterized by their connectivity and network properties"

#### From hub property output datatable, extract interaction count, enhancer/promoter count, and list of coordinates of all enh/pros in each hub into new table
awk -F, 'NR>1 {print $2"\t"$4"\t"$3"\t"$18}' ${WORKDIR}/Hub_properties_control/Control_hsc.csv_cluster_plot_df_reconnected.csv | sort -k2nr > ${WORKDIR}/Control_reconnected_EDGE.tsv
awk -F, 'NR>1 {print $2"\t"$4"\t"$3"\t"$18}' ${WORKDIR}/Hub_properties_case/Case_hsc.csv_cluster_plot_df_reconnected.csv | sort -k2nr > ${WORKDIR}/Case_reconnected_EDGE.tsv

#### Characterize hubs based on the largest genomic interval covered by their enhancer/promoters
INTERVAL_FUNC=${INPUTDIR}/make_contiguous_hub_intervals.py
python $INTERVAL_FUNC ${WORKDIR}/Control_reconnected_EDGE.tsv ${WORKDIR}/Control_hub_by_interval.txt
python $INTERVAL_FUNC ${WORKDIR}/Case_reconnected_EDGE.tsv ${WORKDIR}/Case_hub_by_interval.txt
echo "Control and Case hubs successfully categorized by the contiguous genomic interval covered by their enhancer/promoters"

#### [Optional] Create .bed file of hub genomic intervals to allow for visualization in (IGV) gene browser to help evaluate
# if contact frequency cutoff for filtering input EE/PP/EP interactions was valid
awk 'NR>1 {print $1"\t"$2"\t"$3}' ${WORKDIR}/Control_hub_by_interval.txt > ${DATADIR}/Control_hub_by_interval_IGV.bed
awk 'NR>1 {print $1"\t"$2"\t"$3}' ${WORKDIR}/Case_hub_by_interval.txt > ${DATDIR}/Case_hub_by_interval_IGV.bed

#### Final output 1: 2 separate lists of control and case hubs (respectively) that are categorized by genomic interval, interaction count, and enhancer/promoter count
#### ** NOTE THAT THESE FILES MUST EACH BE FILTERED PRIOR TO DOWNSTREAM ANALYSIS IN ORDER TO ISOLATE TRUE HUBS WITH INTERACTION COUNT > 5 **
#### Remove header and sort hub interval lists by interaction count to create control and case hub lists
# that can be used for all downstream analyses examining these conditions in isolation from one another.
# Sorting hubs in this manner also faciliates pairing of the Control and Case hubs in subsequent step to create a union hub list
awk 'NR>1 {print $0}' ${WORKDIR}/Control_hub_by_interval.txt | sort -k5nr > ${INPUTDIR}/example_input_data/Control_hub_by_interval_sort.txt
awk 'NR>1 {print $0}' ${WORKDIR}/Case_hub_by_interval.txt | sort -k5nr > ${INPUTDIR}/example_input_data/Case_hub_by_interval_sort.txt

#### Filter hub lists such that only true hubs with interaction count > 5 are kept; do NOT use this as input to differential hub calling
awk '{if($5 > 5){print $0}}' ${WORKDIR}/Control_hub_by_interval_sort.txt > ${DATADIR}/Control_hub_by_interval_sort_filt.txt
awk '{if($5 > 5){print $0}}' ${WORKDIR}/Case_hub_by_interval_sort.txt > ${DATADIR}/Case_hub_by_interval_sort_filt.txt

#### (Optional addition for final output 1): Annotate control and case hub lists with expressed genes in each condition
ANNOTATE_FUNC=${INPUTDIR}/hub_interval_to_gene_node.py
# Annotate union hubs with lists of genes that are expressed in either case or control conditions (from RNA-seq data) that are contained within hub genomic intervals
python $ANNOTATE_FUNC ${DATADIR}/Control_hub_by_interval_sort.txt $EXPRESSED_GENES_CONTROL ${DATADIR}/Control_hub_by_interval_sort_exp_gene.txt single
python $ANNOTATE_FUNC ${DATADIR}/Case_hub_by_interval_sort.txt $EXPRESSED_GENES_CASE ${DATADIR}/Case_hub_by_interval_sort_exp_gene.txt single

#### Generate union list of hubs between control and case hub lists, eliminate hubs with less than 6 edges in both 
# control AND case conditions from union list, and sort output by chromosomal interval
UNION_FUNC=${INPUTDIR}/make_union_hublist_comprehensive_node.py
NO="no_gene.txt"
python $UNION_FUNC ${WORKDIR}/Control_hub_by_interval_sort.txt ${WORKDIR}/Case_hub_by_interval_sort.txt $OUT_NAME$NO
sort -k5nr -k8nr $OUT_NAME$NO > ${OUT_NAME}.sort.${NO}

#### Final output 2: One unionlist of all hubs from control and case conditions mapped onto one another and categorized by interaction count, enhancer/promoter count,
# and (expressed) genes as well as other attributes in each condition
ANNOTATE_FUNC=${INPUTDIR}/hub_interval_to_gene_node.py
# Annotate union hubs with lists of genes that are expressed in either case or control conditions (from RNA-seq data) that are contained within hub genomic intervals
python $ANNOTATE_FUNC ${OUT_NAME}.sort.${NO} $EXPRESSED_GENES_UNION ${OUT_NAME}.expgene.txt union

