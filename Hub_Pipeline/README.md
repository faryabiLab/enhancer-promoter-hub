# Spatial Enhancer-Promoter Hub Identification
Pipeline used to 1) filter and annotate Hi-C or SMC1 HiChIP contact frequency data with 
enhancers and promoters and 2) cluster this data to identify enhancer-promoter hubs
(and differential hubs), which are networks of regulatory elements that spatially 
interact within the nucleus.

## Hi-C Interaction Filtering 
Uses Hi-C chromatin contact frequency, H3K27ac ChIP-seq peaks, ATAC-seq peaks, TSS, and 
RNA-seq gene expression to identify spatial interactions between enhancers (i.e. 
accessible H3K27ac peaks) and promoters (i.e. actively transcribed, accessible TSSes). 
Spatial interactions between valid regulatory elements are then assigned a normalized 
contact frequency score and filtered to yield a datatable of valid spatial interactions
that can be used downstream to identify hubs. Example data is provided to illustrate 
input data formatting requirements. 

## SMC1 HiChIP Interaction Filtering
Uses SMC1 HiChIP chromatin contact frequency, H3K27ac ChIP-seq peaks, TSS, and 
RNA-seq gene expression to identify spatial interactions between enhancers (i.e.
H3K27ac peaks) and promoters (i.e. actively transcribed TSSes). Spatial interactions 
between valid regulatory elements are then assigned a normalized contact frequency 
score and filtered to yield a datatable of valid spatial interactions
that can be used downstream to identify hubs. Example data is provided to illustrate 
input data formatting requirements. 

## Hub Pipeline for Spatial Enahncer-Promoter Hub Detection 
Once valid spatial interactions between regulatory elements are determined, hub_pipeline.sh
can be run to perform hierarhical spectral clustering in order to construct enhancer-promoter
hubs, which can be subsequently categorized and filtered by their number of interactions
and regulatory element nodes as well as the expressed genes contained within their intervals. 
This script also contains the option to compare hubs across two separate models on the basis of 
their genomic interval and interaction counts in order to identify differential hubs in silico. 
Finally, calculate_hyperconnected_hubs.R can be used to identify hyperinteracting (i.e. hyperconnected) hubs from a file of valid hubs within a cell model on the basis of interaction count. 
