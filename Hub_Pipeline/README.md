## Spatial Enhancer-Promoter Hub Identification
Pipeline used to 1) filter and annotate Hi-C or SMC1 HiChIP contact frequency data with 
enhancers and promoters, 2) cluster this data to identify enhancer-promoter hubs, which are 
networks of spatially interacting regulatory elements within the nucleus, and 3) compare hubs 
across conditions.

## Hi-C Interaction Filtering 
Uses Hi-C chromatin contact frequency, H3K27ac ChIP-seq peaks, ATAC-seq peaks, TSS, and 
RNA-seq gene expression to identify spatial interactions between enhancers (i.e. 
accessible H3K27ac peaks) and promoters (i.e. actively transcribed, accessible TSSes). 
Spatial interactions between valid regulatory elements are then assigned a normalized 
contact frequency score and filtered to yield a datatable of valid spatial interactions
that can be used to identify hubs. Example data is provided to illustrate 
input data formatting requirements. 

## SMC1 HiChIP Interaction Filtering
Uses SMC1 HiChIP data, H3K27ac ChIP-seq peaks, TSS, and RNA-seq gene expression to 
identify spatial interactions between enhancers (i.e. H3K27ac peaks) and promoters 
(i.e. actively transcribed TSSes). A HiChIP interaction caller is used to detect 
significant SMC1 HiChIP interactions, and significant spatial interactions between valid 
regulatory elements are then determined and used to identify hubs. Example data is provided to 
illustrate input data formatting requirements. 

## Hub Pipeline for Spatial Enhancer-Promoter Hub Detection 
Once valid spatial interactions between regulatory elements are determined, hub_pipeline.sh
(for two conditions) or hub_pipeline_single.sh (for a single condition) can be used to cluster 
interactions and construct enhancer-promoter hubs, which are categorized by within-hub spatial 
interaction counts, regulatory element counts, and the genes broadly contained within the hubs. 
The hub_pipeline.sh script also offers the option to compare hubs across two separate models 
on the basis of their genomic overlap and interaction counts in order to identify differential hubs 
in silico. Finally, calculate_hyperconnected_hubs.R can be used to identify hyperinteracting hubs 
on the basis of interaction count. 
