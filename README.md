# scRNA-seq-for-UCB-HSPC-iNK
Title: Large-scale generation of iNK and CAR iNK cells from CD34+ hematopoietic stem and progenitor cells for adoptive immunotherapy  
bioRxiv: [https://www.biorxiv.org/content/10.1101/2025.01.07.631650v1](https://www.biorxiv.org/content/10.1101/2024.07.30.605741v1)  

Summary:  
Raw data of droplet-based scRNA-seq （10 × Genomics） can be downloaded and processed by CellRanger software package (version 6.0.1), then aggregated by the ‘aggr’ function and subjected to Seurat (version 3.2.3) for further analysis.  

Raw data of scRNA-seq:  
The scRNA-seq data of UCB-NK, ESC-iNK, iPSC-iNK, iNK (UCB HSPC-iNK) and CD19 CAR-iNK (UCB HSPC-CD19 CAR-iNK) have been deposited in the GSA public database. The accession number of ESC-iNK and iPSC-iNK is HRA001609 (published in "Lateral plate mesoderm cell-based organoid system for NK cell regeneration from human pluripotent stem cells "([https://pubmed.ncbi.nlm.nih.gov/36344493/](https://pubmed.ncbi.nlm.nih.gov/36344493/)). The accession number of UCB-NK, iNK and CD19-CAR-iNK is HRA007978 (unpublished in this work). 

# scRNA seq analysis pipeline
## Projection of UCB-NK (Naural UCB NK), iNK (induced NK from human UCB CD34+ HSPC) and CD19 CAR-iNK (CD19 CAR-iNK from human UCB CD34+ HSPC)
    01_projection_UCBNK_and_iNK.R

## Identification of DEGS and GO terms enriched between (ESC-iNK vs iNK, iPSC-iNK VS iNK, UCB-NK VS iNK, CD19 CAR-iNK VS iNK)
    02_DEanalysis_GOenrichment.R

## Visulization of gene expression
    03_genes_expression.R
    
## Projection of UCB-NK, ESC-iNK (hESC derived iNK), iPSC-iNK (iPSC derived iNK), iNK and CD19 CAR-iNK
    04_projection_UCBNK_and_iNK_ESCiPS_iNK.R
