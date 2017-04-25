RhesusLncRNA

This is the readme file.

This repository contains the data and code for the literature "Annotation, and cluster analysis of spatiotemporal- and sex-related lncRNA expression in primate brain". To have an overview of the data in the literature, we uploaded the modules of lncRNAs and mRNAs by WGCNA analysis, including the detail information of gene expression level in each module. The co-expression between lncRNAs and mRNAs were also uploaded.

At the same time, The code that was used to generate the result was also uploaded. All code was written in the R or Perl programming language. WGCNA program was written in R, and co-expression program was written in Perl.

The following description is the detail information for eac file and folder.

#####
The "LncRNA_modules" foler
This folder contains the lncRNA modules by WGCNA method descriped in the study. The detail description for each file is as following:

All_DELnc_WGCNA_modules.xlsx: This file contains the lncRNA expression information and the corresponding module number and color.

DElnc_WGCNA_module_eigengenes.xlsx: This file contains the eigengenes value for each module.

Dendrogram_and_module_colors.pdf: This file is the dendrogram heatmap for all lncRNAs.

Eigengene_dendrogram_heatmap.pdf: This file is the dendrogram heatmap of eigengenes for all modules.

ME*bar_plot.pdf: These files are the bar plot of eigengenes for all lncRNA modules.


#####
The "mRNA_modules" folder
This folder contains the similar information to the "LncRNA_modules" folder, except that it contains the mRNA modules information.


#####
The "LncRNA_homology_analysis" folder
This folder contains the homology analysis result of rhesus lncRNA with other species. The detail description for each file is as following:

LncRNA_homology_result.xlsx: This file contains the blast result of rhesus lncRNAs with other species. The other species contains mouse, human, rhesus and gorilla from noncode database (http://www.noncode.org/). The file format is the tab separated blast result format.


#####
The "mRNA-lncRNA_co-expression_network" folder
This folder contains the co-expression analysis result between lncRNAs and mRNAs. The detail description for each file is as following:
mRNA_LncRNA_network_cor0.7_p0.01_filter_Clu.txt: This file is the filtered co-expression results, which contains the Pearson correlation coefficient and the p-value. The other two columns are the mRNA and lncRNA module number, respectively.


#####
The "Code" folder
This folder contains the perl and R script that could be called to generate the analysis result. The detail description for each file is as following:

cor.pl: This file is the main correlation coefficient calculation script, and if would call the "cor.r" program in this folder. There is a description in the program when you call this script.

cor.r: This is the correlation coefficient calculation script that is called by "cor.pl".

filter_target.pl: This file is the result filtering script that following the "cor.pl" program. There is a description in the program when you call this script.

WGCNA.r: This file is the co-expression analysis code using the WGCNA module. There is a description in the program when you call this script.


#####
All_sample_LncRNA_exp_RPKM.xlsx: This file is the expression level of all lncRNAs.


#####
All_sample_merged_RPKM_ed.xlsx: This file is the expression level of all genes.
