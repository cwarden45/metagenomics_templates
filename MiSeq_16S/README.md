### Reference Preparation ###

a) RDPclassifier - already trained against RDP reference

### Order to Run Scripts ###

1) preprocess_RDPclassifier.py

2) merged_length_QC_stats.R

3) create_count_table_RDPclassifier.R

4) qc.R

5) differential_abundance.R


### Dependencies (some optional) ###

Biopython: http://biopython.org/wiki/Biopython

PEAR: http://sco.h-its.org/exelixis/web/software/pear/

RDPclassifier: https://sourceforge.net/projects/rdp-classifier/

### Parameter Values ###
| Parameter | Value|
|---|---|
|comp_name | Name of differential abundance comparison (used to name output file)
|plot_groups | Names of columns in *sample_description_file* to be plotted in QC and differential abundance plots.  Use commas to plot multiple groups|
|deg_groups|Names of columns in *sample_description_file* to be plotted in QC and differential abundance plots.  Use commas to include multiple variables|
|treatment_group|Treatment group for primary variable; enter *continuous* for a continuous variable.|
|Classification_Folder|Folder for per-sample count and abundance values|
|Reads_Folder|Path to Raw Reads|
|RDPclassifier_Jar|Full Path to RDPclassifier .jar File|
|Java_Mem|Memory allocation for Java (RDPclassifier)|
|Threads|Number of Reads (PEAR)|
|Classifier|Method to assign genus-level classifications.  Can be *RDPclassifier*, *mothur*, or *BWA*|
|pvalue_method|Method to Calculate P-value.  Can be *edgeR*, *limma-voom*, *DESeq2*,|
|fdr_method|Method to Calculate FDR.  Can be *BH* (Benjamini and Hochberg),*q-value*, or *q-lfdr*|
|sample_description_file|Name of Sample Description File|
|total_counts_file|Name of File to Contain Total Read Counts|
|classified_stats_file|Name of File to Contain Classification Rate Information|
|cluster_distance| Distance metric for dendrogram.  Can be *Euclidean* or *Pearson_Dissimilarity*|
|abundance_file|Name of File to Contain Percent Abundance Values|
|counts_file|Name of File to Contain Read Counts Per Genus|
|abundance_cutoff|Minimum reliable abundance level|
|minimum_fraction_expressed|Minimum fraction of samples with expression above *abundance_cutoff*|
|pvalue_cutoff|Maximum p-value to consider a genus as having differential abundance|
|fdr_cutoff|Maximum FDR to consider a genus as having differential abundance|
