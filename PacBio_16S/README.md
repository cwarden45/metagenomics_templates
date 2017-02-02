### Reference Preparation ###

**a) RDPclassifier** - already trained against RDP reference

**b) mothur** - indices are automaticaly produced, but you need to provide mothur-formatted reference sequences.  See the [MiSeq SOP](http://www.mothur.org/wiki/MiSeq_SOP) for more details, which provides links to download SILVA and RDP formatted databases.

**NOTE**: Template assumes you use SILVA sequences with RDP taxonomy

**NOTE**: The templates provided work with mothur version 1.38.0.  Some earlier versions of mothur provide classifications in a different format, and are therefore incompatible with these templates.  This particular version mothur is available in this [Docker image](https://hub.docker.com/r/cwarden45/metagenomics/).

**c) BWA** - prepare reference sequences using functions to re-train RDP classifier, then prepare BWA index.  See `RDP_BWA_filter_seqs.py` for example.  FASTA sequences with taxonomy information in header are available from [RDPclassifier sourceforge page](https://sourceforge.net/projects/rdp-classifier/files/RDP_Classifier_TrainingData/).

**d) BLAST** - FASTA reference should be prepared using 'makeblastdb'.  Testing performed using RDP sequences (same as BWA), using version of BLAST+ (v2.4.0) available in this [Docker image](https://hub.docker.com/r/cwarden45/metagenomics/) (not legacy BLAST v2.2.22 used by QIIME)).

### Order to Run Scripts ###

1) `full_run_bax2bam.py` + `create_CCS_and_FASTQ.py` (or `run_pbbarcode_template.py` + `reformat_barcoded_CCS_FASTQ.py`)

PacBio dependencies available in this [Docker image](https://hub.docker.com/r/cwarden45/general-pacbio/).

2) `preprocess_RDPclassifier_or_BWA.py` or `preprocess_mothur.py`

3) `length_QC_stats.R` (one way to determine minimum and maximium merged length)

-If using full-length sequences, I would recommend setting Min.Length=1450 and Max.Length=1550

4) `run_classifier.py`

5) `create_count_table.py`

6) `qc.R`

7) `differential_abundance.R`


### Dependencies (some optional) ###

pbccs: https://github.com/PacificBiosciences/pbccs

ConsensusTools: https://github.com/PacificBiosciences/SMRT-Analysis/wiki/ConsensusTools-v2.3.0-Documentation

pbbarcode: https://github.com/PacificBiosciences/pbbarcode

Biopython: http://biopython.org/wiki/Biopython

RDPclassifier: https://sourceforge.net/projects/rdp-classifier/

mothur: http://www.mothur.org/

BWA: http://bio-bwa.sourceforge.net/

BLAST+: https://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/

samtools: http://samtools.sourceforge.net/

metagenomeSeq: https://bioconductor.org/packages/release/bioc/html/metagenomeSeq.html

limma: http://bioconductor.org/packages/release/bioc/html/limma.html

heatmap.3: https://github.com/obigriffith/biostar-tutorials/blob/master/Heatmaps/heatmap.3.R

heatmap.3 example: https://www.biostars.org/p/18211/

### Parameter Values ###
| Parameter | Value|
|---|---|
|comp_name | Name of differential abundance comparison (used to name output file)
|plot_groups | Names of columns in *sample_description_file* to be plotted in QC and differential abundance plots.  Use commas to plot multiple groups|
|deg_groups|Names of columns in *sample_description_file* to be plotted in QC and differential abundance plots.  Use commas to include multiple variables|
|treatment_group|Treatment group for primary variable; enter *continuous* for a continuous variable.|
|Classification_Folder|Folder for per-sample count and abundance values|
|Reads_Folder|Path to CCS FASTQ Files|
|Raw_CCS_Folder|Path to CCS Unaligned BAM|
|Min_CCS_Cycles|Minimum cycles to define CCS read (recommend at least 5, default in `ccs` is 3)|
|Raw_Code_PC|Path to output folder for most results|
|Result_Folder|Path to output folder for selected, final results|
|Java_Mem|Memory allocation for Java (used by RDPclassifier)|
|Threads|Number of Reads (used by PEAR, BWA, mothur)|
|Classifier|Method to assign genus-level classifications.  Can be *RDPclassifier*, *mothur*, *BWA*, or *BLAST*|
|RDPclassifier_Jar|Full Path to RDPclassifier .jar File|
|mothur_ref|Path to mothur-formatted reference sequence|
|mothur_tax|Path to mothur-formatted taxonomy file|
|BWA_Ref|Path to indexed BWA reference (also used for BLAST)|
|pvalue_method|Method to Calculate P-value.  Can be *limma-ab*, *limma-counts*, or *metagenomeSeq* (*limma-ab* uses abundance percentages, *limma-counts* uses limma-voom)|
|fdr_method|Method to Calculate FDR.  Can be *BH* (Benjamini and Hochberg),*q-value*, or *q-lfdr*|
|sample_description_file|Name of Sample Description File|
|total_counts_file|Name of File to Contain Total Read Counts|
|classified_stats_file|Name of File to Contain Classification Rate Information|
|cluster_distance| Distance metric for dendrogram.  Can be *Euclidean*, *Pearson_Dissimilarity*, *Spearman_Dissimilarity*, or *Kendall_Dissimilarity*|
|abundance_file|Name of File to Contain Percent Abundance Values|
|counts_file|Name of File to Contain Read Counts Per Genus|
|abundance_cutoff|Maximum abundance must be above this level to be considered for differential abundance analysis|
|rounding_factor|Rounding factor for abundance fold-change value (should be below `abundance_cutoff` if you want to detect genera whose maximum abundance is near cutoff)|
|fc_cutoff|Minimum linear fold-change difference between groups to be considered differentially expressed|
|pvalue_cutoff|Maximum p-value to consider a genus as having differential abundance|
|fdr_cutoff|Maximum FDR to consider a genus as having differential abundance|
|sec_pvalue_cutoff|If comparing two gene lists, p-value threshold for list you want to filter out|
|sec_fdr_cutoff|If comparing two gene lists, FDR threshold for list you want to filter out|
|interaction| Method for comparing an interaction of two variables.  Can be *model*, *filter-overlap*, or *no*|
|secondary_trt| If comparing two gene lists, this is treatment group for the list that you want to filter out; enter *continuous* for a continuous variable and a correlation will be provided instead of a fold-change value (also converts second variable from factor to numeric, even if interaction is set to *no*)|
