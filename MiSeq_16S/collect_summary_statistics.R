sample.file = "sample_description.txt"
fastqc.file = "FastQC_total_reads.txt"
length.count.file = "total_read_counts.txt"
classified.stat.file = "classification_stats.txt"

output.file = "combined_statistics.txt"

sample.table = read.table(sample.file, head=T, sep="\t")

fastqc.table = read.table(fastqc.file, head=T, sep="\t")
fastqc.table = fastqc.table[match(sample.table$sampleID,fastqc.table$SampleID),]

length.table = read.table(length.count.file, head=T, sep="\t")
totalID = as.character(length.table$Sample)
#totalID = gsub("_L001_R1_001","",totalID)
length.table = length.table[match(sample.table$sampleID,totalID),]

class.table = read.table(classified.stat.file, head=T, sep="\t")
classID = as.character(class.table$SampleID)
#classID = gsub("_L001_R1_001","",classID)
class.table = class.table[match(sample.table$sampleID,classID),]

output.table = data.frame(Sample=sample.table$sampleID,
							userID = sample.table$userID,
							total.reads = fastqc.table$TotalReads,
							mothur.merged.qc.count = length.table$Total.Reads,
							length.pass.rate = paste(length.table$Percent.Length.Filtered,"%",sep=""),
							percent.classified = class.table$Percent.Classified,
							classifed.reads = class.table$Classified.Reads)
write.table(output.table, output.file, quote=F, sep="\t", row.names=F)