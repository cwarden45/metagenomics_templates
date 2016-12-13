param.table = read.table("parameters.txt", header=T, sep="\t")
sample.description.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
classifier = as.character(param.table$Value[param.table$Parameter == "Classifier"])
classification.folder = as.character(param.table$Value[param.table$Parameter == "Classification_Folder"])
read.folder = as.character(param.table$Value[param.table$Parameter == "Reads_Folder"])
min.cycles = as.character(param.table$Value[param.table$Parameter == "Min_CCS_Cycles"])
plot.groups = unlist(strsplit(as.character(param.table$Value[param.table$Parameter == "plot_groups"]), split=","))

fixed.color.palatte = c("green","orange","purple","cyan","pink","maroon","yellow","grey","red","blue","black","darkgreen","thistle1","tan","orchid1",colors())

xrange = c(1300,1700)
yrange = c(0,0.3)

sample.table = read.table(sample.description.file, sep="\t", header=T)
sampleID = as.character(sample.table$sampleID)
groupID = as.character(sample.table$Group)

length.upper5= rep(NA, length(sampleID))
length.upper.quantile = rep(NA, length(sampleID))
length.med = rep(NA, length(sampleID))
length.mean = rep(NA, length(sampleID))
length.lower.quantile = rep(NA, length(sampleID))
length.lower5 = rep(NA, length(sampleID))

length.files = file.path(read.folder,sampleID)
if ((classifier == "RDPclassifier")|(classifier == "BWA")){
	print(classifier)
	length.files = paste(length.files, ".ccs.",min.cycles,"x.read_length",sep="")
} else if(classifier == "mothur"){
	print(classifier)
	length.files = paste(length.files, ".ccs.",min.cycles,"x.good.read_length",sep="")
} else{
	stop(paste("Need to write code for classifier: ",classifier,sep=""))
}
		
group.levels = levels(as.factor(groupID))
num.sample.types = length(group.levels)
color.palette = fixed.color.palatte[1:num.sample.types]

labelColors = rep("black",times=length(length.files))
for (i in 1:length(groupID)){
	labelColors[groupID == as.character(group.levels[i])] = color.palette[i]
}#end for (i in 1:length(target.groups))
		
for (group in plot.groups){
	dist.file  = paste("filtered_CCS_length_dist_by_",group,".png",sep="")
	png(dist.file)
	for (i in 1:length(length.files)){
		length.table = read.table(length.files[i],sep="\t", head=T)
		length.quantiles = quantile(length.table$Length, c(0.05,0.25,0.5,0.75, 0.95))
		length.lower5[i] = length.quantiles[1]
		length.lower.quantile[i] = length.quantiles[2]
		length.med[i] = length.quantiles[3]
		length.upper.quantile[i] = length.quantiles[4]
		length.upper5[i] = length.quantiles[5]
		length.mean[i] = mean(length.quantiles)
				
		den <- density(length.table$Length, na.rm=T, from=xrange[1], to=xrange[2])
		if(i == 1){
			plot(den$x, den$y, type="l", xlab ="Filtered CCS Length",
					ylab = "Density",col=labelColors[i], ylim = yrange)
		}else{
			lines(den$x, den$y, type = "l", col=labelColors[i], new=FALSE)
		}#end else
	}#end for (i in 1:length(target.files))
	legend("topleft",legend=group.levels, col=color.palette, lwd=2)
	dev.off()
			
}#end for (group in plot.groups)

output.table = data.frame(sample=sampleID,source=length.files,
							length.lower5=length.lower5, length.lower.quantile=length.lower.quantile, length.med=length.med,
							length.mean=length.mean, length.upper.quantile=length.upper.quantile, length.upper5=length.upper5)
write.table(output.table,"filtered_CCS_length_quantiles.txt",sep="\t", row.names=F, quote=F)