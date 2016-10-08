param.table = read.table("parameters.txt", header=T, sep="\t")
sample.description.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
classifier = as.character(param.table$Value[param.table$Parameter == "Classifier"])
classification.folder = as.character(param.table$Value[param.table$Parameter == "Classification_Folder"])
read.folder = as.character(param.table$Value[param.table$Parameter == "Reads_Folder"])
plot.groups = unlist(strsplit(as.character(param.table$Value[param.table$Parameter == "plot_groups"]), split=","))

fixed.color.palatte = c("green","orange","purple","cyan","pink","maroon","yellow","grey","red","blue","black","darkgreen","thistle1","tan","orchid1",colors())

xrange = c(0,600)
yrange = c(0,1)

if ((classifier == "RDPclassifier")|(classifier == "BWA")){
	print(classifier)
	merged.read.folder = file.path(read.folder, "PEAR_Merged")
	
	sample.table = read.table(sample.description.file, sep="\t", header=T)
	fileID = as.character(sample.table$sampleID)
	shortID = as.character(sample.table$userID)
	targetID = as.character(sample.table$Target)
	groupID = as.character(sample.table$Group)
	
	length.upper5= rep(NA, length(fileID))
	length.upper.quantile = rep(NA, length(fileID))
	length.med = rep(NA, length(fileID))
	length.mean = rep(NA, length(fileID))
	length.lower.quantile = rep(NA, length(fileID))
	length.lower5 = rep(NA, length(fileID))
	
	targets = levels(as.factor(targetID))
	
	for (target in targets){
		print(target)
		
		target.files = file.path(merged.read.folder,fileID[targetID == target])
		target.files = paste(target.files, ".assembled.read_length",sep="")
		target.sample = shortID[targetID == target]
		target.groups = groupID[targetID == target]
		
		group.levels = levels(as.factor(target.groups))
		num.sample.types = length(group.levels)
		color.palette = fixed.color.palatte[1:num.sample.types]

		labelColors = rep("black",times=length(target.files))
		for (i in 1:length(target.groups)){
			labelColors[target.groups == as.character(group.levels[i])] = color.palette[i]
		}#end for (i in 1:length(target.groups))
		
		for (group in plot.groups){
			dist.file  = paste(target,"_PEAR_merged_length_dist_by_",group,".png",sep="")
			png(dist.file)
			for (i in 1:length(target.files)){
				target.table = read.table(target.files[i],sep="\t", head=T)
				length.quantiles = quantile(target.table$Length, c(0.05,0.25,0.5,0.75, 0.95))
				length.lower5[shortID == target.sample[i]] = length.quantiles[1]
				length.lower.quantile[shortID == target.sample[i]] = length.quantiles[2]
				length.med[shortID == target.sample[i]] = length.quantiles[3]
				length.upper.quantile[shortID == target.sample[i]] = length.quantiles[4]
				length.upper5[shortID == target.sample[i]] = length.quantiles[5]
				length.mean[shortID == target.sample[i]] = mean(length.quantiles)
				
				den <- density(target.table$Length, na.rm=T, from=xrange[1], to=xrange[2])
				if(i == 1){
					plot(den$x, den$y, type="l", xlab ="Merged Length", main = target,
							ylab = "Density",col=labelColors[i], ylim = yrange)
				}else{
					lines(den$x, den$y, type = "l", col=labelColors[i], new=FALSE)
				}#end else
			}#end for (i in 1:length(target.files))
			legend("topleft",legend=group.levels, col=color.palette, lwd=2)
			dev.off()
			
		}#end for (group in plot.groups)
		
	}#end for (target in targets)

	output.table = data.frame(fileID=fileID, shortID=shortID,targetID=targetID,
							length.lower5=length.lower5, length.lower.quantile=length.lower.quantile, length.med=length.med,
							length.mean=length.mean, length.upper.quantile=length.upper.quantile, length.upper5=length.upper5)
	write.table(output.table,"PEAR_merged_length_quantiles.txt",sep="\t", row.names=F, quote=F)
	
} else if(classifier == "mothur"){
	print(classifier)
	merged.read.folder = file.path(read.folder, "mothur_merged")
	
	sample.table = read.table(sample.description.file, sep="\t", header=T)
	fileID = as.character(sample.table$sampleID)
	shortID = as.character(sample.table$userID)
	targetID = as.character(sample.table$Target)
	groupID = as.character(sample.table$Group)
	
	length.upper5= rep(NA, length(fileID))
	length.upper.quantile = rep(NA, length(fileID))
	length.med = rep(NA, length(fileID))
	length.mean = rep(NA, length(fileID))
	length.lower.quantile = rep(NA, length(fileID))
	length.lower5 = rep(NA, length(fileID))
	
	targets = levels(as.factor(targetID))
	
	for (target in targets){
		print(target)
		
		target.files = file.path(merged.read.folder,fileID[targetID == target])
		target.files = paste(target.files, ".trim.contigs.good.read_length",sep="")
		target.sample = shortID[targetID == target]
		target.groups = groupID[targetID == target]
		
		group.levels = levels(as.factor(target.groups))
		num.sample.types = length(group.levels)
		color.palette = fixed.color.palatte[1:num.sample.types]

		labelColors = rep("black",times=length(target.files))
		for (i in 1:length(target.groups)){
			labelColors[target.groups == as.character(group.levels[i])] = color.palette[i]
		}#end for (i in 1:length(target.groups))
		
		for (group in plot.groups){
			dist.file  = paste(target,"_mothur_contig_length_dist_by_",group,".png",sep="")
			png(dist.file)
			for (i in 1:length(target.files)){
				target.table = read.table(target.files[i],sep="\t", head=T)
				length.quantiles = quantile(target.table$Length, c(0.05,0.25,0.5,0.75, 0.95))
				length.lower5[shortID == target.sample[i]] = length.quantiles[1]
				length.lower.quantile[shortID == target.sample[i]] = length.quantiles[2]
				length.med[shortID == target.sample[i]] = length.quantiles[3]
				length.upper.quantile[shortID == target.sample[i]] = length.quantiles[4]
				length.upper5[shortID == target.sample[i]] = length.quantiles[5]
				length.mean[shortID == target.sample[i]] = mean(length.quantiles)
				
				den <- density(target.table$Length, na.rm=T, from=xrange[1], to=xrange[2])
				if(i == 1){
					plot(den$x, den$y, type="l", xlab ="Merged Length", main = target,
							ylab = "Density",col=labelColors[i], ylim = yrange)
				}else{
					lines(den$x, den$y, type = "l", col=labelColors[i], new=FALSE)
				}#end else
			}#end for (i in 1:length(target.files))
			legend("topleft",legend=group.levels, col=color.palette, lwd=2)
			dev.off()
			
		}#end for (group in plot.groups)
		
	}#end for (target in targets)

	output.table = data.frame(fileID=fileID, shortID=shortID,targetID=targetID,
							length.lower5=length.lower5, length.lower.quantile=length.lower.quantile, length.med=length.med,
							length.mean=length.mean, length.upper.quantile=length.upper.quantile, length.upper5=length.upper5)
	write.table(output.table,"mothur_contig_length_quantiles.txt",sep="\t", row.names=F, quote=F)
}
} else{
	stop(paste("Need to write code for classifier: ",classifier,sep=""))
}