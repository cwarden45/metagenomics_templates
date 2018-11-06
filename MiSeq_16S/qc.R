colLab = function(n, labelColors, clusMember) { 
   if(is.leaf(n)) { 
       a <- attributes(n) 
	   #print(a)
       # clusMember - vector of sample names (ordered to match label color.palette)
       # labelColors - a vector of color.palette for the above grouping 
       labCol <- labelColors[clusMember == a$label]
	   #print(labCol)
       attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol) 
   } 
   n 
}#end def colLab

reformat.label = function(name, kept.characters=20){
	#new.name = paste("...",unlist(substr(name, nchar(name)-kept.characters,nchar(name))),sep="")
	assignment.info = unlist(strsplit(name,split=";"))
	new.name = assignment.info[length(assignment.info)]
	return(new.name)
}#end def reformat.label

param.table = read.table("parameters.txt", header=T, sep="\t")
sample.description.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
classification.file = as.character(param.table$Value[param.table$Parameter == "abundance_file"])
classifier = as.character(param.table$Value[param.table$Parameter == "Classifier"])
cluster.distance = as.character(param.table$Value[param.table$Parameter == "cluster_distance"])
plot.groups = unlist(strsplit(as.character(param.table$Value[param.table$Parameter == "plot_groups"]), split=","))

sample.table = read.table(sample.description.file, sep="\t", header=T)
userID = as.character(sample.table$userID)

ab.table = read.table(classification.file, sep="\t", header=T)
assignment = ab.table$Assignment
ab.mat = ab.table[,match(userID, names(ab.table))]
ab.mean = apply(ab.mat, 1, mean, na.rm = T)
ab.mat= ab.mat[order(ab.mean, decreasing = TRUE),]
ab.mat[is.na(ab.mat)] = 0
assignment = assignment[order(ab.mean, decreasing = TRUE)]

#barplot
png(paste(classifier,"_top_barplot.png",sep=""))
class.colors = rep("gray",times=length(assignment))
top.colors = c("green","blue","red","brown","orange","purple","cyan","yellow","pink","maroon","cornflowerblue","burlywood","darkseagreen1","gold","deeppink","yellowgreen","violet","tomato","springgreen","olivedrab")
ab.max = apply(ab.mat, 1, max, na.rm = T)
colored.assignments = assignment[ab.max > 1]
if(length(colored.assignments) > 20){
colored.assignments = colored.assignments[1:20]
}
for (i in 1:length(colored.assignments)){
	class.colors[assignment == colored.assignments[i]] = top.colors[i]
}
	
par(mfrow=c(1,2),mar=c(7,2,2,2),oma=c(7,2,2,2))
barplot(as.matrix(ab.mat), names.arg=names(ab.mat), las=2, col=class.colors, ylim = c(0,100))

plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n',xlab="",ylab="")
if(classifier == "BWA"){
	legend.labels = sapply(as.character(colored.assignments), reformat.label)
	legend("left",as.character(legend.labels),
			col=top.colors[1:length(colored.assignments)], pch=15, cex=0.9, xpd=T, inset=-0.1)

}else{
	legend("left",as.character(colored.assignments),
		col=top.colors[1:length(colored.assignments)], pch=15, cex=0.9, xpd=T, inset=-0.1)
}
dev.off()

#plot per Group
groupIDs = levels(sample.table$Group)

for (groupID in groupIDs){
	print(groupID)
	donor.table = sample.table[sample.table$Group == groupID,]
	
	donor.userID = as.character(donor.table$userID)

	if(nrow(donor.table) > 1){
		donor.userID = as.character(donor.table$userID)
		
		donor.assignment = ab.table$Assignment
		donor.ab.mat = ab.table[,match(donor.userID, names(ab.table))]
		donor.ab.mean = apply(donor.ab.mat, 1, mean, na.rm = T)
		donor.ab.mat= donor.ab.mat[order(donor.ab.mean, decreasing = TRUE),]
		donor.ab.mat[is.na(donor.ab.mat)] = 0
		donor.assignment = donor.assignment[order(donor.ab.mean, decreasing = TRUE)]
		
		donor.class.colors = rep("gray",times=length(donor.assignment))
		for (i in 1:length(colored.assignments)){
			donor.class.colors[donor.assignment == colored.assignments[i]] = top.colors[i]
		}

		png(paste(classifier,"_top_barplot_",groupID,".png",sep=""))
		par(mfrow=c(1,2),mar=c(10,2,2,2),oma=c(4,2,2,2))
		barplot(as.matrix(donor.ab.mat), names.arg=donor.userID, las=2,
					col=donor.class.colors, ylim = c(0,100), main=groupID)

		plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n',xlab="",ylab="")
		legend("left",as.character(colored.assignments),
				col=top.colors[1:length(colored.assignments)], pch=15, cex=0.9, xpd=T, inset=-0.1)

		dev.off()
	}else{
		donor.userID = as.character(donor.table$userID)
		donor.ab = ab.table[,match(donor.userID, names(ab.table))]
		names(donor.ab)=ab.table$Assignment
		donor.ab = donor.ab[!is.na(donor.ab)&(donor.ab > 1)]
		
		pie.colors = rep("gray", length(donor.ab))
		
		for (i in 1:length(donor.ab)){
			temp.genus = names(donor.ab)[i]
			if (temp.genus %in% colored.assignments)
				pie.colors[i]=top.colors[colored.assignments == temp.genus]
		}
		
		png(paste(classifier,"_top_barplot_",groupID,".png",sep=""))
		par(mar=c(2,7,2,10))
		pie(donor.ab, col=pie.colors, main=donor.userID)
		dev.off()
	}
}#end for for (groupID in groupIDs)

for (group in plot.groups){
	fixed.color.palatte = c("green","orange","purple","cyan","pink","maroon","yellow","grey","red","blue","black","darkgreen","thistle1","tan","orchid1",colors())

	print(group)
	temp.mat = ab.mat[,!is.na(sample.table[,group])]
	print(dim(temp.mat))
	qc.grp = sample.table[,group]
	qc.grp = qc.grp[!is.na(sample.table[,group])]
	clusterID = userID[!is.na(sample.table[,group])]
	
	pca.values <- prcomp(na.omit(data.matrix(temp.mat)))
	pc.values <- data.frame(pca.values$rotation)
	variance.explained <- (pca.values $sdev)^2 / sum(pca.values $sdev^2)
	pca.table <- data.frame(PC = 1:length(variance.explained), percent.variation = variance.explained, t(pc.values))

	pca.text.file = paste(group,"_pca_values.txt",sep="")
	write.table(pca.table, pca.text.file, quote=F, row.names=F, sep="\t")
	
	groups = levels(as.factor(as.character(qc.grp)))
	num.sample.types = length(groups)
	color.palette <- fixed.color.palatte[1:length(groups)]

	labelColors = rep("black",times=ncol(temp.mat))
	for (i in 1:length(groups)){
		labelColors[qc.grp == as.character(groups[i])] = color.palette[i]
	}#end for (i in 1:length(groups))

	pca.file = paste("pca_by_",group,".png",sep="")
	png(file=pca.file)
	plot(pc.values$PC1, pc.values$PC2, col = labelColors, xlab = paste("PC1 (",round(100* variance.explained[1] , digits = 2),"%)", sep = ""),
			ylab = paste("PC2 (",round(100* variance.explained[2] , digits = 2),"%)", sep = ""), pch=19)
	legend("topleft",legend=groups,col=color.palette,  pch=19)
	dev.off()


	cluster.file = paste("cluster_by_",group,"_",cluster.distance,".png",sep="")
	if(cluster.distance == "Euclidean"){
		dist1 <- dist(as.matrix(t(temp.mat)))
	}else if (cluster.distance == "Pearson_Dissimilarity"){
		cor.mat = cor(as.matrix(temp.mat))
		dis.mat = 1 - cor.mat
		dist1 <- as.dist(dis.mat)
	}else if (cluster.distance == "Spearman_Dissimilarity"){
		cor.mat = cor(as.matrix(temp.mat), method="spearman")
		dis.mat = 1 - cor.mat
		dist1 <- as.dist(dis.mat)
	}else if (cluster.distance == "Kendall_Dissimilarity"){
		cor.mat = cor(as.matrix(temp.mat), method="kendall")
		dis.mat = 1 - cor.mat
		dist1 <- as.dist(dis.mat)
	}else{
		stop("cluster_distance must be 'Euclidean' or 'Pearson_Dissimilarity'")
	}
	clusMember <- groups
	hc <- hclust(dist1)
	dend1 <- as.dendrogram(hc)
	png(file = cluster.file)
	dend1 <- dendrapply(dend1, colLab, labelColors=labelColors, clusMember=clusterID) 
	a <- attributes(dend1) 
	attr(dend1, "nodePar") <- c(a$nodePar, lab.col = labelColors) 
	op <- par(mar = par("mar") + c(0,0,0,10)) 
	plot(dend1, horiz=T)
	par(op) 
	dev.off()
	
}#end for (group in plot.groups)
