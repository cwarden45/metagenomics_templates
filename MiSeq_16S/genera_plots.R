plot.group = "Group"
color.group = "Group"

expression.file = "abundance.txt"
meta.file = "sample_description.txt"
genes = c("Lactobacillus","Escherichia/Shigella")

#you'll probably still want to modify mtext parameters...otherwise, most modifications can be done above

fixed.color.palatte = c("green","orange","purple","cyan","pink","maroon","yellow","grey","red","blue","black",colors())

meta.table = read.table(meta.file, sep="\t", head=T)
indexID = meta.table$userID
plot.groupID = meta.table[,plot.group]
color.groupID = meta.table[,color.group]

groupColors = rep("gray",times=length(color.groupID))
groups = levels(as.factor(as.character(color.groupID)))
for(i in 1:length(groups)){
	groupColors[color.groupID == groups[i]] = fixed.color.palatte[i]
}

expression.table = read.table(expression.file, sep="\t", head=T)
expression.mat = expression.table[,match(indexID, names(expression.table))]
all.genes = as.character(expression.table$Assignment)

for (gene in genes){
	print(gene)
	plot.gene = gsub("/",".",gene)
	gene.expr = as.numeric(expression.mat[all.genes == gene,])

	output.file = paste("genera_plot_",plot.gene,"_plot_by_",plot.group,"_color_by_",color.group,".png",sep="")
	png(output.file)
	par(mar=c(15,5,2,2))
	plot(as.numeric(plot.groupID), gene.expr, ylim=c(0,100), xaxt="n",
		xlab="", ylab="Genera Abundance", pch=16, col=groupColors, cex=1, type="p", main=gene)
	mtext(levels(as.factor(as.character(plot.groupID))), side=1, at =1:length(levels(as.factor(as.character(plot.groupID)))), las=2, line=2)
	legend("bottom", ncol=2, legend=groups, col = fixed.color.palatte[1:length(groups)],
			xpd=T, inset=-0.7, pch=16)
	dev.off()
}#end for (gene in genes)