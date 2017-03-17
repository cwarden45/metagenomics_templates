avgGroupExpression = function (geneExpr, groups) {
	avg.expr = tapply(geneExpr, groups, mean)
	return(avg.expr)
}#end def avgGroupExpression

ratio2fc = function(value)
{
	if(value >= 1){
		return(value)
	} else {
		return (-1/value)
	}
}#end def ratio2fc

calc.gene.cor = function(arr, indep.var)
{	
	na.count = length(arr[!is.na(arr)])
	if((na.count >= 3) & (sd(arr) != 0)){
		gene.cor.coef = cor(arr,indep.var)
	} else {
		gene.cor.coef = NA
	}
	return(gene.cor.coef)
}#end def calc.gene.cor

count.na.values = function(arr)
{
	return(length(arr[is.na(arr)]))
}#end def count.na.values

library(gplots)
fixed.color.palatte = c("green","orange","purple","cyan","pink","maroon","yellow","grey","red","blue","black","darkgreen","thistle1","tan","orchid1",colors())

param.table = read.table("parameters.txt", header=T, sep="\t")
comp.name=as.character(param.table$Value[param.table$Parameter == "comp_name"])
min.abundance= as.numeric(as.character(param.table$Value[param.table$Parameter == "abundance_cutoff"]))
fc.rounding.factor= as.numeric(as.character(param.table$Value[param.table$Parameter == "rounding_factor"]))
fc.cutoff= as.numeric(as.character(param.table$Value[param.table$Parameter == "fc_cutoff"]))
pvalue.cutoff = as.numeric(as.character(param.table$Value[param.table$Parameter == "pvalue_cutoff"]))
fdr.cutoff = as.numeric(as.character(param.table$Value[param.table$Parameter == "fdr_cutoff"]))
pvalue.cutoff2 = as.numeric(as.character(param.table$Value[param.table$Parameter == "sec_pvalue_cutoff"]))
fdr.cutoff2 = as.numeric(as.character(param.table$Value[param.table$Parameter == "sec_fdr_cutoff"]))
deg.groups = unlist(strsplit(as.character(param.table$Value[param.table$Parameter == "deg_groups"]), split=","))
plot.groups = unlist(strsplit(as.character(param.table$Value[param.table$Parameter == "plot_groups"]), split=","))
trt.group = as.character(param.table$Value[param.table$Parameter == "treatment_group"])
trt.group2 = as.character(param.table$Value[param.table$Parameter == "secondary_trt"])
interaction.flag = as.character(param.table$Value[param.table$Parameter == "interaction"])
interaction.flag[interaction.flag == "none"]="no"
pvalue.method = as.character(param.table$Value[param.table$Parameter == "pvalue_method"])
fdr.method = as.character(param.table$Value[param.table$Parameter == "fdr_method"])
output.folder = as.character(param.table$Value[param.table$Parameter == "Raw_Code_PC"])
user.folder = as.character(param.table$Value[param.table$Parameter == "Result_Folder"])
sample.description.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
counts.file = as.character(param.table$Value[param.table$Parameter == "counts_file"])
abundance.file= as.character(param.table$Value[param.table$Parameter == "abundance_file"])

setwd(output.folder)

sample.description.table = read.table(sample.description.file, sep="\t", header=T)
longID = sample.description.table$sampleID
sample.label = sample.description.table$userID

deg.group.table = sample.description.table[,deg.groups]
if (length(deg.groups) == 1){
	deg.meta = sample.description.table[!is.na(deg.group.table),]
} else {
	deg.grp.na.counts = apply(deg.group.table, 1, count.na.values)
	deg.meta = sample.description.table[deg.grp.na.counts == 0,]
}

ab.table = read.table(abundance.file, head=T, sep="\t")
assignments = as.character(ab.table$Assignment)
ab.mat= ab.table[,match(sample.label,names(ab.table))]
ab.mat[is.na(ab.mat)] = 0
ab.mat = matrix(as.numeric(unlist(ab.mat)), ncol=length(sample.label))
rownames(ab.mat) = as.character(assignments)
colnames(ab.mat) = as.character(sample.label)

counts.table = read.table(counts.file, head=T, sep="\t")
counts = counts.table[,match(sample.label,names(counts.table))]
counts[is.na(counts)] = 0
counts = matrix(as.numeric(unlist(counts)), ncol=length(sample.label))
rownames(counts) = as.character(assignments)
colnames(counts) = as.character(sample.label)

print(dim(ab.mat))
ab.max = apply(ab.mat, 1, max, na.rm = T)
assignments = assignments[ab.max > min.abundance]
ab.mat = ab.mat[match(assignments, rownames(ab.mat)),]
counts = counts[match(assignments, rownames(counts)),]
print(dim(ab.mat))

ab.mat = ab.mat + fc.rounding.factor

if(length(plot.groups) == 1){
	print("Averaging Abundance for One Variable (for plot.groups)")
	grp = sample.description.table[,plot.groups]
} else if ((length(plot.groups) == 2)&(interaction.flag == "no")){
	print("Averaging Abundance for First Variable (for plot.groups)")
	grp = sample.description.table[,plot.groups[1]]
} else if (length(plot.groups) == 2){
	print("Averaging Abundance for Interaction Variable (for plot.groups)")
	grp = paste(sample.description.table[,plot.groups[1]],sample.description.table[,plot.groups[2]],sep=":")
} else {
	stop("Code only compatible with 2 variables (with or without a 3rd interaction variable")
}

groupIDs = as.character(levels(as.factor(grp)))
average.ab = data.frame(t(apply(ab.mat, 1, avgGroupExpression, groups = grp)))
if(length(groupIDs) == 1){
	average.ab = t(average.ab)
} else {
	average.ab = average.ab
}
colnames(average.ab) = paste("avg.abundance", sub("-",".",groupIDs), sep=".")

#removed undefined group IDs (so, you can visualize samples that aren't in your comparison)
if(length(deg.groups) == 1){
	var1 = sample.description.table[,deg.groups]
	deg.counts = counts[,!is.na(var1)]
	deg.ab = ab.mat[,!is.na(var1)]
	
	ab.max = apply(deg.ab, 1, max, na.rm = T)
	assignments = assignments[ab.max > min.abundance]
	ab.mat = ab.mat[match(assignments, rownames(ab.mat)),]
	counts = counts[match(assignments, rownames(counts)),]
	deg.counts = deg.counts[match(assignments, rownames(deg.counts)),]
	deg.ab = deg.ab[match(assignments, rownames(deg.ab)),]
	average.ab = average.ab[match(assignments, rownames(deg.ab)),]
	print(dim(ab.mat))
	
	var1 = var1[!is.na(var1)]
} else if (length(deg.groups) == 2){
	if(interaction.flag == "filter-overlap"){
		var1 = sample.description.table[,deg.groups[1]]
		var2 = sample.description.table[,deg.groups[2]]
		deg.samples = !is.na(var1)&!is.na(var2)
		deg.counts = counts[,deg.samples]
		deg.ab = ab.mat[,deg.samples]

		ab.max = apply(deg.ab, 1, max, na.rm = T)
		assignments = assignments[ab.max > min.abundance]
		ab.mat = ab.mat[match(assignments, rownames(ab.mat)),]
		counts = counts[match(assignments, rownames(counts)),]
		deg.counts = deg.counts[match(assignments, rownames(deg.counts)),]
		deg.ab = deg.ab[match(assignments, rownames(deg.ab)),]
		average.ab = average.ab[match(assignments, rownames(average.ab)),]
		print(dim(ab.mat))

		var1 = var1[deg.samples]
		var2 = var2[deg.samples]
			
		prim.deg.grp = var1[var2 == trt.group2]
		prim.counts = counts[,var2 == trt.group2]
		prim.ab = ab.mat[,var2 == trt.group2]

		sec.deg.grp = var1[var2 != trt.group2]
		sec.counts = counts[,var2 != trt.group2]
		sec.ab = ab.mat[,var2 != trt.group2]

	} else{
		var1 = sample.description.table[,deg.groups[1]]
		var2 = sample.description.table[,deg.groups[2]]
		deg.samples = !is.na(var1)&!is.na(var2)
		deg.counts = counts[,deg.samples]
		deg.ab = ab.mat[,deg.samples]

		ab.max = apply(deg.ab, 1, max, na.rm = T)
		assignments = assignments[ab.max > min.abundance]
		ab.mat = ab.mat[match(assignments, rownames(ab.mat)),]
		counts = counts[match(assignments, rownames(counts)),]
		deg.counts = deg.counts[match(assignments, rownames(deg.counts)),]
		deg.ab = deg.ab[match(assignments, rownames(deg.ab)),]
		average.ab = average.ab[match(assignments, rownames(average.ab)),]
		print(dim(ab.mat))

		var1 = var1[deg.samples]
		var2 = var2[deg.samples]
	}
} else {
	stop("Code currently doesn't support more than 2 group model for DEG (with or without interaction)")
}

if(length(deg.groups) == 1){
	print("Averaging Abundance for One Variable (for deg.groups)")
	contrast.grp = var1
} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
	print("Averaging Abundance for First Variable (for deg.groups)")
	contrast.grp = var1
} else if (length(deg.groups) == 2){
	print("Averaging Abundance for Interaction Variable (for deg.groups)")
	contrast.grp = paste(var1,var2,sep=":")
} else {
	stop("Code only compatible with 2 variables (with or without a 3rd interaction variable")
}

if (trt.group == "continuous"){
	contrast.grp = as.numeric(contrast.grp)
	
	gene.cor = apply(deg.ab, 1, calc.gene.cor, indep.var=contrast.grp)
	
	fc.table = data.frame(cor=gene.cor)
} else {
	groupIDs = as.character(levels(as.factor(contrast.grp)))
	contrast.ab = data.frame(t(apply(deg.ab, 1, avgGroupExpression, groups = contrast.grp)))
	colnames(contrast.ab) = paste("avg.abundance", sub("-",".",groupIDs), sep=".")
}#end else

if((interaction.flag == "no") & (trt.group != "continuous")){
	print("Calculating fold-change for primary variable")
	trt.expr = contrast.ab[,paste("avg.abundance", sub("-",".",trt.group), sep=".")]
	cntl.expr = contrast.ab[,paste("avg.abundance", sub("-",".",groupIDs[groupIDs != trt.group]), sep=".")]

	fc = trt.expr / cntl.expr
	fc = round(sapply(fc, ratio2fc), digits=2)
	fc.table = data.frame(fold.change=fc)
} else if ((interaction.flag == "model")|(interaction.flag == "filter-overlap")){
	if ((trt.group == "continuous")&(trt.group2 == "continuous")){
		print("Calculating correlation for secondary variable")
		sec.contrast.grp = as.numeric(var2)
		
		gene.cor2 = apply(deg.ab, 1, calc.gene.cor, indep.var=sec.contrast.grp)
		
		fc.table = data.frame(prim.cor=gene.cor, sec.cor = gene.cor2)
	} else if (trt.group == "continuous"){
		print("Fold-change / correlation cutoff not used for mixed variable analysis")
		print("NOTE: 'Up-Regulated' R output refers to assignments that vary with FDR and p-value cutoffs")
		print("However, fold-change / correlation values for each separate variable are still provided")

		sec.groupIDs = var2
		sec.groups = as.character(levels(as.factor(sec.groupIDs)))
		sec.contrast.ab = data.frame(t(apply(deg.ab, 1, avgGroupExpression, groups = sec.groupIDs)))
		colnames(sec.contrast.ab) = paste("avg.abundance", sub("-",".",groupIDs), sep=".")
		sec.trt.expr = sec.contrast.ab[,paste("avg.abundance", sub("-",".",trt.group2), sep=".")]
		sec.cntl.expr = sec.contrast.ab[,paste("avg.abundance", sub("-",".",sec.groups[sec.groups != trt.group2]), sep=".")]

		sec.fc = sec.trt.expr  / sec.cntl.expr
		sec.fc = round(sapply(sec.fc, ratio2fc), digits=2)
		
		fc.table = data.frame(prim.cor=gene.cor, sec.fc = sec.fc)
	} else if (trt.group2 == "continuous"){	
		print("Fold-change / correlation cutoff not used for mixed variable analysis")
		print("NOTE: 'Up-Regulated' R output refers to assignments that vary with FDR and p-value cutoffs")
		print("However, fold-change / correlation values for each separate variable are still provided")

		prim.groupIDs = var1
		prim.groups = as.character(levels(as.factor(prim.groupIDs)))
		prim.contrast.ab = data.frame(t(apply(deg.ab, 1, avgGroupExpression, groups = prim.groupIDs)))
		colnames(prim.contrast.ab) = paste("avg.abundance", sub("-",".",prim.groups), sep=".")
		prim.trt = trt.group
		prim.cntl = prim.groups[prim.groups != trt.group]
		prim.trt.expr = prim.contrast.ab[,paste("avg.abundance", sub("-",".",prim.trt), sep=".")]
		prim.cntl.expr = prim.contrast.ab[,paste("avg.abundance", sub("-",".",prim.cntl), sep=".")]

		prim.fc = prim.trt.expr  /  prim.cntl.expr
		prim.fc = round(sapply(prim.fc, ratio2fc), digits=2)
		
		sec.contrast.grp = as.numeric(var2)
		gene.cor2 = apply(deg.ab, 1, calc.gene.cor, indep.var=sec.contrast.grp)
		
		fc.table = data.frame(prim.fc=prim.fc, sec.cor = gene.cor2)
	} else {
		print("Calculating fold-change table for primary variables (within subsets of secondary variable)")
		prim.groups = paste(var1,var2,sep=":")
		prim.trt = paste(trt.group,trt.group2,sep=":")
		prim.cntl = paste(prim.groups[prim.groups != trt.group],trt.group2,sep=":")
		prim.trt.expr = contrast.ab[,paste("avg.abundance", sub("-",".",prim.trt), sep=".")]
		prim.cntl.expr = contrast.ab[,paste("avg.abundance", sub("-",".",prim.cntl), sep=".")]

		prim.fc = prim.trt.expr  /  prim.cntl.expr
		prim.fc = round(sapply(prim.fc, ratio2fc), digits=2)

		sec.groups = as.character(levels(as.factor(sample.description.table[,deg.groups[2]])))
		sec.trt = paste(trt.group, sec.groups[sec.groups != trt.group2], sep=":")
		sec.cntl = paste(prim.groups[prim.groups != trt.group], sec.groups[sec.groups != trt.group2], sep=":")
		sec.trt.expr = contrast.ab[,paste("avg.abundance", sub("-",".",sec.trt), sep=".")]
		sec.cntl.expr = contrast.ab[,paste("avg.abundance", sub("-",".",sec.cntl), sep=".")]

		sec.fc = sec.trt.expr  / sec.cntl.expr
		sec.fc = round(sapply(sec.fc, ratio2fc), digits=2)

		#don't provide in table, but use for direction assignment
		overall.fc = prim.fc - sec.fc
		
		fc.table = data.frame(fc1 = prim.fc, fc2=sec.fc)
		colnames(fc.table) = c(paste("fold.change",trt.group,":",trt.group2,sep="."),
								paste("fold.change",trt.group,":",sec.groups[sec.groups != trt.group2], sep="."))
	}#end else
}else if(trt.group == "continuous"){
	print("Skipping fold-change calculation for continuous variable")
}else{
		stop("interaction must be \"no\", \"model\", or \"filter-overlap\"")
}#end else

rep.check = 1
for (i in 1:length(deg.groups)){
	deg.group = deg.groups[i]
	
	if((i == 1) & (trt.group != "continuous")){
		deg.group.values = as.character(deg.meta[,deg.group])
		min.reps = min(table(deg.group.values))
		if (min.reps < 2){
			rep.check=0
			print("There are not at least 2 samples per-group in order to calculate p-value.")
			print("In the future, please make sure you at least have duplicate samples.")
		}#end if (min.reps < 2)
	} else if ((i == 2) & (trt.group2 != "continuous")){
		deg.group.values = as.character(deg.meta[,deg.group])
		min.reps = min(table(deg.group.values))
		if (min.reps < 2){
			rep.check=0
			print("There are not at least 2 samples per-group in order to calculate p-value.")
			print("In the future, please make sure you at least have duplicate samples.")
		}#end if (min.reps < 2)
	} else if (i > 2){
		stop("Workflow currently doesn't support use of more than 2 variables")
	}
}#end for (deg.group in deg.groups)

if(rep.check == 1){
	#start p-value calculation

	if (pvalue.method == "metagenomeSeq"){
		library(metagenomeSeq)
		
		if ((length(deg.groups) == 2)&(interaction.flag == "filter-overlap")){
			print("metagenomeSeq, Two-Step Analysis")
			stop("Need to add code!")
		} else {
			if (length(deg.groups) == 1){
				print("metagenomeSeq with 1 variable")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}
				pheno = data.frame(var1=var1)
				rownames(pheno) = colnames(deg.counts)
				
				features = data.frame(genus = assignments)
				rownames(features) = rownames(deg.counts)
				
				mrObj = newMRexperiment(deg.counts, phenoData=AnnotatedDataFrame(pheno), featureData=AnnotatedDataFrame(features))
				scalePercentile = cumNormStatFast(mrObj)
				mrObj = cumNorm(mrObj, p = scalePercentile)
				design = model.matrix(~var1)
				fit = fitFeatureModel(mrObj,design)
				test.pvalue = as.numeric(fit$pvalues)
			} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
				print("metagenomeSeq with 2 variables")
				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				} else{
					var1 = as.factor(as.character(var1))
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				} else{
					var2 = as.factor(as.character(var2))
				}

				pheno = data.frame(var1=var1, var2=var2)
				rownames(pheno) = colnames(deg.counts)
				
				features = data.frame(genus = assignments)
				rownames(features) = rownames(deg.counts)
				
				mrObj = newMRexperiment(deg.counts, phenoData=AnnotatedDataFrame(pheno), featureData=AnnotatedDataFrame(features))
				scalePercentile = cumNormStatFast(mrObj)
				mrObj = cumNorm(mrObj, p = scalePercentile)
				design = model.matrix(~var1+var2)
				#get error message if trying to use fitFeatureModel
				settings = zigControl(maxit = 10, verbose = TRUE)
				fit = fitZig(obj = mrObj, mod = design, useCSSoffset = FALSE,
								control = settings)
				pvalue.mat = fit$eb$p.value
				test.pvalue = as.numeric(pvalue.mat[,2])
			} else if ((length(deg.groups) == 2)&(interaction.flag == "model")){
				print("metagenomeSeq with 2 variables plus interaction")
				stop("Need to add code!")
			}
		}#end else
	} else if (pvalue.method == "limma-ab"){
		library(limma)
		
		if ((length(deg.groups) == 2)&(interaction.flag == "filter-overlap")){
			print("limma, Two-Step Analysis")

			if (trt.group == "continuous"){
				prim.deg.grp = as.numeric(prim.deg.grp)
			}
			
			design = model.matrix(~prim.deg.grp)
			fit = lmFit(prim.ab, design)
			fit = eBayes(fit)
			pvalue.mat = data.frame(fit$p.value)
			prim.pvalue = pvalue.mat[,2]
			

			if (trt.group2 == "continuous"){
				sec.deg.grp = as.numeric(sec.deg.grp)
			}
			
			design = model.matrix(~sec.deg.grp)
			fit = lmFit(sec.ab,design)
			fit = eBayes(fit)
			pvalue.mat = data.frame(fit$p.value)
			sec.pvalue = pvalue.mat[,2]
			
		} else {
			if (length(deg.groups) == 1){
				print("limma with 1 variable")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}
				design = model.matrix(~var1)
				fit = lmFit(deg.ab,design)
				fit = eBayes(fit)
				pvalue.mat = data.frame(fit$p.value)
				test.pvalue = pvalue.mat[,2]
			} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
				print("limma with 2 variables")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				} else{
					var1 = as.factor(as.character(var1))
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				} else{
					var2 = as.factor(as.character(var2))
				}
				design = model.matrix(~var1 + var2)
				fit = lmFit(deg.ab,design)
				fit = eBayes(fit)
				pvalue.mat = data.frame(fit$p.value)
				test.pvalue = pvalue.mat[,2]
			} else if ((length(deg.groups) == 2)&(interaction.flag == "model")){
				print("limma with 2 variables plus interaction")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				}
				design = model.matrix(~var1*var2 + var1 + var2)
				fit = lmFit(deg.ab,design)
				fit = eBayes(fit)
				pvalue.mat = data.frame(fit$p.value)
				test.pvalue = pvalue.mat[,4]
			}
		}#end else
	} else if (pvalue.method == "limma-counts"){
		library(edgeR)
		library(limma)
		
		if ((length(deg.groups) == 2)&(interaction.flag == "filter-overlap")){
			print("limma-voom, Two-Step Analysis")

			if (trt.group == "continuous"){
				prim.deg.grp = as.numeric(prim.deg.grp)
			}
			
			y <- DGEList(counts=prim.counts, genes=assignments)
			png(paste(comp.name,"prim_voom_plot.png",sep="_"))
			design <- model.matrix(~prim.deg.grp)
			v <- voom(y,design,plot=TRUE)
			dev.off()
			fit <- lmFit(v,design)
			fit <- eBayes(fit)
			pvalue.mat = data.frame(fit$p.value)
			prim.pvalue = pvalue.mat[,2]
			

			if (trt.group2 == "continuous"){
				sec.deg.grp = as.numeric(sec.deg.grp)
			}
			
			y <- DGEList(counts=sec.counts, genes=assignments)
			design <- model.matrix(~sec.deg.grp)
			png(paste(comp.name,"sec_voom_plot.png",sep="_"))
			v <- voom(y,design,plot=TRUE)
			dev.off()
			fit <- lmFit(v,design)
			fit <- eBayes(fit)
			pvalue.mat = data.frame(fit$p.value)
			sec.pvalue = pvalue.mat[,2]
			
		} else {
			y <- DGEList(counts=deg.counts, genes=assignments)
			if (length(deg.groups) == 1){
				print("limma-voom with 1 variable")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}
				design <- model.matrix(~var1)
				png(paste(comp.name,"voom_plot.png",sep="_"))
				v <- voom(y,design,plot=TRUE)
				dev.off()
				fit <- lmFit(v,design)
				fit <- eBayes(fit)
				pvalue.mat = data.frame(fit$p.value)
				test.pvalue = pvalue.mat[,2]
			} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
				print("limma-voom with 2 variables")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				} else{
					var1 = as.factor(as.character(var1))
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				} else{
					var2 = as.factor(as.character(var2))
				}
				design <- model.matrix(~var1 + var2)
				png(paste(comp.name,"voom_plot.png",sep="_"))
				v <- voom(y,design,plot=TRUE)
				dev.off()
				fit <- lmFit(v,design)
				fit <- eBayes(fit)
				pvalue.mat = data.frame(fit$p.value)
				test.pvalue = pvalue.mat[,2]
			} else if ((length(deg.groups) == 2)&(interaction.flag == "model")){
				print("limma-voom with 2 variables plus interaction")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				}
				design <- model.matrix(~var1*var2 + var1 + var2)
				png(paste(comp.name,"voom_plot.png",sep="_"))
				v <- voom(y,design,plot=TRUE)
				dev.off()
				fit <- lmFit(v,design)
				fit <- eBayes(fit)
				pvalue.mat = data.frame(fit$p.value)
				test.pvalue = pvalue.mat[,4]
			}
		}#end else
	} else{
		stop("pvalue_method must be \"limma-ab\", \"limma-counts\", or \"metagenomeSeq\"")
	}
} else{
	test.pvalue = rep(1,times=length(assignments))
	prim.pvalue = rep(1,times=length(assignments))
	sec.pvalue = rep(1,times=length(assignments))
}#end else

if (trt.group == "continuous"){
	upID = "Increased Abundance"
	downID = "Decreased Abundance"
} else {
	upID = paste(trt.group," Up",sep="")
	downID = paste(trt.group," Down",sep="")	
}

if (interaction.flag == "no"){
	if (fdr.method == "BH"){
		fdr = p.adjust(test.pvalue, "fdr")
	} else if (fdr.method == "q-value"){
		library(qvalue)
		qobj <- qvalue(p = test.pvalue)
		fdr = qobj$qvalue
		png(paste(pvalue.method,"_qvalue_plot.png",sep=""))
		qHist = hist(qobj)
		print(qHist)
		dev.off()
	} else if (fdr.method == "q-lfdr"){
		library(qvalue)
		qobj <- qvalue(p = test.pvalue)
		fdr = qobj$lfdr
		png(paste(pvalue.method,"_qvalue_plot.png",sep=""))
		qHist = hist(qobj)
		print(qHist)
		dev.off()
	} else {
		stop("fdr_method must be \"BH\", \"q-value\", or \"q-lfdr\"")
	}
	status = rep("No Change", times=length(fdr))
	if (trt.group == "continuous"){
		status[(gene.cor >= 0) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = upID
		status[(gene.cor <= 0) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = downID
		pvalue.table = data.frame(p.value = test.pvalue, FDR = fdr)
	} else{
		status[(fc >= fc.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = upID
		status[(fc <= -fc.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = downID
		pvalue.table = data.frame(p.value = test.pvalue, FDR = fdr)
	}#end else
} else{
	trt.group = prim.trt
	if(interaction.flag == "model"){
		if (fdr.method == "BH"){
			fdr = p.adjust(test.pvalue, "fdr")
		} else if (fdr.method == "q-value"){
			library(qvalue)
			qobj <- qvalue(p = test.pvalue)
			fdr = qobj$qvalue
			png(paste(pvalue.method,"_qvalue_plot.png",sep=""))
			qHist = hist(qobj)
			print(qHist)
			dev.off()
		} else if (fdr.method == "q-lfdr"){
			library(qvalue)
			qobj <- qvalue(p = test.pvalue)
			fdr = qobj$lfdr
			png(paste(pvalue.method,"_qvalue_plot.png",sep=""))
			qHist = hist(qobj)
			print(qHist)
			dev.off()
		} else {
			stop("fdr_method must be \"BH\", \"q-value\", or \"q-lfdr\"")
		}
		status = rep("No Change", times=length(fdr))
		if ((trt.group == "continuous")&(trt.group2 == "continuous")){
			status[(gene.cor.int >= 0) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = upID
			status[(gene.cor.int <= 0) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = downID
			pvalue.table = data.frame(p.value = test.pvalue, FDR = fdr)
		} else if ((trt.group != "continuous")&(trt.group2 != "continuous")){
			status[(overall.fc >= fc.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = upID
			status[(overall.fc <= -fc.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = downID
			pvalue.table = data.frame(p.value = test.pvalue, FDR = fdr)
		} else {
			upID = "Variable Abundance"
			status[(test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = upID
			pvalue.table = data.frame(p.value = test.pvalue, FDR = fdr)
		}#end else
	} else if (interaction.flag == "filter-overlap"){
		if (fdr.method == "BH"){
			fdr = p.adjust(prim.pvalue, "fdr")
		} else if (fdr.method == "q-value"){
			library(qvalue)
			qobj <- qvalue(p = prim.pvalue)
			fdr = qobj$qvalue
			png(paste(pvalue.method,"_qvalue_plot.png",sep=""))
			qHist = hist(qobj)
			print(qHist)
			dev.off()
		} else if (fdr.method == "q-lfdr"){
			library(qvalue)
			qobj <- qvalue(p = prim.pvalue)
			fdr = qobj$lfdr
			png(paste(pvalue.method,"_qvalue_plot.png",sep=""))
			qHist = hist(qobj)
			print(qHist)
			dev.off()
		} else {
			stop("fdr_method must be \"BH\", \"q-value\", or \"q-lfdr\"")
		}
		pass1.status = rep("No Change", times=length(fdr))
		if (trt.group == "continuous"){
			pass1.status[(gene.cor >= 0) & (prim.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = paste(trt.group," Up",sep="")
			pass1.status[(gene.cor <= 0) & (prim.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = paste(trt.group," Down",sep="")
		} else{
			pass1.status[(prim.fc >= fc.cutoff) & (prim.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = paste(trt.group," Up",sep="")
			pass1.status[(prim.fc <= -fc.cutoff) & (prim.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = paste(trt.group," Down",sep="")
		}#end else

		print(paste("Primary Up-Regulated: ",length(pass1.status[pass1.status == paste(trt.group," Up",sep="")]),sep=""))
		print(paste("Primary Down-Regulated: ",length(pass1.status[pass1.status == paste(trt.group," Down",sep="")]),sep=""))

		if (fdr.method == "BH"){
			sec.fdr = p.adjust(sec.pvalue, "fdr")
		} else if (fdr.method == "q-value"){
			library(qvalue)
			qobj <- qvalue(p = sec.pvalue)
			sec.fdr = qobj$qvalue
			png(paste(pvalue.method,"_qvalue_plot.png",sep=""))
			qHist = hist(qobj)
			print(qHist)
			dev.off()
		} else if (fdr.method == "q-lfdr"){
			library(qvalue)
			qobj <- qvalue(p = sec.pvalue)
			sec.fdr = qobj$lfdr
			png(paste(pvalue.method,"_qvalue_plot.png",sep=""))
			qHist = hist(qobj)
			print(qHist)
			dev.off()
		} else {
			stop("fdr_method must be \"BH\", \"q-value\", or \"q-lfdr\"")
		}		

		pass2.status = rep("No Change", times=length(fdr))
		if (trt.group2 == "continuous"){
			pass2.status[(gene.cor2 >= 0) & (sec.pvalue <= pvalue.cutoff2) & (sec.fdr <= fdr.cutoff2)] = paste(trt.group," Up",sep="")
			pass2.status[(gene.cor2 <= 0) & (sec.pvalue <= pvalue.cutoff2) & (sec.fdr <= fdr.cutoff2)] = paste(trt.group," Down",sep="")
		} else{
			pass2.status[(sec.fc >= fc.cutoff) & (sec.pvalue <= pvalue.cutoff2) & (sec.fdr <= fdr.cutoff2)] = paste(trt.group," Up",sep="")
			pass2.status[(sec.fc <= -fc.cutoff) & (sec.pvalue <= pvalue.cutoff2) & (sec.fdr <= fdr.cutoff2)] = paste(trt.group," Down",sep="")
		}#end else

		print(paste("Secondary Up-Regulated: ",length(pass2.status[pass2.status == paste(trt.group," Up",sep="")]),sep=""))
		print(paste("Secondary Down-Regulated: ",length(pass2.status[pass2.status == paste(trt.group," Down",sep="")]),sep=""))
			
		pvalue.table = data.frame(prim.pvalue = prim.pvalue, prim.FDR = fdr,
										sec.pvalue=sec.pvalue, sec.fdr=sec.fdr)
			
		status = rep("No Change", times=length(fdr))
		status[(pass1.status == paste(trt.group," Up",sep="")) & (pass2.status == "No Change")] = upID
		status[(pass1.status == paste(trt.group," Down",sep="")) & (pass2.status == "No Change")] = downID
	} else{
		stop("interaction must be \"no\", \"model\", or \"filter-overlap\"")
	}#end else
}#end else

print(paste("Up-Regulated: ",length(status[status == upID]),sep=""))
print(paste("Down-Regulated: ",length(status[status == downID]),sep=""))

if (interaction.flag == "filter-overlap"){
	pvalue.method = paste(pvalue.method,"two-step_filtered",sep="_")
}

if(rep.check == 1){
	deg.table = data.frame(genus = assignments,	average.ab, fc.table, pvalue.table, status = status)
} else {
	deg.table = data.frame(genus = assignments,	average.ab, fc.table, status = status)	
}#end else

deg.file = paste(comp.name,"_",pvalue.method,"_differential_abundance_fdr_",fdr.cutoff,"_pval_",pvalue.cutoff,".txt",sep="")
deg.file = gsub(":",".",deg.file)
write.table(deg.table, file=deg.file, row.names=F, quote=F, sep="\t")

final.deg.file = paste(user.folder,"/",comp.name,"_differential_abundance_stats.txt",sep="")
write.table(deg.table, file=final.deg.file, row.names=F, quote=F, sep="\t")

temp.ab = ab.mat
temp.ab = temp.ab[status != "No Change", ]
deg.assignments = assignments[status != "No Change"]

if(length(deg.assignments) > 1){
	if(length(plot.groups) > 1){
		source("heatmap.3.R")
		grp1 = as.character(sample.description.table[,plot.groups[1]])
		grp2 = as.character(sample.description.table[,plot.groups[2]])
		group.levels = c(levels(as.factor(grp1)),levels(as.factor(grp2)))

		color.palette <- fixed.color.palatte[1:length(group.levels)]
		labelColors1 = rep("black",times=length(sample.label))
		for (i in 1:length(group.levels)){
			labelColors1[grp1 == as.character(group.levels[i])] = color.palette[i]
		}#end for (i in 1:length(group.levels))
		labelColors2 = rep("black",times=length(sample.label))
		for (i in 1:length(group.levels)){
			labelColors2[grp2 == as.character(group.levels[i])] = color.palette[i]
		}#end for (i in 1:length(group.levels))
		
		temp.ab = t(temp.ab)
		if(length(deg.assignments) < 25){
			colnames(temp.ab) = deg.assignments
		} else {
			colnames(temp.ab) = rep("", length(deg.assignments))
		}
		rownames(temp.ab) = sample.label

		column_annotation <- as.matrix(deg.assignments)
		colnames(column_annotation) <- c("")

		row_annotation <- data.frame(label1 = labelColors1, label2 = labelColors2)
		row_annotation = as.matrix(t(row_annotation))
		rownames(row_annotation) <- c(plot.groups)

		heatmap.file <- paste(comp.name,"_",pvalue.method,"_differential_abundance_fdr_",fdr.cutoff,"_pval_",pvalue.cutoff,".png",sep="")
		heatmap.file = gsub(":",".",heatmap.file)
		png(file = heatmap.file)
		heatmap.3(temp.ab, col=colorpanel(33, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
					RowSideColors=row_annotation, trace="none", margins = c(20,20),RowSideColorsSize=4, dendrogram="both")
		legend("topright", legend=group.levels,
							col=color.palette,
							pch=15, cex=0.7)
		dev.off()
		
		if(interaction.flag != "no"){
			temp.fc.table = as.matrix(fc.table)
			if (((trt.group == "continuous") & (trt.group2 == "continuous")) | ((trt.group != "continuous") & (trt.group2 != "continuous"))){
				temp.fc.table = temp.fc.table[,-ncol(temp.fc.table)]
			}
			temp.fc.table = temp.fc.table[status != "No Change", ]
			if(length(deg.assignments) < 25){
				rownames(temp.fc.table) = deg.assignments
			} else {
				rownames(temp.fc.table) = rep("",times=length(deg.assignments))
			}
			colnames(temp.fc.table) = gsub(".:.",":",gsub("fold.change.","",colnames(temp.fc.table)))
		
			temp.fc.table[temp.fc.table < -10] = -10
			temp.fc.table[temp.fc.table > 10] = 10
		
			heatmap.file <- paste("fold_change_",comp.name,"_",pvalue.method,"_differential_abundance_fdr_",fdr.cutoff,"_pval_",pvalue.cutoff,".png",sep="")
			heatmap.file = gsub(":",".",heatmap.file)
			png(file = heatmap.file)
			heatmap.2(temp.fc.table, col=colorpanel(33, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
						trace="none", margins = c(20,20), cexCol=1.5)
			dev.off()
		}#end if(interaction.flag != "no")
		
	} else {
		group.levels = levels(as.factor(sample.description.table[,plot.groups]))

		color.palette <- fixed.color.palatte[1:length(group.levels)]
		labelColors = rep("black",times=length(sample.label))
		for (i in 1:length(group.levels)){
			labelColors[grp == as.character(group.levels[i])] = color.palette[i]
		}#end for (i in 1:length(group.levels))

		temp.ab = t(temp.ab)
		if(length(deg.assignments) < 25){
			colnames(temp.ab) = deg.assignments
		} else {
			colnames(temp.ab) = rep("", length(deg.assignments))
		}
		rownames(temp.ab) = sample.label
		
		heatmap.file = paste(comp.name,"_",pvalue.method,"_differential_abundance_fdr_",fdr.cutoff,"_pval_",pvalue.cutoff,".png",sep="")
		heatmap.file = gsub(":",".",heatmap.file)
		png(file = heatmap.file)
		heatmap.2(temp.ab, col=colorpanel(33, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
					 RowSideColors=labelColors, trace="none", margins = c(20,20))
		legend("topright", group.levels, col=color.palette, pch=15, cex=0.8)
		dev.off()
	}#end else
}#end if(length(deg.assignments) > 1)
