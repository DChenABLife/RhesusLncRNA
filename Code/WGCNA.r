
rm(list=ls())		#remove all the objects in the work space

Argus <- commandArgs(trailingOnly = TRUE)	#get the arguments.
if (length(Argus) < 2) {
	stop('Your input arguments is wrong!\n
		args1:\t the expression data file;\n 		
		args2:\t the out directory.\n\n'
	)
}
#args2:\t the annotation data file;\n
#args3:\t the Sample inf data file;\n

library('WGCNA')
library('gplots')

ExpFile <- Argus[1]
#AnnoFile <- Argus[2]
#Blocks <- as.numeric(Argus[3])
OutDir <- Argus[2]

CurrentDir <- getwd()

#The following settingis important, do not omit.
options(stringsAsFactors = FALSE)

femData = read.delim(file=ExpFile,header=T)		#read the Expression file, separated by tab
DataExpr = as.data.frame(t(femData[,-1]))		#remove the annotation column and transpose the frame
names(DataExpr) = femData[,1]					#add names to the Expression data
rownames(DataExpr) = names(femData)[-1]
nGenes = ncol(DataExpr)
probes = names(DataExpr)
nSamples = ncol(DataExpr)

GSG = goodSamplesGenes(DataExpr, verbose = 3)	# Iterative filtering of samples and genes with too many missing entries

DataExprG = DataExpr[GSG$goodSamples, GSG$goodGenes]		#remove the offending genes and samples from the data
DataExprT <- t(DataExprG)
# SampleTree = flashClust(dist(DataExprG), method = "average")		#there is no function named flashClust

# #sizeGrWindow(12,9)
# pdf(file = paste(OutDir,"/","sampleClustering.pdf",sep=""), width = 12, height = 9)
# par(cex = 0.6)
# par(mar = c(0,4,2,0))
# #plot the sample tree, and save as pdf format
# plot(sampleTree,
# 	main = "Sample clustering to detect outliers",sub="",
# 	xlab="", cex.lab = 1.5,
# 	cex.axis = 1.5, cex.main = 2
# )
# abline(h=15, col="red",lwd=2.5)
# dev.off()

#Loading clinical trait data

# Choose a set of soft-thresholding powers
powers = c(c(1:10),seq(from = 12, to=20, by=2))
# Call the networktopology analysis function
sft = pickSoftThreshold(DataExprG, powerVector = powers, verbose = 5)

pdf(file=paste(OutDir,"/","Soft_thresholding_detect.pdf",sep=''),width = 12, height = 9)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	xlab="Soft Threshold (power)",
	ylab="Scale Free Topology Model Fit,signed R^2",type="n",
	main = paste("Scale independence")
)
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
	xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
	main = paste("Mean connectivity")
)
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
cex=cex1,col="red")
dev.off()

#sft$powerEstimate				#there is something wrong in thr following code.
# if (is.na(sft$powerEstimate)) {
	# powerEstimate <- 30
# } else {
	# powerEstimate <- sft$powerEstimate
# }
net = blockwiseModules(DataExprG,
	power = sft$powerEstimate,		#the appropriate soft-thresholding power
	minModuleSize = 30,
	#blocks = ifelse(Blocks == 0, "", Blocks),
	maxBlockSize = 5000,
	reassignThreshold = 0, mergeCutHeight = 0.25,
	numericLabels = TRUE, pamRespectsDendro = FALSE,
	saveTOMs = TRUE,
	saveTOMFileBase = "BlockTOM",
	verbose = 3
)

moduleColors = labels2colors(net$colors)
geneTree = net$dendrograms[[1]]
pdf(file = paste(OutDir,"/","Dendrogram_and_module_colors.pdf",sep=""), width = 27, height = 6)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
	"Module\n colors",
	dendroLabels = FALSE, hang = 0.03,
	addGuide = TRUE, guideHang = 0.05,
	cex.main = 3, cex.axis = 2, cex.lab = 2,
	cex.colorLabels = 2
)
dev.off()

data <- data.frame(names(DataExprG),DataExprT,net$colors,moduleColors)
write.table(data,file = "WGCNA_module.txt", row.names = FALSE,quote=F,sep='\t')

HubGenes <- chooseTopHubInEachModule(DataExprG,moduleColors)			#select the hub gene
write.table (HubGenes,file = "HubGenes_of_each_module.xls",quote=F,sep='\t')

MEs0 = moduleEigengenes(DataExprG, moduleColors)$eigengenes			#Calculate module eigengenes
MEs = orderMEs(MEs0)			#Put close eigenvectors next to each other
rownames(MEs) <- rownames(DataExprG)
write.table(MEs,file = 'WGCNA_module_eigengenes.xls', quote = F,sep = '\t')		#row.names = FALSE,

for (i in 1:ncol(MEs)) {
	Name <- colnames(MEs[i])
	digit <- as.matrix(MEs[i])
	Color <- sub('ME','',Name)
	pdf(file=paste(Name,"_bar_plot.pdf",sep=""),width=9,height=6)
	barplot2(digit[,1],main='Module bar plot',
		ylim = c(min(digit[,1])-0.1,max(digit[,1])+0.1),			#min(digit[,1]),max(digit[,1])
		font.axis=2,cex.main=3,
		col = Color,
		las=3
	)
	dev.off()
}

#Plot the dendrogram and heatmap of module eigengene
pdf(file="Eigengene_dendrogram_heatmap.pdf",width = 12, height = 12)
par(cex = 2.0)
plotEigengeneNetworks(MEs,"",marDendro = c(0,8,2,4),marHeatmap = c(6,8,2,4), cex.lab = 1.6,xLabelsAngle = 90)
dev.off()

# Calculate topological overlap anew: this could
# be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(DataExprG, power = 6)
# Transform dissTOM with a power to make moderately strong
# connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function
#sizeGrWindow(9,9)
pdf(file='Network_heatmap_plot_all_genes.pdf',width=12,height=12)
TOMplot(plotTOM, geneTree, moduleColors,
main = "Network heatmap plot, all genes")
dev.off()

# # names (colors) of the modules
# modNames = substring(names(MEs), 3)
# # Module membership
# gmm = bicorAndPvalue(DataExprG, MEs)
# geneModuleMembership = as.data.frame(gmm$bicor)
# MMPvalue = as.data.frame(gmm$p)
# names(geneModuleMembership) = paste("MM", modNames, sep="")
# names(MMPvalue) = paste("p.MM", modNames, sep="")

# Annot <- read.delim(file=AnnoFile,header=T)
# probes2annot = match(probes,Annot$substanceBXH)

# geneInfo0 = data.frame(substanceBXH = probes[probes2annot],
# 	geneSymbol = Annot$gene_symbol[probes2annot],
# 	LocusLinkID = Annot$LocusLinkID[probes2annot],
# 	moduleColor = moduleColors,
# 	#geneTraitSignificance,
# 	#GSPvalue
# )
# modOrder = order(-abs(cor(MEs, weight, use = "p")))
# for (mod in 1:ncol(geneModuleMembership))
# {
# 	oldNames = names(geneInfo0)
# 	geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
# 	MMPvalue[, modOrder[mod]]);
# 	names(geneInfo0) = c(oldNames,
# 	paste("MM.", modNames[modOrder[mod]], sep=""),
# 	paste("p.MM.", modNames[modOrder[mod]], sep=""))
# }
# geneOrder = order(geneInfo0$moduleColor)
# geneInfo = geneInfo0[geneOrder, ]
# write.delim(geneInfo, file = "geneInfo.txt", row.names = FALSE)
#singleBlockMEs = moduleEigengenes(DataExpr, moduleColors)$eigengenes		# the module eigengenes