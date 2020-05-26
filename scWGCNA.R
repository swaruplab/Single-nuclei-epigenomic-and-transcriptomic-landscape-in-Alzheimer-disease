####Code by Vivek Swarup, PhD UC Irvine #############
###contact vswarup@uci.edu################
#####################################################
##Single-nuclei Consensus Weighted Gene Co-expression Network Analysis (scWGCNA)
##ROSMAP data was downloaded and processed as mentioned in Morabito et al., https://doi.org/10.1101/695221
## Mathys et al data was downloaded from Synapse (syn18485175; doi:10.7303/syn18485175)
## Our snRNA-seq and Bulk-tissue RNA-seq data is available at https://www.synapse.org/#!Synapse:syn22079621/
#####################################################

library(Seurat)
library(WGCNA)
enableWGCNAThreads()
library(flashClust)

##Get the metacells for microglial (MG) from the snRNA-seq integrated data
metacell_seurat=readRDS('MG_joint_metacells_50.rds')
targets.MG=metacell_seurat@meta.data
group=factor(targets.MG$Diagnosis,c('Control','AD'))
datExpr.Cluster <- as.data.frame(GetAssayData(metacell_seurat, assay='RNA', slot='data')[VariableFeatures(metacell_seurat),])
datExpr.Cluster=as.data.frame(t(datExpr.Cluster))
rm(metacell_seurat)

ensembl=read.csv('hg38_GeneSymbol_ENSG.csv')
gnS=intersect(colnames(datExpr.Cluster),ensembl$Gene.name)
cat(length(gnS),'\n')

datExpr.Cluster=datExpr.Cluster[,match(gnS,colnames(datExpr.Cluster))]
ensembl1=ensembl[match(gnS,ensembl$Gene.name),]
colnames(datExpr.Cluster)=ensembl1$Gene.stable.ID
cat(length(gnS),'\n')

##load Bulk-tissue RNA-seq Data

load('/home/vivek/AD2019/cWGCNA/UCI_PFC_expression_metaData.rda')
load('/home/vivek/AD2019/cWGCNA/ROSMAP_DLPFC_expression_metaData.rda')
datExpr.ROSMAP=as.data.frame(t(normExpr.ROSMAP))

gnS=intersect(colnames(datExpr.Cluster),intersect(colnames(datExpr.ROSMAP),colnames(datExpr.UCI)))
cat(length(gnS),'\n')

datExpr.ROSMAP =datExpr.ROSMAP[,match(gnS,colnames(datExpr.ROSMAP))]
datExpr.UCI =datExpr.UCI[,match(gnS,colnames(datExpr.UCI))]
datExpr.Cluster =datExpr.Cluster [,match(gnS,colnames(datExpr.Cluster))]


# Make a multi-Expression data list containing all the data

nSets=3
setLabels=c("UCI","ROSMAP","Cluster.NucSeq.MG")
shortLabels=setLabels

multiExpr=vector(mode="list",length=nSets)

multiExpr[[1]] = list(data=as.data.frame(datExpr.UCI)) # UCI Data
names(multiExpr[[1]]$data)=colnames(datExpr.UCI)
rownames(multiExpr[[1]]$data)=rownames(datExpr.UCI)

multiExpr[[2]] = list(data=as.data.frame(datExpr.ROSMAP)) # ROSMAP Data
names(multiExpr[[2]]$data)=colnames(datExpr.ROSMAP)
rownames(multiExpr[[2]]$data)=rownames(datExpr.ROSMAP)

multiExpr[[3]] = list(data=as.data.frame(datExpr.Cluster)) #New Brain
names(multiExpr[[3]]$data)=colnames(datExpr.Cluster)
rownames(multiExpr[[3]]$data)=rownames(datExpr.Cluster)

checkSets(multiExpr) # check data size


multiMeta=list(UCI=list(data=targets.UCI),ROSMAP=list(data=targets.ROSMAP),Cluster.NucSeq.MG=list(data=targets.MG))

save(list=ls(),file="Consensus_UCI_ROSMAP_ClusterMG.rda")


## Network Construction

# Choose a set of soft-thresholding powers
powers = c(seq(1,10,by=1), seq(12,30, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
verbose = 100,networkType="signed",corFnc="bicor")[[2]]);


# Plot the results:
pdf("1_Power.pdf", height=10, width=18)

colors = c("blue", "red","black")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "mean connectivity",
"Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
for (col in 1:length(plotCols))
{
ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
}
}
# Plot the quantities in the chosen columns vs. the soft thresholding power

par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
if (set==1)
{
plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
main = colNames[col]);
addGrid();
}
if (col==1)
{
text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
labels=powers,cex=cex1,col=colors[set]);
} else
text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
labels=powers,cex=cex1,col=colors[set]);
if (col==1)
{
legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
} else
legend("topright", legend = setLabels, col = colors, pch = 20) ;
}
dev.off()


save(list=ls(),file="Consensus_UCI_ROSMAP_ClusterMG.rda")


## Consensus WGCNA

load("Consensus_UCI_ROSMAP_ClusterMG.rda")

softPower=18 ## power was chosen based on figure above

net=blockwiseConsensusModules(multiExpr, blocks = NULL,
                                         #maxBlockSize = 30000, ## This should be set to a smaller size if the user has limited RAM
                                         randomSeed = 12345,
                                         corType = "pearson", ## no use for bicor
                                         power = softPower,
                                         consensusQuantile = 0.3,
                                         networkType = "signed",
                                         TOMType = "unsigned",
                                         TOMDenom = "min",
                                         scaleTOMs = TRUE, scaleQuantile = 0.8,
                                         sampleForScaling = TRUE, sampleForScalingFactor = 1000,
                                         useDiskCache = TRUE, chunkSize = NULL,
                                         deepSplit = 4,
                                         pamStage=FALSE,
                                         detectCutHeight = 0.995, minModuleSize = 50,
                                         mergeCutHeight = 0.2,
                                         saveConsensusTOMs = TRUE,
                                         consensusTOMFilePattern = "ConsensusTOM-block.%b.rda")


save(list=ls(),file="Consensus_UCI_ROSMAP_ClusterMG.rda")

consMEs = net$multiMEs;
moduleLabels = net$colors;
# Convert the numeric labels to color labels
moduleColors = labels2colors(moduleLabels)
table(moduleColors)
consTree = net$dendrograms[[1]];
##UCI
MEs=moduleEigengenes(multiExpr[[1]]$data, colors = moduleColors, nPC=1)$eigengenes
MEs=orderMEs(MEs)

datExpr.UCI=multiExpr[[1]]$data

meInfo<-data.frame(rownames(datExpr.UCI), MEs)
colnames(meInfo)[1]= "SampleID"
KMEs<-signedKME(datExpr.UCI, MEs,outputColumnName = "kME",corFnc = "bicor")

ensembl=ensembl[na.omit(match(colnames(datExpr.UCI),ensembl$Gene.stable.ID)),]
geneInfo=as.data.frame(cbind(ensembl$Gene.stable.ID,ensembl$Gene.name,moduleColors, KMEs))

# merged gene symbol column
colnames(geneInfo)[1]= "Ensembl.Gene.ID"
colnames(geneInfo)[2]= "GeneSymbol"
colnames(geneInfo)[3]= "Initially.Assigned.Module.Color"

write.csv(geneInfo,file=paste('geneInfoSigned_UCI.csv',sep=''))


PCvalues.UCI=MEs
####ROSMAP
MEs=moduleEigengenes(multiExpr[[2]]$data, colors = moduleColors, nPC=1)$eigengenes
MEs=orderMEs(MEs)

datExpr.ROSMAP=multiExpr[[2]]$data

meInfo<-data.frame(rownames(datExpr.ROSMAP), MEs)
colnames(meInfo)[1]= "SampleID"
KMEs<-signedKME(datExpr.ROSMAP, MEs,outputColumnName = "kME",corFnc = "bicor")

ensembl=ensembl[na.omit(match(colnames(datExpr.ROSMAP),ensembl$Gene.stable.ID)),]
geneInfo=as.data.frame(cbind(ensembl$Gene.stable.ID,ensembl$Gene.name,moduleColors, KMEs))

# merged gene symbol column
colnames(geneInfo)[1]= "Ensembl.Gene.ID"
colnames(geneInfo)[2]= "GeneSymbol"
colnames(geneInfo)[3]= "Initially.Assigned.Module.Color"

write.csv(geneInfo,file=paste('geneInfoSigned_ROSMAP.csv',sep=''))


PCvalues.ROSMAP=MEs
PCvalues.ROSMAP=PCvalues.ROSMAP[,match(colnames(PCvalues.UCI),colnames(PCvalues.ROSMAP))]


MEs.snRNAseq=moduleEigengenes(multiExpr[[3]]$data, colors = moduleColors, nPC=1)$eigengenes
MEs.snRNAseq=orderMEs(MEs.snRNAseq)
PCvalues.snRNAseq=MEs.snRNAseq
PCvalues.snRNAseq=PCvalues.snRNAseq[,match(colnames(PCvalues.UCI),colnames(PCvalues.snRNAseq))]


pdf('ME_trajectory_Plot.MGClusters.pdf',width=26,height=12,useDingbats=F)
##
 toplot=t(PCvalues.UCI)
 toplot.ROSMAP=t(PCvalues.ROSMAP)
 toplot.snRNAseq=t(PCvalues.snRNAseq)

 cols=substring(colnames(MEs),3,20)
 par(mfrow=c(2,5))
 par(mar=c(12,6,4,2))

 for (i in 1:nrow(toplot)) {
   group.UCI=factor(targets.UCI$Path.Classification,c('Control','Early-Pathology AD','Late-Pathology AD'))
   l1=summary(lm(toplot[i,]~group.UCI))
   Early.pval=signif(l1$coefficients[2,4],2)
   Late.pval=signif(l1$coefficients[3,4],2)
   boxplot(toplot[i,]~group.UCI,col=c('blue','lightgreen','red'),ylab="ME Value",main=paste(rownames(toplot)[i],"\n","Early.Path pvalue=",Early.pval," Late.path pvalue=",Late.pval,sep=""),xlab=NULL,las=2,outpch = NA,xact="n")
   stripchart(toplot[i,]~group.UCI,vertical=T,method="jitter",pch=21,col='maroon',bg='bisque',add=T)
   axis(1,labels=F)
   axis(2,labels=F)

  boxplot(toplot[i,]~factor(targets.UCI$Sex,c('M','F')),col=cols[i],ylab="ME",main=rownames(toplot)[i],xlab=NULL,las=2)

  group.ROSMAP=factor(targets.ROSMAP$Path.Classification,c('Control','Early-Pathology AD','Late-Pathology AD'))
  l1=summary(lm(toplot.ROSMAP[i,]~group.ROSMAP))
  Early.pval=signif(l1$coefficients[2,4],2)
  Late.pval=signif(l1$coefficients[3,4],2)
  boxplot(toplot.ROSMAP[i,]~group.ROSMAP,col=c('blue','lightgreen','red'),ylab="ME Value",main=paste(rownames(toplot)[i],"\n","Early.Path pvalue=",Early.pval," Late.path pvalue=",Late.pval,sep=""),xlab=NULL,las=2,outpch = NA,xact="n")
  stripchart(toplot.ROSMAP[i,]~group.ROSMAP,vertical=T,method="jitter",pch=21,col='maroon',bg='bisque',add=T)
  axis(1,labels=F)
  axis(2,labels=F)

 boxplot(toplot.ROSMAP[i,]~as.factor(targets.ROSMAP$msex.x),col=cols[i],ylab="ME",main=rownames(toplot)[i],xlab=NULL,las=2)

  boxplot(toplot.snRNAseq[i,]~group,col=cols[i],ylab="ME",main=paste(rownames(toplot)[i],'snRNAseq'),xlab=NULL,las=2)


 }

dev.off()





pdf("SignedDendro_Consensus.pdf",height=10, width=15)
plotDendroAndColors(consTree, moduleColors, "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
main = "Consensus gene dendrogram and module colors")
dev.off()

##Final saving
save(list=ls(),file="Consensus_UCI_ROSMAP_ClusterMG.rda")
