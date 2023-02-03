setwd("E:/Osteoclast_singlecell/GSE125088")

library("Seurat")
library("monocle")
library("ggplot2")
library("reshape2")
library("stringr")
library("ggpubr")

gse125088 <- Read10X(data.dir = "E:/Osteoclast_singlecell/GSE125088/data")
gse125088 <- CreateSeuratObject(gse125088)

gse125088 <- PercentageFeatureSet(gse125088, pattern = 'mt-', col.name = 'Mitochondrial')

ribos <- readLines("C:/Users/asp1d/Desktop/CommonLists/mouse_ribogenes.txt")
ribos <- intersect(ribos, row.names(gse125088@assays$RNA@data))
gse125088 <- PercentageFeatureSet(gse125088, features = ribos, col.name = 'Ribosomal')

gse125088 <- subset(gse125088, subset = Mitochondrial < 10)
gse125088 <- subset(gse125088, subset = nFeature_RNA < 6000)
gse125088 <- subset(gse125088, subset = nFeature_RNA > 1500)

gse125088 <- NormalizeData(gse125088)
gse125088 <- FindVariableFeatures(gse125088, nfeatures = 4000)
gse125088 <- ScaleData(gse125088)
gse125088 <- RunPCA(gse125088)
ElbowPlot(gse125088)

gse125088 <- RunUMAP(gse125088, dims = 1:30, reduction = 'pca')
FeaturePlot(gse125088, features = c("Acp5", "Ctsk", "Dcstamp", "Nfatc1"))

gse125088 <- FindNeighbors(gse125088, reduction = 'umap', dims = 1:2)
gse125088 <- FindClusters(gse125088)

posost <- subset(gse125088, idents = c('2', '8', '9', '22', '29', '32'))
posost <- FindVariableFeatures(posost, nfeatures = 4000)
posost <- ScaleData(posost)
posost <- RunPCA(posost)
posost <- RunUMAP(posost, dims = 1:30, reduction = 'pca')

posost <- FindNeighbors(posost, reduction = 'umap', dims = 1:2)
posost <- FindClusters(posost)

osteo125088 <- subset(posost, idents = c('14', '19'))

osteo125088 <- FindVariableFeatures(osteo125088, nfeatures = 4000)
osteo125088 <- ScaleData(osteo125088)

###############
####monocle####
###############
V3toMonocle <- function(x){
fdata           <-  data.frame(gene_short_name = rownames(x@assays$RNA@data))
rownames(fdata) <- rownames(x@assays$RNA@data)
fd <- new('AnnotatedDataFrame', data = fdata) 
pd <- new('AnnotatedDataFrame', data = x@meta.data) 
newCellDataSet(x@assays$RNA@data,phenoData = pd, featureData = fd, 
                           lowerDetectionLimit = 0.01)
}

readGMT <- function (file) 
{
  if (!grepl("\\.gmt$", file)[1]) {
    stop("Pathway information must be a .gmt file")
  }
  geneSetDB = readLines(file)
  geneSetDB = strsplit(geneSetDB, "\t")
  names(geneSetDB) = sapply(geneSetDB, "[", 1)
  geneSetDB = lapply(geneSetDB, "[", -1:-2)
  geneSetDB = lapply(geneSetDB, function(x) {
    x[which(x != "")]
  })
  return(geneSetDB)
}


########
##run monocle
#######

ostmono <- V3toMonocle(osteo125088)
ostmono <- setOrderingFilter(ostmono, ordering_genes = osteo125088@assays$RNA@var.features)
ostmono <- estimateSizeFactors(ostmono)
ostmono <- estimateDispersions(ostmono)
ostmono <- reduceDimension(ostmono, reduction_method = 'DDRTree', norm_method = 'none')
ostmono  <- orderCells(ostmono)

plot_cell_trajectory(ostmono, color_by = 'seurat_clusters', markers = 'Acp5')


reactome <- readGMT("C:/Users/asp1d/Desktop/CommonLists/Pathways_Reactome.gmt")
ostavg <- AverageExpression(osteo125088)
ostavg <- ostavg$RNA
ostavg <- ostavg[rowSums(ostavg) > 0, ]
exprgenes <- row.names(ostavg)

reactome <- lapply(reactome, function(x) str_to_title(x))
reactome <- lapply(reactome, function(x) intersect(x, exprgenes))
lenreact <-  lapply(reactome, function(x) length(x))
reactome <- reactome[lenreact > 19]
lenreact <- lenreact[lenreact > 19]
reactome <- reactome[lenreact < 500]
lenreact <- lenreact[lenreact < 500]

compmono <- lapply(reactome, function(reactome) 
  setOrderingFilter(ostmono, ordering_genes = reactome)) #list of monocles

compmono <- lapply(compmono, function(compmono) 
  reduceDimension(compmono, reduction_method = "DDRTree", norm_method = "none", auto_param_selection = F))

compmono <- lapply(compmono, function(compmono) orderCells(compmono))

dftime <- lapply(compmono, function(compmono) compmono$Pseudotime)
dftime <- as.data.frame(dftime)
row.names(dftime) <- WhichCells(osteo125088)

cortime   <- lapply(1:length(reactome), function(x) cor(ostmono$Pseudotime, dftime[[x]])) #pseudotime correlations
dualcomp  <- as.data.frame(cortime)
dualcomp <- t(dualcomp)
dualcomp <- as.data.frame(dualcomp)
colnames(dualcomp)[1] = "Correlation_Order"
row.names(dualcomp) <- names(reactome)
lenreact <- as.data.frame(t(as.data.frame(lenreact)))
dualcomp$SetSize <- lenreact$V1


p1 <- plot_cell_trajectory(ostmono, color_by = 'Pseudotime', show_branch_points = F)
p1$layers[[1]]$aes_params$alpha <- 0.7
p1$layers[[2]]$aes_params$alpha <- 0.7
p1$layers[[2]]$aes_params$stroke <- 1

##correlation with pseudotime
pick6 <- c("Acp5", "Atp6v0d2", "Nfatc1", "Ctsk", "Dcstamp", "Mmp9")
ostsca <- osteo125088@assays$RNA@data
ostsca <- as.data.frame(t(ostsca))
ostsca <- ostsca[pick6]
ostsca$Pseudotime <- ostmono$Pseudotime
ostsca <- melt(ostsca, id.vars = 'Pseudotime')
colnames(ostsca)[2] <- 'Gene'
colnames(ostsca)[3] <- 'Expression'

ggplot(ostsca, aes(Pseudotime, Expression, fill = Gene))+
  geom_point(shape = 21, alpha = 0.6, stroke = 0.6)+
  geom_smooth(method = 'loess', aes(color = Gene))+
  theme_classic()+
  theme(text = element_text(size = 14, face = 'bold'))


