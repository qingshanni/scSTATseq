##########
########
######this is the re-analysis of the dataset presented in Tsukaki et al 
###Nature Metabolism 2020 of in vitro stimulated (RANKL+MCSF) BMDMs
#####includes sampling at three timepoints

setwd("E:/Osteoclast_singlecell/GSE147174/")

library("plyr")
library("Seurat")
library("monocle")
library("ggplot2")
library("reshape2")
library("ggpubr")
library("harmony")
library("stringr")


#########
###first read in the raw files downloaded from GSE147174
#########
day0  <- Read10X("E:/Osteoclast_singlecell/GSE147174/data/Day0")
day1  <- Read10X("E:/Osteoclast_singlecell/GSE147174/data/Day1")
day3  <- Read10X("E:/Osteoclast_singlecell/GSE147174/data/Day3")


###annotations are given as ENSMUSG instead of gene names
##so we need to do an adjustment first


biomart <- read.csv("C:/Users/asp1d/Desktop/CommonLists/mouse_biomart.txt", 
                    header = T, sep = "\t", as.is = T)
biomart <- biomart[!duplicated(biomart$Gene.stable.ID), ]

d0d1 <- intersect(row.names(day0), row.names(day1))
biomart <- biomart[biomart$Gene.stable.ID %in% d0d1, ]


day0 <- CreateSeuratObject(day0, project = 'D0')
day1 <- CreateSeuratObject(day1, project = 'D1')
day3 <- CreateSeuratObject(day3, project = 'D3')

joint <- merge(day0, c(day1, day3), 
               add.cell.ids = c("D0", "D1", "D3"))

remove(day0, day1, day3)
gc()


joint$orig.ident <- as.factor(joint$orig.ident)
table(joint$orig.ident)

#D0     D1     D3 
#737280 737280 737280 

#this is actually a unfiltered matrix of all possible barcodes, without any cutoff
#very cumbersome for analysis
#so first we add a minimum filter

joint <- subset(joint, subset = nFeature_RNA > 200)

table(joint$orig.ident)

#D0   D1   D3 
#2385 6188 5708 

joint <- joint@assays$RNA@counts
joint <- joint[rowSums(joint) > 1, ]

jmeta <- joint@meta.data
###annotations are given as ENSMUSG instead of gene names
##so we need to do an adjustment first


biomart <- read.csv("C:/Users/asp1d/Desktop/CommonLists/mouse_biomart.txt", 
                    header = T, sep = "\t", as.is = T)
biomart <- biomart[!duplicated(biomart$Gene.stable.ID), ]

dgenes <- row.names(joint)
biomart <- biomart[biomart$Gene.stable.ID %in% dgenes, ]
dgenes <- intersect(row.names(joint), biomart$Gene.stable.ID)

joint <- joint[dgenes, ]
joint <- as.data.frame(joint)

ncol(joint) #14281

joint$Gene.stable.ID <- row.names(joint)
joint <- join(joint, biomart, by = 'Gene.stable.ID')
joint <- joint[order(rowSums(joint[1:14281]), decreasing = T), ]
joint <- joint[!duplicated(joint$Gene.name), ]
joint <- joint[!is.na(joint$Gene.name), ]
row.names(joint) <- joint$Gene.name
joint <- joint[1:14281]
joint <- joint[row.names(joint) != "", ]

joint <- CreateSeuratObject(joint, meta.data = jmeta)

joint <- PercentageFeatureSet(joint, pattern = 'mt-', 
                              col.name = 'Mitochondrial')

##ribosome gene list from miyazaki website
ribos <- readLines("C:/Users/asp1d/Desktop/CommonLists/mouse_ribogenes.txt")
ribos <- intersect(ribos, row.names(joint@assays$RNA@data))
joint <- PercentageFeatureSet(joint, features = ribos, 
                              col.name = 'Ribosomal')

joint$orig.ident <- as.factor(joint$orig.ident)
table(joint$orig.ident)

#D0   D1   D3 
#2385 6188 5708 

VlnPlot(joint, features = c("nFeature_RNA", "Mitochondrial", "Ribosomal"))

joint <- subset(joint, subset = nFeature_RNA > 500)
joint <- subset(joint, subset = nFeature_RNA < 6000)
joint <- subset(joint, subset = Mitochondrial < 10)

table(joint$orig.ident)
#D0   D1   D3 
#2058 2687 2052 

joint <- NormalizeData(joint)
joint <- FindVariableFeatures(joint, nfeatures = 4000)
joint <- ScaleData(joint)
joint <- RunPCA(joint)
joint <- RunHarmony(joint, group.by.vars = 'orig.ident', dims.use = 1:30, 
                    block.size = 0.025)
joint <- RunUMAP(joint, reduction = 'harmony', dims = 1:30)

DimPlot(joint, group.by = 'orig.ident')


FeaturePlot(joint, features = c("Csf1r", "Mrc1" , "Acp5", "Ctsk",
                                "Nfatc1", "Tnfrsf11a"), ncol = 3,
            cols = viridis(10), min.cutoff = 'q05', max.cutoff = 'q95')&NoAxes()


joint <- FindNeighbors(joint, dims = 1:2, reduction = 'umap')
joint <- FindClusters(joint)
DimPlot(joint, label = T, label.size = 4)

osteo <- subset(joint, idents = c('7', '17', '8', '0', '6', '11', '5'))
saveRDS(osteo, "osteo_seurat.rds")
saveRDS(joint, "Joint_seurat.rds")


#####
###realm of possible osteoclasts
osteo <- FindVariableFeatures(osteo, nfeatures = 4000)
osteo <- ScaleData(osteo)
osteo <- RunPCA(osteo)

osteo <- RunHarmony(osteo, group.by.vars = 'orig.ident', dims.use = 1:30, 
                    block.size = 0.025)
osteo <- RunUMAP(osteo, reduction = 'harmony', dims = 1:30)

DimPlot(osteo, group.by = 'orig.ident')
FeaturePlot(osteo, features = c("Acp5", "Ctsk", "Nfatc1", "Tnfrsf11a"))

osteo$Stash_higherlevel_clusters <- osteo$seurat_clusters
osteo <- FindNeighbors(osteo, dims = 1:2, reduction = 'umap')
osteo <- FindClusters(osteo)
DimPlot(osteo, label = T, label.size = 4)

osteo147174 <- subset(osteo, idents = c(4, 6, 7, 10, 17, 3, 13, 14, 16))
osteo147174 <- FindVariableFeatures(osteo147174, nfeatures = 4000)


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

ostmono <- V3toMonocle(osteo147174)
ostmono <- setOrderingFilter(ostmono, ordering_genes = osteo147174@assays$RNA@var.features)
ostmono <- estimateSizeFactors(ostmono)
ostmono <- estimateDispersions(ostmono)
ostmono <- reduceDimension(ostmono, reduction_method = 'DDRTree', norm_method = 'none')
ostmono  <- orderCells(ostmono)

plot_cell_trajectory(ostmono, color_by = 'Pseudotime', markers = 'Ctsk')


reactome <- readGMT("C:/Users/asp1d/Desktop/CommonLists/Pathways_Reactome.gmt")
ostavg <- AverageExpression(osteo147174)
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
row.names(dftime) <- WhichCells(osteo147174)

cortime   <- lapply(1:length(reactome), function(x) cor(ostmono$Pseudotime, dftime[[x]])) #pseudotime correlations
dualcomp  <- as.data.frame(cortime)
dualcomp <- t(dualcomp)
dualcomp <- as.data.frame(dualcomp)
colnames(dualcomp)[1] = "Correlation_Order"
row.names(dualcomp) <- names(reactome)
lenreact <- as.data.frame(t(as.data.frame(lenreact)))
dualcomp$SetSize <- lenreact$V1

write.csv(dualcomp, "pseudotime_correlations_gse145477.csv")
write.csv(dftime, "pseudotime_table_gse145477.csv")

p1 <- plot_cell_trajectory(ostmono, color_by = 'Pseudotime', show_branch_points = F)
p1$layers[[1]]$aes_params$alpha <- 0.7
p1$layers[[2]]$aes_params$alpha <- 0.7
p1$layers[[2]]$aes_params$stroke <- 1

##correlation with pseudotime
pick6 <- c("Acp5", "Atp6v0d2", "Nfatc1", "Ctsk", "Dcstamp", "Mmp9")
ostsca <- osteo147174@assays$RNA@data
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

