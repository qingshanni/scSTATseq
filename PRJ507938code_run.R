prj507938 <- read.csv("E:/Osteoclast_singlecell/PRJNA507938/TPM_salmon_all_gene_symbol.txt", 
                      header = T, sep = "\t", row.names = 1)
srameta <- read.csv("E:/Osteoclast_singlecell/PRJNA507938/SraRunTable.txt", 
                    header = T)

mm10 <- read.csv("D:/scSTATseq_paper/comparisons_gene_stats/mouse_genome.txt", 
                 header = T, as.is = T, sep = "\t")
mm10 <- mm10[!duplicated(mm10$Gene.name), ]
nrow(mm10)
detg <- intersect(mm10$Gene.name, row.names(prj507938))

prj507938 <- prj507938[detg, ]

detcell <- intersect(names(prj507938), srameta$Run)
srameta <- srameta[srameta$Run %in% detcell, ]

##seurat
prj507938 <- CreateSeuratObject(prj507938)
prjmeta <- prj507938@meta.data
row.names(srameta) <- srameta$Run
prjmeta$Run <- row.names(prjmeta)
prjmeta <- join(prjmeta, srameta, by = 'Run')
row.names(prjmeta) <- prjmeta$Run

prj507938 <- AddMetaData(prj507938, metadata = prjmeta)
prj507938$SortType <- paste0(prj507938$Tissue, "_", prj507938$Phenotype)

prj507938 <- PercentageFeatureSet(prj507938, pattern = 'mt-', 
                                  col.name = 'Mitochondrial')

rpls <- grep(row.names(prj507938@assays$RNA@counts), pattern ="Rpl", value = T)
rpss <- grep(row.names(prj507938@assays$RNA@counts), pattern ="Rps", value = T)

prj507938 <- PercentageFeatureSet(prj507938, features = union(rpls, rpss), 
                                  col.name = 'Ribosomal')


##now filter
prj507938 <- subset(prj507938, subset = Mitochondrial < 10)
prj507938 <- subset(prj507938, subset = nFeature_RNA > 2000)
prj507938 <- subset(prj507938, subset = nFeature_RNA < 12000)


prj507938 <- NormalizeData(prj507938)
prj507938 <- FindVariableFeatures(prj507938, nfeatures = 4000)
prj507938 <- ScaleData(prj507938)
prj507938 <- RunPCA(prj507938)
prj507938 <- RunUMAP(prj507938, reduction = 'pca', dims = 1:30)

DimPlot(prj507938, group.by = 'orig.ident')


FeaturePlot(prj507938, features = c("Csf1r", "Mrc1" , "Acp5", "Ctsk",
                                "Nfatc1", "Tnfrsf11a"), ncol = 3,
            cols = viridis(10), min.cutoff = 'q05', max.cutoff = 'q95')&NoAxes()


DimPlot(prj507938, group.by = 'SortType')

pdf("umap_sorted_populations.pdf", height = 5, width = 7)
DimPlot(prj507938, group.by = 'SortType', pt.size = 2)&NoAxes()
dev.off()

prj507938 <- FindNeighbors(prj507938, dims = 1:2, reduction = 'umap')
prj507938 <- FindClusters(prj507938)
DimPlot(prj507938, label = T, label.size = 4)

osteo <- subset(prj507938, idents = c('1', '2', '4', '5'))
osteo <- FindVariableFeatures(osteo, nfeatures = 4000)

saveRDS(osteo, "osteoclast_precursor_seurat.rds")
saveRDS(prj507938, "full_seurat.rds")


###
ostmono <- V3toMonocle(osteo)
ostmono <- setOrderingFilter(ostmono, ordering_genes = osteo@assays$RNA@var.features)
ostmono <- estimateSizeFactors(ostmono)
ostmono <- estimateDispersions(ostmono)
ostmono <- reduceDimension(ostmono, reduction_method = 'DDRTree', norm_method = 'none')
ostmono  <- orderCells(ostmono)

plot_cell_trajectory(ostmono, color_by = 'Pseudotime', markers = 'Ctsk')


reactome <- readGMT("C:/Users/asp1d/Desktop/CommonLists/Pathways_Reactome.gmt")
ostavg <- AverageExpression(osteo)
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
row.names(dftime) <- WhichCells(osteo)

cortime   <- lapply(1:length(reactome), function(x) cor(ostmono$Pseudotime, dftime[[x]])) #pseudotime correlations
dualcomp  <- as.data.frame(cortime)
dualcomp <- t(dualcomp)
dualcomp <- as.data.frame(dualcomp)
colnames(dualcomp)[1] = "Correlation_Order"
row.names(dualcomp) <- names(reactome)
lenreact <- as.data.frame(t(as.data.frame(lenreact)))
dualcomp$SetSize <- lenreact$V1

write.csv(dualcomp, "pseudotime_correlations_prj507938.csv")
write.csv(dftime, "pseudotime_table_prj507938.csv")

###
p1 <- plot_cell_trajectory(ostmono, color_by = 'SortType', show_branch_points = F)
p1$layers[[1]]$aes_params$alpha <- 0.7
p1$layers[[2]]$aes_params$alpha <- 0.7
p1$layers[[2]]$aes_params$stroke <- 1