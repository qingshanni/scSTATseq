setwd("D:/refT")
###read in libraries

library("plyr")
library("Seurat")
library("viridis")
library("reshape2")
library("ggplot2")
library("matrixStats")

###load prepped comparisons
scstat <- readRDS("D:/hek293_reference/scstat_cd8t.rds")
celseq2 <- readRDS("D:/hek293_reference/clean_celseq2/celseq.cd8t.raw.rds")
chromium <- readRDS("D:/hek293_reference/clean_chromium/chromium.cd8t.raw.rds")
qtzseq2 <- readRDS("D:/hek293_reference/clean_qtz2/quartz.cd8t.raw.rds")
smart2 <- readRDS("D:/hek293_reference/clean_smrt2/smrt2.cd8t.raw.rds")
smart3 <- readRDS("D:/hek293_reference/clean_smrt3/smart3.cd8t.raw.rds")


###individually filter each dataset
#scstat
VlnPlot(scstat, features = 
          c("nFeature_RNA", 'Mitochondrial', 'Ribosomal')
        , ncol = 4)

scstat.F.median <- median(scstat$nFeature_RNA)
scstat.F.mad <- 1.4826*median(abs(scstat$nFeature_RNA - scstat.F.median))
scstat.M.median <- median(scstat$Mitochondrial)
scstat.M.mad <- 1.4826*median(abs(scstat$Mitochondrial - scstat.M.median))
scstat.R.median <- median(scstat$Ribosomal)
scstat.R.mad <- 1.4826*median(abs(scstat$Ribosomal - scstat.R.median))

scstat <- subset(scstat, subset = nFeature_RNA < scstat.F.median+2.5*scstat.F.mad)
scstat <- subset(scstat, subset = nFeature_RNA > scstat.F.median-2.5*scstat.F.mad)
scstat <- subset(scstat, subset = Mitochondrial < 10)
scstat <- subset(scstat, subset = Ribosomal < scstat.R.median+2.5*scstat.R.mad)
scstat <- subset(scstat, subset = Ribosomal > scstat.R.median-2.5*scstat.R.mad)

levels(scstat$orig.ident) <- "scSTATseq"

#celseq2
VlnPlot(celseq2, features = 
          c("nFeature_RNA", 'Mitochondrial', 'Ribosomal')
        , ncol = 4, group.by = 'orig.ident')

celseq2.F.median <- median(celseq2$nFeature_RNA)
celseq2.F.mad <- 1.4826*median(abs(celseq2$nFeature_RNA - celseq2.F.median))
celseq2.M.median <- median(celseq2$Mitochondrial)
celseq2.M.mad <- 1.4826*median(abs(celseq2$Mitochondrial - celseq2.M.median))
celseq2.R.median <- median(celseq2$Ribosomal)
celseq2.R.mad <- 1.4826*median(abs(celseq2$Ribosomal - celseq2.R.median))


celseq2 <- subset(celseq2, subset = nFeature_RNA < celseq2.F.median+2.5*celseq2.F.mad)
celseq2 <- subset(celseq2, subset = nFeature_RNA > celseq2.F.median-2.5*celseq2.F.mad)
celseq2 <- subset(celseq2, subset = Mitochondrial < celseq2.M.median+2.5*celseq2.M.mad)
celseq2 <- subset(celseq2, subset = Mitochondrial > celseq2.M.median-2.5*celseq2.M.mad)
celseq2 <- subset(celseq2, subset = Ribosomal < celseq2.R.median+2.5*celseq2.R.mad)
celseq2 <- subset(celseq2, subset = Ribosomal > celseq2.R.median-2.5*celseq2.R.mad)

levels(celseq2$orig.ident) <- "CELseq2"

#chromium
VlnPlot(chromium, features = 
          c("nFeature_RNA", 'Mitochondrial', 'Ribosomal')
        , ncol = 4, group.by = 'orig.ident')

chromium.F.median <- median(chromium$nFeature_RNA)
chromium.F.mad <- 1.4826*median(abs(chromium$nFeature_RNA - chromium.F.median))
chromium.M.median <- median(chromium$Mitochondrial)
chromium.M.mad <- 1.4826*median(abs(chromium$Mitochondrial - chromium.M.median))
chromium.R.median <- median(chromium$Ribosomal)
chromium.R.mad <- 1.4826*median(abs(chromium$Ribosomal - chromium.R.median))


chromium <- subset(chromium, subset = nFeature_RNA < chromium.F.median+2.5*chromium.F.mad)
chromium <- subset(chromium, subset = nFeature_RNA > chromium.F.median-2.5*chromium.F.mad)
chromium <- subset(chromium, subset = Mitochondrial < chromium.M.median+2.5*chromium.M.mad)
chromium <- subset(chromium, subset = Mitochondrial > chromium.M.median-2.5*chromium.M.mad)
chromium <- subset(chromium, subset = Ribosomal < chromium.R.median+2.5*chromium.R.mad)
chromium <- subset(chromium, subset = Ribosomal > chromium.R.median-2.5*chromium.R.mad)

levels(chromium$orig.ident) <- "Chromium"

#qtzseq2
VlnPlot(qtzseq2, features = 
          c("nFeature_RNA", 'Mitochondrial', 'Ribosomal')
        , ncol = 4, group.by = 'orig.ident')

qtzseq2.F.median <- median(qtzseq2$nFeature_RNA)
qtzseq2.F.mad <- 1.4826*median(abs(qtzseq2$nFeature_RNA - qtzseq2.F.median))
qtzseq2.M.median <- median(qtzseq2$Mitochondrial)
qtzseq2.M.mad <- 1.4826*median(abs(qtzseq2$Mitochondrial - qtzseq2.M.median))
qtzseq2.R.median <- median(qtzseq2$Ribosomal)
qtzseq2.R.mad <- 1.4826*median(abs(qtzseq2$Ribosomal - qtzseq2.R.median))


qtzseq2 <- subset(qtzseq2, subset = nFeature_RNA < qtzseq2.F.median+2.5*qtzseq2.F.mad)
qtzseq2 <- subset(qtzseq2, subset = nFeature_RNA > qtzseq2.F.median-2.5*qtzseq2.F.mad)
qtzseq2 <- subset(qtzseq2, subset = Mitochondrial < qtzseq2.M.median+2.5*qtzseq2.M.mad)
qtzseq2 <- subset(qtzseq2, subset = Mitochondrial > qtzseq2.M.median-2.5*qtzseq2.M.mad)
qtzseq2 <- subset(qtzseq2, subset = Ribosomal < qtzseq2.R.median+2.5*qtzseq2.R.mad)
qtzseq2 <- subset(qtzseq2, subset = Ribosomal > qtzseq2.R.median-2.5*qtzseq2.R.mad)

levels(qtzseq2$orig.ident) <- "QTZseq2"

#smart2
VlnPlot(smart2, features = 
          c("nFeature_RNA", 'Mitochondrial', 'Ribosomal')
        , ncol = 4, group.by = 'orig.ident')

smart2.F.median <- median(smart2$nFeature_RNA)
smart2.F.mad <- 1.4826*median(abs(smart2$nFeature_RNA - smart2.F.median))
smart2.M.median <- median(smart2$Mitochondrial)
smart2.M.mad <- 1.4826*median(abs(smart2$Mitochondrial - smart2.M.median))
smart2.R.median <- median(smart2$Ribosomal)
smart2.R.mad <- 1.4826*median(abs(smart2$Ribosomal - smart2.R.median))


smart2 <- subset(smart2, subset = nFeature_RNA < smart2.F.median+2.5*smart2.F.mad)
smart2 <- subset(smart2, subset = nFeature_RNA > smart2.F.median-2.5*smart2.F.mad)
smart2 <- subset(smart2, subset = Mitochondrial < smart2.M.median+2.5*smart2.M.mad)
smart2 <- subset(smart2, subset = Mitochondrial > smart2.M.median-2.5*smart2.M.mad)
smart2 <- subset(smart2, subset = Ribosomal < smart2.R.median+2.5*smart2.R.mad)
smart2 <- subset(smart2, subset = Ribosomal > smart2.R.median-2.5*smart2.R.mad)

levels(smart2$orig.ident) <- "SMARTseq2"

#smart3
VlnPlot(smart3, features = 
          c("nFeature_RNA", 'Mitochondrial', 'Ribosomal')
        , ncol = 4, group.by = 'orig.ident')

smart3.F.median <- median(smart3$nFeature_RNA)
smart3.F.mad <- 1.4826*median(abs(smart3$nFeature_RNA - smart3.F.median))
smart3.M.median <- median(smart3$Mitochondrial)
smart3.M.mad <- 1.4826*median(abs(smart3$Mitochondrial - smart3.M.median))
smart3.R.median <- median(smart3$Ribosomal)
smart3.R.mad <- 1.4826*median(abs(smart3$Ribosomal - smart3.R.median))


smart3 <- subset(smart3, subset = nFeature_RNA < smart3.F.median+2.5*smart3.F.mad)
smart3 <- subset(smart3, subset = nFeature_RNA > smart3.F.median-2.5*smart3.F.mad)
smart3 <- subset(smart3, subset = Mitochondrial < smart3.M.median+2.5*smart3.M.mad)
smart3 <- subset(smart3, subset = Mitochondrial > smart3.M.median-2.5*smart3.M.mad)
smart3 <- subset(smart3, subset = Ribosomal < smart3.R.median+2.5*smart3.R.mad)
smart3 <- subset(smart3, subset = Ribosomal > smart3.R.median-2.5*smart3.R.mad)

levels(smart3$orig.ident) <- "SMARTseq3"


########
####
###
##
#Now we combine dataframes together to draw our figure
##
###
####
#######

joint <- merge(celseq2, c(chromium, qtzseq2, scstat, smart2, smart3))
VlnPlot(joint, features = 'nFeature_RNA', group.by = 'orig.ident')

jmeta <- joint@meta.data
colnames(jmeta)[1] <- "Library"
colnames(jmeta)[3] <- "Number of Genes Detected"
jmeta <- jmeta[1:6]

pdf("violinxboxplot_nfeature_RNA.pdf", height = 6, width = 11)
ggplot(jmeta, aes(Library, `Number of Genes Detected`, fill = Library))+
  geom_boxplot(alpha = 0.7, notch = T, width = 0.2)+
  geom_violin(alpha = 0.4)+
  theme_classic()+ylim(0,10000)+
  theme(text = element_text(size = 20, face = 'bold'))
dev.off()

pdf("violinxboxplot_Mitochondrial.pdf", height = 6, width = 10)
ggplot(jmeta, aes(Library, Mitochondrial, fill = Library))+
  geom_boxplot(alpha = 0.7, notch = T, width = 0.2)+
  geom_violin(alpha = 0.4)+
  theme_classic()+ylim(0,60)+
  theme(text = element_text(size = 16, face = 'bold'))
dev.off()

pdf("violinxboxplot_Ribosomal.pdf", height = 6, width = 10)
ggplot(jmeta, aes(Library, Ribosomal, fill = Library))+
  geom_boxplot(alpha = 0.7, notch = T, width = 0.2)+
  geom_violin(alpha = 0.4)+
  theme_classic()+ylim(0,50)+
  theme(text = element_text(size = 16, face = 'bold'))
dev.off()


####
AverageDetectionRate_V3 <- function(
    object,
    thresh.min = 0
) {
  ident.use <- object@active.ident
  data.all <- data.frame(row.names = rownames(x = object@assays$RNA@data))
  for (i in sort(x = unique(x = ident.use))) {
    temp.cells <- WhichCells(object = object, ident = i)
    data.temp <- apply(
      X = object@assays$RNA@data[, temp.cells],
      MARGIN = 1,
      FUN = function(x) {
        return(sum(x > thresh.min)/length(x = x))
      }
    )
    data.all <- cbind(data.all, data.temp)
    colnames(x = data.all)[ncol(x = data.all)] <- i
  }
  colnames(x = data.all) <- sort(x = unique(x = ident.use))
  return(data.all)
}

####calculate detection curves
joint <- SetIdent(joint, value = 'orig.ident')

jdet <- AverageDetectionRate_V3(joint)
write.csv(jdet, "tcell_jointdetctionrate.csv")


##finangle
jdetlong <- read.csv("D:/hek293_reference/hCD8T_finalfig/long.csv", 
                     header =T, row.names = 1)
jdetlong <- melt(jdetlong, id.vars = 'Rank')
colnames(jdetlong)[2] <- 'Library'
colnames(jdetlong)[3] <- 'DetectionRate'
jdetlong <- jdetlong[jdetlong$Rank < 20001, ]

pdf("dropoutdistribution_tcell.pdf", height = 6, width = 10)
ggplot(jdetlong, aes(Rank, DetectionRate, color = Library))+
  geom_line(size = 1.7, alpha = 0.7)+theme_classic()+
  theme(text = element_text(size = 24, face = 'bold'))
dev.off()


###now calculate the EPV
jmeta.out <- joint@meta.data
joint.count <- joint@assays$RNA@data
joint.count <- as.matrix(joint.count)
libraries <- levels(as.factor(jmeta.out$orig.ident))
spl.lib <- lapply(libraries, function(x)
  row.names(jmeta.out[jmeta.out$orig.ident == x, ]))

lib.epv <- lapply(spl.lib, function(x)
  (rowSds(joint.count[,x])/rowMeans(joint.count[,x])) -  ##CV
    (sqrt(rowMeans(joint.count[,x]))/rowMeans(joint.count[,x])) #poisson
  )
names(lib.epv) <- libraries
lib.epv <- as.data.frame(lib.epv)
write.csv(lib.epv, 'hCD8T_library_epvvalues.csv')

lib.epv$Mean <- rowMeans(lib.epv)
lib.epv <- lib.epv[lib.epv$Mean != 'NaN', ]
lib.epv <- lib.epv[1:6]

lib.epv$Gene <- row.names(lib.epv)
lib.epv <- melt(lib.epv, id.vars = 'Gene')

colnames(lib.epv)[1] <- 'Gene_Count'
colnames(lib.epv)[2] <- 'Library'
colnames(lib.epv)[3] <- 'EPV'

##option A
ggplot(lib.epv, aes(Library, EPV, fill = Library))+
  geom_boxplot(alpha = 0.7, notch = T, width = 0.2)+
  geom_violin(alpha = 0.4)+
  theme_classic()+ylim(-10,10)+
  theme(text = element_text(size = 16, face = 'bold'))


###option B
pdf("EPVdistribution_tcell.pdf", height = 6, width = 10)
ggplot(lib.epv, aes(EPV, fill = Library))+
  geom_freqpoly(aes(color = Library), size = 1.25,
                alpha = 0.7, position= 'identity', bins = 500)+
  geom_histogram(alpha = 0.2, position= 'identity', bins = 500)+
  theme_classic()+xlim(-2,2)+labs(y = 'Gene Count')+
  theme(text = element_text(size = 24, face = 'bold'))


####