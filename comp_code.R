library(plyr)

grch38 <- read.csv("D:/scSTATseq_paper/comparisons_gene_stats/human_genome.txt", 
                   header = T, as.is  =T, sep = "\t")
grch38 <- grch38[!duplicated(grch38$Gene.name), ]
nrow(grch38)
#[1] 40221

mm10 <- read.csv("D:/scSTATseq_paper/comparisons_gene_stats/mouse_genome.txt", 
                 header = T, as.is = T, sep = "\t")
mm10 <- mm10[!duplicated(mm10$Gene.name), ]
nrow(mm10)
#[1] 53174

scstat <- read.csv("D:/scSTATseq_paper/comparisons_gene_stats/raw264_scstat_detvexp.csv", 
                   header = T, as.is = T)

smrt2 <- read.csv("D:/scSTATseq_paper/comparisons_gene_stats/raw264_smart2_detvexp.csv", 
                  header = T, as.is = T)

mutgenes <- intersect(scstat$X, smrt2$X)
mutgenes <- intersect(mutgenes, mm10$Gene.name)

mm10 <- mm10[mm10$Gene.name %in% mutgenes,  ]
scstat <- scstat[scstat$X %in% mutgenes, ]
smrt2 <- smrt2[smrt2$X %in% mutgenes, ]

colnames(scstat)[1] <- "Gene.name"
colnames(scstat)[2] <- "Count_scstat"
colnames(scstat)[3] <- "Detection_scstat"
colnames(scstat)[4] <- "Expr_scstat"

colnames(smrt2)[1] <- "Gene.name"
colnames(smrt2)[2] <- "Detection_smrt2"
colnames(smrt2)[3] <- "Count_smrt2"
colnames(smrt2)[4] <- "Expr_smrt2"

m1 <- join(mm10, scstat, by = 'Gene.name')
m1 <- join(m1, smrt2, by = 'Gene.name')

colnames(m1)[5] <- "GC_Content"
colnames(m1)[10] <- "Transcript_Length"
colnames(m1)[6] <- "Transcript_Count"

ggplot(m1, aes(Detection_scstat, GC_Content))+
  geom_point(shape = 21, fill = 'salmon', alpha = 0.5, size = 1.5)+
  geom_density_2d_filled(alpha = 0.3)+
  theme_classic()

#pdf("gene_GC_content_detection.pdf", height = 5, width = 8)
ggplot(m1, aes(Detection_scstat, GC_Content))+
  geom_point(shape = 21, fill = 'salmon', alpha = 0.5, size = 1.5)+
  geom_density_2d_filled(alpha = 0.3)+
  theme_classic()+theme(text = element_text(size = 14, face = 'bold'))


ggplot(m1, aes(Gene.type, Detection_smrt2))+
  geom_violin(fill = 'skyblue', alpha = 0.8)+
  geom_boxplot(outlier.size = 0, varwidth = T, alpha =  0.5)