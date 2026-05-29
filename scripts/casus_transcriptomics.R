setwd("~/Rprogrammas/transcriptomics/casus_transcriptomics")
rm(list = ls())
if(!is.null(dev.list())) dev.off()

BiocManager::install("Rsamtools")
BiocManager::install("Rsubread")
library(Rsamtools)
library(Rsubread)

buildindex(
  basename = 'ref_homo_sapiens',
  reference = 'GCF_000001405.40_GRCh38.p14_genomic.fna',
  memory = 8000,
  indexSplit = TRUE)

align.normaal1 <- align(index = "ref_homo_sapiens", readfile1 = "SRR4785819_1_subset40k.fastq", readfile2 = "SRR4785819_2_subset40k.fastq", output_file = "normaal1.BAM")
align.normaal2 <- align(index = "ref_homo_sapiens", readfile1 = "SRR4785820_1_subset40k.fastq", readfile2 = "SRR4785820_2_subset40k.fastq", output_file = "normaal2.BAM")
align.normaal3 <- align(index = "ref_homo_sapiens", readfile1 = "SRR4785828_1_subset40k.fastq", readfile2 = "SRR4785828_2_subset40k.fastq", output_file = "normaal3.BAM")
align.normaal4 <- align(index = "ref_homo_sapiens", readfile1 = "SRR4785831_1_subset40k.fastq", readfile2 = "SRR4785831_2_subset40k.fastq", output_file = "normaal4.BAM")
align.rheuma1 <- align(index = "ref_homo_sapiens", readfile1 = "SRR4785979_1_subset40k.fastq", readfile2 = "SRR4785979_2_subset40k.fastq", output_file = "rheuma1.BAM")
align.rheuma2 <- align(index = "ref_homo_sapiens", readfile1 = "SRR4785980_1_subset40k.fastq", readfile2 = "SRR4785980_2_subset40k.fastq", output_file = "rheuma2.BAM")
align.rheuma3 <- align(index = "ref_homo_sapiens", readfile1 = "SRR4785986_1_subset40k.fastq", readfile2 = "SRR4785986_2_subset40k.fastq", output_file = "rheuma3.BAM")
align.rheuma4 <- align(index = "ref_homo_sapiens", readfile1 = "SRR4785988_1_subset40k.fastq", readfile2 = "SRR4785988_2_subset40k.fastq", output_file = "rheuma4.BAM")

# Bestandsnamen van de monsters
samples <- c('normaal1', 'normaal2', 'normaal3', 'normaal4', 'rheuma1', 'rheuma2', 'rheuma3', 'rheuma4')

# Voor elk monster: sorteer en indexeer de BAM-file
# Sorteer BAM-bestanden
lapply(samples, function(s) {sortBam(file = paste0(s, '.BAM'), destination = paste0(s, '.sorted'))
})
# Indexeer de gesorteerde BAM-file
lapply(samples, function(s) {indexBam(file = paste0(s, '.sorted.bam'))
})

# Count matrix
count_vector <- c("normaal1.BAM", "normaal2.BAM", "normaal3.BAM", "normaal4.BAM", "rheuma1.BAM", "rheuma2.BAM", "rheuma3.BAM", "rheuma4.BAM")
count_matrix <- featureCounts(
  files = count_vector,
  annot.ext = "genomic.gtf",
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE, 
  GTF.featureType = "gene",
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE
)

counts <- count_matrix$counts
head(counts)

colnames(counts) <- c("normaal1", "normaal2", "normaal3", "normaal4", "rheuma1", "rheuma2", "rheuma3", "rheuma4")
head(counts)

write.csv(counts, "rheuma_countmatrix.csv")

# statistiek en analyse
counts <- read.table("count_matrix_RA.txt", header = TRUE)
colnames(counts) <- c("normaal1", "normaal2", "normaal3", "normaal4", "rheuma1", "rheuma2", "rheuma3", "rheuma4")

BiocManager::install("DESeq2")
BiocManager::install("KEGGREST")
BiocManager::install("EnhancedVolcano")
BiocManager::install("pathview")

library(DESeq2)
library(KEGGREST)
library(EnhancedVolcano)
library(pathview)

treatment <- c("control", "control", "control", "control", "rheuma", "rheuma", "rheuma", "rheuma")
treatment_table <- data.frame(treatment)

rownames(treatment_table) <- c("normaal1", "normaal2", "normaal3", "normaal4", "rheuma1", "rheuma2", "rheuma3", "rheuma4")

# Maak DESeqDataSet aan
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = treatment_table,
                              design = ~ treatment)

# Voer analyse uit
dds <- DESeq(dds)
resultaten <- results(dds)

write.table(resultaten, file = 'Stat_Analyse.csv', row.names = TRUE, col.names = TRUE)

sum(resultaten$padj < 0.05 & resultaten$log2FoldChange > 1, na.rm = TRUE)
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange < -1, na.rm = TRUE)

# Sorteren
hoogste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = TRUE), ]
laagste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = FALSE), ]
laagste_p_waarde <- resultaten[order(resultaten$padj, decreasing = FALSE), ]

# Plotten
EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj')

# Plot opslaan
dev.copy(png, 'VolcanoplotCASUS.png', 
         width = 8,
         height = 10,
         units = 'in',
         res = 500)
dev.off()

# Go analyse voorbereiding
BiocManager::install('goseq')
BiocManager::install('geneLenDataBase')
library(goseq)
library(geneLenDataBase)
library(tidyverse)

supportedOrganisms() %>% filter(str_detect(Genome, "hg19"))

#significante data er uit halen
sigData <- as.integer(!is.na(resultaten$padj) & resultaten$padj < 0.01)
names(sigData) <- rownames(resultaten)

head(sigData)

#PWF gebruiken

pwf <- nullp(sigData, "hg19", "geneSymbol", bias.data = resultaten$padj)
goResults <- goseq(pwf, "hg19","geneSymbol", test.cats=c("GO:BP"))

#Package installeren en inladen
install.packages('GOplot')
library(GOplot)


goResults %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")

resultaten[1] <- NULL
resultaten[2:5] <- NULL

pathview(
  gene.data = resultaten,
  pathway.id = "hsa05323",  
  species = "hsa",          
  gene.idtype = "KEGG",     
)

keggLink("hsa", "pathway:hsa05323")
