###############################################################################
# Transcriptomics Reumatoïde Artritis
# Differentiële Expressie, GO Enrichment en KEGG Pathway Analysis
# Reumatoïde Artritis vs Gezonde Vrouwen
###############################################################################

###############################################################################
# 1. SET WORKING DIRECTORY
###############################################################################

setwd("C:/Users/mette/OneDrive - NHL Stenden/BML2/Periode 4/Casus Reuma Artritis/Data_RA_raw/Data_RA_raw/")
getwd()

###############################################################################
# 2. READ ALIGNMENT
###############################################################################

# Packages installeren
BiocManager::install("Rsubread")

# Package inladen 
library(Rsubread)

###############################################################################
# 2.1 Build Reference Genome Index
###############################################################################

buildindex(
  basename = "ref_reuma",
  reference = "GCF_000001405.40_GRCh38.p14_genomic.fna",
  memory = 4000,
  indexSplit = TRUE)

###############################################################################
# 2.2 Align Reads to Reference Genome
###############################################################################

# Controle samples
align.cntrl1 <- align(
  index = "ref_reuma",
  readfile1 = "SRR4785819_1_subset40k.fastq",
  readfile2 = "SRR4785819_2_subset40k.fastq",
  output_file = "cntrl1.BAM")

align.cntrl2 <- align(
  index = "ref_reuma",
  readfile1 = "SRR4785820_1_subset40k.fastq",
  readfile2 = "SRR4785820_2_subset40k.fastq",
  output_file = "cntrl2.BAM")

align.cntrl3 <- align(
  index = "ref_reuma",
  readfile1 = "SRR4785828_1_subset40k.fastq",
  readfile2 = "SRR4785828_2_subset40k.fastq",
  output_file = "cntrl3.BAM")

align.cntrl4 <- align(
  index = "ref_reuma",
  readfile1 = "SRR4785831_1_subset40k.fastq",
  readfile2 = "SRR4785831_2_subset40k.fastq",
  output_file = "cntrl4.BAM")

# Reumatoïde artritis samples
align.reuma1 <- align(
  index = "ref_reuma",
  readfile1 = "SRR4785979_1_subset40k.fastq",
  readfile2 = "SRR4785979_2_subset40k.fastq",
  output_file = "reuma1.BAM")

align.reuma2 <- align(
  index = "ref_reuma",
  readfile1 = "SRR4785980_1_subset40k.fastq",
  readfile2 = "SRR4785980_2_subset40k.fastq",
  output_file = "reuma2.BAM")

align.reuma3 <- align(
  index = "ref_reuma",
  readfile1 = "SRR4785986_1_subset40k.fastq",
  readfile2 = "SRR4785986_2_subset40k.fastq",
  output_file = "reuma3.BAM")

align.reuma4 <- align(
  index = "ref_reuma",
  readfile1 = "SRR4785988_1_subset40k.fastq",
  readfile2 = "SRR4785988_2_subset40k.fastq",
  output_file = "reuma4.BAM")

###############################################################################
# 3. MAKEN VAN COUNT MATRIX
###############################################################################
library(Rsubread)

samples <- c(
  "cntrl1.BAM", "cntrl2.BAM", "cntrl3.BAM", "cntrl4.BAM",
  "reuma1.BAM", "reuma2.BAM", "reuma3.BAM", "reuma4.BAM")

count_matrix <- featureCounts(
  files = samples,
  annot.ext = "genomic.gtf",
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE,
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE)

counts <- count_matrix$counts

colnames(counts) <- c(
  "cntrl1", "cntrl2", "cntrl3", "cntrl4",
  "reuma1", "reuma2", "reuma3", "reuma4")

write.csv(counts, "reuma_countmatrix.csv")

###############################################################################
# 4. ANALYSE VAN DIFFERENTIËLE GENEXPRESSIE (DESeq2)
###############################################################################

BiocManager::install("DESeq2")
BiocManager::install("KEGGREST")
BiocManager::install("EnhancedVolcano")
BiocManager::install("pathview")

library(DESeq2)
library(KEGGREST)
library(EnhancedVolcano)
library(pathview)

count_matrix_RA <- read.table(
  "count_matrix_RA.txt",
  row.names = 1)

treatment <- c(
  "control", "control", "control", "control",
  "patiënt", "patiënt", "patiënt", "patiënt")

treatment_table <- data.frame(treatment)

rownames(treatment_table) <- c(
  "cntrl1", "cntrl2", "cntrl3", "cntrl4",
  "reuma1", "reuma2", "reuma3", "reuma4")

colnames(count_matrix_RA) <- rownames(treatment_table)

# 4.4 DESeqDataSet aanmaken
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_RA,
  colData = treatment_table,
  design = ~ treatment)

# 4.5 Analyse uitvoeren 
dds <- DESeq(dds)

resultaten <- results(dds)

write.table(
  resultaten,
  file = "Resultaten_Casus.csv",
  row.names = TRUE,
  col.names = TRUE)

###############################################################################
# 5. VOLCANO PLOT
###############################################################################

EnhancedVolcano(
  resultaten,
  lab = rownames(resultaten),
  x = "log2FoldChange",
  y = "padj")

###############################################################################
# 6. GENE ONTOLOGY (GO) ENRICHMENT ANALYSIS
###############################################################################

BiocManager::install(c(
  "clusterProfiler",
  "DOSE",
  "enrichplot",
  "AnnotationDbi",
  "GO.db",
  "biomaRt",
  "org.Hs.eg.db"))

install.packages("ggplot2")

library(clusterProfiler)
library(org.Hs.eg.db)

sig_genen <- subset(
  resultaten,
  padj < 0.05 &
    abs(log2FoldChange) > 1)

sig_genes <- rownames(sig_genen)

gene.df <- bitr(
  sig_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db)

ego_bp <- enrichGO(
  gene = gene.df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05)

###############################################################################
# 6.1 GO Visualization
###############################################################################

dotplot(
  ego_bp,
  showCategory = 15,
  font.size = 7)

png(
  "GO_dotplot_RA.png",
  width = 2000,
  height = 1500,
  res = 300)

dotplot(
  ego_bp,
  showCategory = 20)

dev.off()

###############################################################################
# 7. KEGG PATHWAY ANALYSIS
###############################################################################

entrez_ids <- unique(
  na.omit(gene.df$ENTREZID))

kegg_result <- enrichKEGG(
  gene = entrez_ids,
  organism = "hsa")

dotplot(
  kegg_result,
  showCategory = 10)

barplot(
  kegg_result,
  showCategory = 10)

###############################################################################
# 8. PATHVIEW VISUALIZATION
###############################################################################

entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = rownames(sig_genen),
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first")

entrez_ids <- na.omit(entrez_ids)

gene_data <- sig_genen$log2FoldChange
names(gene_data) <- entrez_ids

gene_data <- gene_data[
  !is.na(names(gene_data))
]

kegg_df <- as.data.frame(kegg_result)

top_pathway <- kegg_df$ID[1]

pv_out <- pathview(
  gene.data = gene_data,
  pathway.id = top_pathway,
  species = "hsa")

#inzoomen pathway
View(kegg_df)

kegg_tabel <- kegg_df[, c(
  "ID",
  "Description",
  "GeneRatio",
  "Count",
  "p.adjust")]

write.csv(kegg_tabel,
          "kegg_df.csv",
          row.names = TRUE)

grep("Rheumatoid", kegg_df$Description, value = TRUE)

grep("TNF", kegg_df$Description, value = TRUE)

grep("NOD", kegg_df$Description, value = TRUE)

pathview(
  gene.data = gene_data,
  pathway.id = "hsa05323",
  species = "hsa")

