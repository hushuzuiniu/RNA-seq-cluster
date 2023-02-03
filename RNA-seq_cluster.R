# library the packages we need
library(data.table)
library(Seurat)
library(readxl)
library(org.Hs.eg.db)
library(annotables)
library(tidyverse)
library(tidyr)
library(rtracklayer)
library(biomaRt)
library(pheatmap)
library(topGO)
library(GO.db)
library(biomaRt)
library(Rgraphviz)
###################### tsv version ########################
# read the tsv version counts file
tsv_matrix<- fread("tsv-datatype-0.csv",sep=",")
sample_info<-read.csv("all_sample_info.txt",sep = "\t")

# process counts file, read in Seurat
tsv_matrix<-as.matrix(tsv_matrix)
rownames(tsv_matrix)<-tsv_matrix[,1]
tsv_matrix<-tsv_matrix[,-1]

# create Seurat object
data_seurat <- CreateSeuratObject(counts = tsv_matrix,project = "mydata")
# Normalization???????????
# data_seurat<-NormalizeData(data_seurat)

# process sample_info file, add sample_info to meta.data
sample_info <- data.frame(sample_info)
rownames(sample_info)<-sample_info[,1]
data_seurat <- AddMetaData(data_seurat,sample_info)


Idents(data_seurat) <- data_seurat@meta.data$cluster
data_seurat.markers <- FindAllMarkers(data_seurat, 
                                      only.pos = TRUE, 
                                      min.pct = 0.05, 
                                      logfc.threshold = 0.25)

# write.table(nk_cells_anno_cluster.markers, 
#             file = paste("../results/12_nk_cell_analyses/Righton_singlecellAnalysis_NK_Cells_cluster",
#                          nk_cell_anno_cell_type_unique[i],"tissue_type.txt",sep="."), 
#             sep = '\t',quote = F,col.names = T,row.names = F)

# get DEG counts matrix 
counts_DEG<-data_seurat.markers[which(data_seurat.markers$p_val_adj<0.01),]
counts_DEG<-data_seurat.markers[which(abs(data_seurat.markers$avg_log2FC)>2),]
counts_DEG<-tsv_matrix[which(counts_DEG$gene %in% rownames(tsv_matrix)),]

# make plot heatmap DEG counts matrix
counts_DEG<-data.frame(counts_DEG)
counts_DEG1<-counts_DEG
counts_DEG<-as.data.frame(lapply(counts_DEG,as.numeric))
rownames(counts_DEG)<-rownames(counts_DEG1)
annotation_col1<-sample_info[c('cluster')]

as.data.frame(scale(counts_DEG))

pheatmap(as.data.frame(scale(counts_DEG)),
         show_rownames = FALSE,
         show_colnames = FALSE,
         cluster_cols = TRUE, 
         cluster_rows = FALSE,
         annotation_col = annotation_col1,
         scale = "row")


#######################excel-essembl matrix version#################
essembl_matrix<- fread("excel-essembl-gene_datatype_0.csv",sep=",")
essembl_matrix<-as.matrix(essembl_matrix)
essembl_matrix<-essembl_matrix[-1,]
rownames(essembl_matrix)<-essembl_matrix[,1]
essembl_matrix<-essembl_matrix[,-1]

essembl_matrix_seurat <- CreateSeuratObject(counts = essembl_matrix,project = "essembl_matrix")
essembl_matrix_seurat <- AddMetaData(essembl_matrix_seurat,sample_info)

Idents(essembl_matrix_seurat) <- essembl_matrix_seurat@meta.data$cluster
essembl_matrix_seurat.markers <- FindAllMarkers(essembl_matrix_seurat, 
                                                only.pos = TRUE, 
                                                min.pct = 0.05, 
                                                logfc.threshold = 0.25)

essembl_matrix_DEG<-essembl_matrix_seurat.markers[which(essembl_matrix_seurat.markers$p_val_adj<0.01),]
essembl_matrix_DEG<-essembl_matrix_seurat.markers[which(abs(essembl_matrix_seurat.markers$avg_log2FC)>2),]
essembl_matrix_DEG<-essembl_matrix[which(essembl_matrix_DEG$gene %in% rownames(essembl_matrix)),]
essembl_matrix_DEG<-data.frame(essembl_matrix_DEG)
essembl_matrix_DEG1<-essembl_matrix_DEG
essembl_matrix_DEG<-as.data.frame(lapply(essembl_matrix_DEG,as.numeric))
rownames(essembl_matrix_DEG)<-rownames(essembl_matrix_DEG1)
test<-as.data.frame(scale(essembl_matrix_DEG))

pheatmap(
  test,  # matrix of counts
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  legend = FALSE,
  cluster_cols = TRUE, # change to TRUE to get a dendrogram of samples
  cluster_rows = FALSE,
  annotation_col = annotation_col1
)


#######################excel-gene matrix version#################
gene_matrix<- fread("excel-gene-datatype0.csv",sep=",")
gene_matrix<-as.matrix(gene_matrix)
gene_matrix<-gene_matrix[-1,]
rownames(gene_matrix)<-gene_matrix[,1]
gene_matrix<-gene_matrix[,-1]

gene_matrix_seurat <- CreateSeuratObject(counts = gene_matrix,project = "gene_matrix")
gene_matrix_seurat <- AddMetaData(gene_matrix_seurat,sample_info)

Idents(gene_matrix_seurat) <- gene_matrix_seurat@meta.data$cluster
gene_matrix_seurat.markers <- FindAllMarkers(gene_matrix_seurat, 
                                             only.pos = TRUE, 
                                             min.pct = 0.05, 
                                             logfc.threshold = 0.25)


gene_matrix_DEG<-gene_matrix_seurat.markers[which(gene_matrix_seurat.markers$p_val_adj<0.01),]
gene_matrix_DEG<-gene_matrix_seurat.markers[which(abs(gene_matrix_seurat.markers$avg_log2FC)>2),]
gene_matrix_DEG<-gene_matrix[which(gene_matrix_DEG$gene %in% rownames(gene_matrix)),]
gene_matrix_DEG<-data.frame(gene_matrix_DEG)
gene_matrix_DEG1<-gene_matrix_DEG
gene_matrix_DEG<-as.data.frame(lapply(gene_matrix_DEG,as.numeric))
rownames(gene_matrix_DEG)<-rownames(gene_matrix_DEG1)
as.data.frame(scale(gene_matrix_DEG))

pheatmap(
  as.data.frame(scale(gene_matrix_DEG)),  # matrix of counts
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  legend = FALSE,
  cluster_cols = TRUE, # change to TRUE to get a dendrogram of samples
  cluster_rows = FALSE,
  annotation_col = annotation_col1
)



##################using RNA-seq method for cluster#############################
library(limma)
library(edgeR)
# read files
gene_matrix<- fread("excel-gene-datatype0.csv",sep=",")
sample_info<-read.csv("all_sample_info.txt",sep = "\t")
# using gene matrix do RNA-seq DEG analysis
gene_matrix<-as.matrix(gene_matrix)
rownames(gene_matrix)<-gene_matrix[,1]
gene_matrix<-gene_matrix[-1,-1]
# convert character matrix to numeric matrix
test<-gene_matrix
gene_matrix <- matrix(as.numeric(gene_matrix),    
                      ncol = ncol(gene_matrix))
rownames(gene_matrix)<-rownames(test)
colnames(gene_matrix)<-colnames(test)

# filter sample information
sample_info<-sample_info[which(sample_info$sample_id %in% colnames(gene_matrix)),]
rownames(sample_info)<-sample_info$sample_id

# build the DGEList object
dgList_raw<- DGEList(
  counts = gene_matrix, # counts data
  samples = sample_info, # sample data
  group = sample_info$cluster # specify experimental groups
)

# normalisation
dgList <- calcNormFactors(dgList_raw)

# create design matrix
design <- model.matrix(~ 0+ group, data = dgList$samples)

# significance tests
dgGlm <- estimateDisp(dgList, design, robust = TRUE)
fit <- glmQLFit(dgGlm, design, robust = TRUE)


# get DEG 
et <- exactTest(dgGlm)
et <- topTags(et,n=1000000)
et <- as.data.frame(et)
etSig <- et[which(et$PValue < 0.05 & abs(et$logFC) > 1),]

# get DEG counts matrix to plot heatmap
DEG<-gene_matrix[which(rownames(etSig) %in% rownames(gene_matrix)),]
annotation_col<-dgGlm$samples[c('cluster')]
pheatmap(
  DEG,  # matrix of counts
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  legend = FALSE,
  cluster_cols = TRUE, # change to TRUE to get a dendrogram of samples
  cluster_rows = FALSE,
  scale = 'column',
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
  #annotation_col = annotation_col
)



known_genes<-read.csv('sjdu.Righton.Expression.conclusion.classfication.byGene.xls',sep = "\t")
known_genes<-known_genes[-1,]
genes <- known_genes$Gene
ensembl_to_symbol<-read.csv("ensembl_to_symbol.csv",row.names =1)
genes<-data.frame(genes)
genes<-merge(genes,ensembl_to_symbol,by.x="genes",by.y="gene_id")
genes<-unique(genes[,-2])

selected_DEG<-DEG[which(rownames(DEG) %in% genes$symbol),]


#############
id_conversion <- bitr(unique(rownames(etSig)), fromType = "SYMBOL",
                      toType = c( "ENTREZID"),
                      OrgDb = org.Hs.eg.db)
deg_for_cluster <- merge(etSig,id_conversion,
                         by.x = 'row.names',by.y='SYMBOL')

geneList=etSig$logFC
names(geneList)=rownames(etSig)
install.packages("R.utils")
R.utils::setOption("clusterProfiler.download.method",'auto') 

kegg_test <- enrichKEGG(
  gene = deg_for_cluster$ENTREZID,
  keyType = "kegg",
  organism  = "hsa",
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff  = 0.05
)

GO_test <- enrichGO(
  gene = deg_for_cluster$ENTREZID,
  ont  = "ALL",
  OrgDb = "org.Hs.eg.db",
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05,
  readable = TRUE,
  minGSSize = 10,
  maxGSSize = 500,
)

library(cowplot)
p1 <- barplot(kegg_test, showCategory = 20)
p2 <- dotplot(kegg_test, showCategory = 20)
pdf('kegg_plot.pdf',width = 16,height = 8)
plot_grid(p1,p2,ncol = 2)
dev.off()
pdf('kegg_heatmap.pdf',width = 16,height = 8)
heatplot(kegg_test, foldChange = geneList)
dev.off()
pdf('kegg_cnet.pdf',width = 16,height = 8)
cnetplot(kegg_test, showCategory = 5)
dev.off()


p3 <- barplot(GO_test, showCategory = 20)
p4 <- dotplot(GO_test, showCategory = 20)
pdf('GO_plot.pdf',width = 16,height = 8)
plot_grid(p3,p4,ncol = 2)
dev.off()
pdf('GO_heatmap.pdf',width = 16,height = 8)
heatplot(GO_test, foldChange = geneList)
dev.off()
options(ggrepel.max.overlaps = Inf)

cnetplot(
  GO_test,
  showCategory = 5,
  color_category = "red",
  color_gene = "blue")

x2<-pairwise_termsim(GO_test)
emapplot(x2)


GO_results<-GO_test@result[,c(1,3,6,10)]
rownames(GO_results)<-c(1:793)
colnames(GO_results)<-c("Category","Pathway","pvalue",'Count')
GO_results$log10<-log10(GO_results$pvalue)*(-1)
GO_results<-GO_results[order(GO_results$pvalue),]
GO_results_plot<-GO_results[1:30,]

pdf('GO_dotplot.pdf',width = 18,height = 9)
p5<-ggplot(data = GO_results_plot, mapping = aes(x = log10, y = Pathway)) +
  geom_point(aes(size=Count,color = -1*log10(pvalue)))+
  scale_color_gradient(low = "blue",high = "red")+
  labs(color=expression(-log[10](Pvalue),size = "Count"))+
  xlab(expression(-log[10](Pvalue)))+facet_grid(Category~.,scales = "free")
p6<-ggplot(data = GO_results_plot, mapping = aes(x=log10, y =Pathway,fill = Category)) +
  geom_bar(stat = "identity", width = 0.5, position = position_dodge(0.8),)+
  xlab(expression(-log[10](Pvalue)))+facet_grid(Category~.,scales = "free")+
  scale_fill_manual(values = c("#669933", "#FFCC66","coral"),name = "Category", 
                    labels = c("Biological Process", "Cell Component", "Molecular Function"))
plot_grid(p5,p6,ncol = 2)
dev.off()

##############################################
GO<-read.table("GO_pathway.txt",sep="\t")
colnames(GO)<-GO[1,]
GO<-GO[-1,]   
GO<-cSplit(GO,"Term",sep="~")

GOPlot<-data.frame(GO$Category,GO$Term_2,GO$PValue,GO$Count)
colnames(GOPlot)[1:4]<-c("Category","Pathway","pvalue","Count")
GOPlot$pvalue<-as.numeric(as.character(GOPlot$pvalue))
GOPlot$Count<-as.numeric(as.character(GOPlot$Count))
GOplot2<-GOPlot[order(GOPlot$pvalue),]
GOplot3<-GOplot2[1:30,]

GOplot3$log10<-log10(GOplot3$pvalue)*(-1)

pdf('GO_dotplot1.pdf',width = 16,height = 9)
ggplot(data = GOplot3, mapping = aes(x = log10, y = Pathway)) +
  geom_point(aes(size=Count,color = -1*log10(pvalue)))+
  scale_color_gradient(low = "blue",high = "red")+
  labs(color=expression(-log[10](Pvalue),size = "Count"))+
  xlab(expression(-log[10](Pvalue)))+facet_grid(Category~.,scales = "free")
dev.off()

KEGG<-read.table("kegg_pathway.txt",sep="\t")
colnames(KEGG)<-KEGG[1,]
KEGG<-KEGG[-1,]   
KEGG<-cSplit(KEGG,"Term",sep=":")

KEGGPlot<-data.frame(KEGG$Category,KEGG$Term_2,KEGG$PValue,KEGG$Count)
colnames(KEGGPlot)[1:4]<-c("Category","Pathway","pvalue","Count")
KEGGPlot$pvalue<-as.numeric(as.character(KEGGPlot$pvalue))
KEGGPlot$Count<-as.numeric(as.character(KEGGPlot$Count))
KEGGplot2<-KEGGPlot[order(KEGGPlot$pvalue),]
KEGGplot3<-KEGGplot2[1:30,]

KEGGplot3$log10<-log10(KEGGplot3$pvalue)*(-1)

ggplot(data = KEGGplot3, mapping = aes(x = log10, y = Pathway)) +
  geom_point(aes(size=Count,color = -1*log10(pvalue)))+
  scale_color_gradient(low = "blue",high = "red")+
  labs(color=expression(-log[10](Pvalue),size = "Count"))+
  xlab(expression(-log[10](Pvalue)))+facet_grid(Category~.,scales = "free")









########################convert ensembl id to gene symbol#########################
tsv_matrix<- fread("tsv-datatype-0.csv",sep=",")
# tsv_matrix<- tidyr :: separate(tsv_matrix, Gene,into = c('gene_id' , 'junk'), sep='\\.')
# tsv_matrix<- tsv_matrix[,-2]
# awk '{if(!NF || /^#/){next}}1' gencode.v38.annotation.gtf | 
# awk '{$2=null;$6=null;$7=null;$8=null;print $3"\t"$0}'| awk '{if(/^g/){print $0}}'| 
# awk '{print $11"\t"$9"\t"$7"\t"$2"\t"$4"\t"$5}' |sed 's/"//g'| sed 's/;//g' > humanGTF

humanGTF<- read_delim("humanGTF", delim = "\t",
                      escape_double = FALSE, col_names = FALSE,
                      trim_ws = TRUE) %>% dplyr::select(X1,X3)
colnames(humanGTF) <- c("symbol","gene_id")
humanGTF$gene_id <- str_split(humanGTF$gene_id,"[.]",simplify = T)[,1]
humanGTF<-unique(humanGTF)
# humanGTF1<-humanGTF[which(humanGTF$gene_id %in% tsv_matrix$gene_id),]

# grep ">" Homo_sapiens.GRCh37.cdna.all.fa >Homo_sapiens_info.txt
# Homo_sapiens_info.txt----Homo_sapiens_info.xlsx
Homo_sapiens_info<-read_excel("Homo_sapiens_info.xlsx",sheet = 1)

# ensembl_id<-tsv_matrix[,1]
# ensembl_id<- tidyr :: separate(ensembl_id, Gene,into = c('gene_id' , 'junk'), sep='\\.')
# ensembl_id<-ensembl_id[,-2]

Homo_sapiens_info<-tidyr :: separate(Homo_sapiens_info, ENSG,into = c('gene_id' , 'junk'), sep='\\.')
Homo_sapiens_info<-Homo_sapiens_info[,-3]
convert<-merge(Homo_sapiens_info,humanGTF, by.x="gene_id",by.y="gene_id",all.x=TRUE)
convert<-convert[order(convert$order),]
rownames(convert)<-convert$order
convert<-convert[,-3]
write.csv(convert,file = "ensembl_to_symbol.csv")


















