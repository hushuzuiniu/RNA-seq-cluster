library(ggplot2)
library(pheatmap)
library(data.table)
library(readxl)
#args<-commandArgs(T)
#setwd(args[1])
#setwd("~/Desktop/RNA_seq_cluster/wcy/");
data_train <- fread("~/Desktop/RNA_seq_cluster/wcy/sdju_19_database_final.csv");
data_test = fread('~/Desktop/RNA_seq_cluster/wcy/no_gold_standard_database.csv');
#sample<-read.table("/software/Righton_Software/database/mRNA_Expression_database/Acute_lymphoblastic_leukemia/sjdu.Righton.Expression.samplesheet.xls",header=TRUE)

######## Training
#### Normalization
train_df = data.frame(data_train[-1,-1]);
rownames(train_df) = data_train$Gene[-1];
train_mat = matrix(as.numeric(as.matrix(train_df)),nrow=nrow(train_df));
rownames(train_mat) = data_train$Gene[-1];
##  筛选read count>1的基因
gene_filter = rownames(train_df)[rowMeans(train_mat)>1];
train_mat_fitered = t(train_mat[gene_filter,]); ## 转置
rownames(train_mat_fitered) = colnames(data_train)[-1];
##  筛选参考样本
Total_count = rowSums(train_mat_fitered);
train_mat_total = train_mat_fitered/Total_count;
q3_total = apply(train_mat_total,1,quantile,prob = 0.75);
q3_mean = mean(q3_total);
samp_cand = which.min(abs(q3_total - q3_mean));
train_count_cand = train_mat_total[samp_cand,];
##  筛选代表性基集因
log_fold_diff = apply(train_mat_total,1,function(x){
  res = log2(train_count_cand/x);
  return(res)
})
gene_set1 = apply(log_fold_diff,2,function(x){
  q30 = quantile(x,na.rm = T,prob = 0.3);
  q70 = quantile(x,na.rm = T,prob = 0.7);
  y = ifelse(x>q30 & x<q70,1,0);
  return(y);
})  ##对log_fold_diff进行排序选择其中中间的40%
count_mean = apply(train_mat_total,1,function(x){
  res = (log2(train_count_cand) + log2(x))/2;
  return(res)
})
gene_set2 = apply(count_mean,2,function(x){
  q05 = quantile(x,na.rm = T,prob = 0.05);
  q95 = quantile(x,na.rm = T,prob = 0.95);
  y = ifelse(x>q05 & x<q95,1,0);
  return(y);
})  ##对count_mean进行排序选择其中中间的90%
gene_set_overlap = gene_set1 * gene_set2;  ##取交集

##  计算代表性基因集的log fold的加权平均数
weighted_sum = colSums(t(train_mat_fitered) * 
                         log_fold_diff * gene_set_overlap,na.rm = T)/
  colSums(t(train_mat_fitered) * gene_set_overlap,na.rm=T);
scale_factor = 2^weighted_sum;
scale_factor[is.na(scale_factor)] = 1;
scale_factor_mean = exp(mean(log(scale_factor)));
scale_factor_norm = scale_factor/scale_factor_mean;
#### fix gene_filter,samp_cand,train_count_cand,scale_factor_mean

####  Normalization for test sample
test_df = data.frame(data_test[-1,-1]);
rownames(test_df) = data_test$Gene[-1];
test_mat = matrix(as.numeric(as.matrix(test_df)),nrow=nrow(test_df));
rownames(test_mat) = data_test$Gene[-1];
##  筛选read conut>1的基因
test_mat_fitered = t(test_mat[gene_filter,]); ## 转置
rownames(test_mat_fitered) = colnames(data_test)[-1];
##  筛选参考样本
Total_count = rowSums(test_mat_fitered);
test_mat_total = test_mat_fitered/Total_count;
##  筛选代表性基集因
log_fold_diff_test = apply(test_mat_total,1,function(x){
  res = log2(train_count_cand/x);
  return(res)
})
gene_set1_test = apply(log_fold_diff_test,2,function(x){
  q30 = quantile(x,na.rm = T,prob = 0.3);
  q70 = quantile(x,na.rm = T,prob = 0.7);
  y = ifelse(x>q30 & x<q70,1,0);
  return(y);
})  ##对log_fold_diff进行排序选择其中中间的40%
count_mean_test = apply(test_mat_total,1,function(x){
  res = (log2(train_count_cand) + log2(x))/2;
  return(res)
})
gene_set2_test = apply(count_mean_test,2,function(x){
  q05 = quantile(x,na.rm = T,prob = 0.05);
  q95 = quantile(x,na.rm = T,prob = 0.95);
  y = ifelse(x>q05 & x<q95,1,0);
  return(y);
})  ##对count_mean进行排序选择其中中间的90%
gene_set_overlap_test = gene_set1_test * gene_set2_test;  ##取交集

##  计算代表性基因集的log fold的加权平均数
weighted_sum_test = colSums(t(test_mat_fitered) * 
                         log_fold_diff_test * gene_set_overlap_test,na.rm = T)/
  colSums(t(test_mat_fitered) * gene_set_overlap_test,na.rm=T);
scale_factor_test = 2^weighted_sum_test;
scale_factor_test[is.na(scale_factor_test)] = 1;
scale_factor_norm_test = scale_factor_test/scale_factor_mean;

####
group = as.character(as.matrix(data_train[1,-1]));
data_train_norm = train_mat_total*scale_factor_norm;
data_test_norm = test_mat_total*scale_factor_norm_test;

save(data_train_norm,data_test_norm,group,
     gene_filter,samp_cand,train_count_cand,scale_factor_mean,
     file = '../数据挖掘/Normalization.RData')

#heatmap
#hmExp=log10(y$pseudo.counts+0.001)
hmExp=log10(t(data_train_norm)+0.001)
hmMat=as.matrix(hmExp)
if(FALSE){
  ann_colors =list(clinical=c("Ph-like"         ="#a50026",
                              "ETV6-RUNX1"       ="#d73027", 
                              "Low_hyperdiploid"="#f46d43",
                              "MEF2D"            ="#fdae61", 
                              "TCF3-PBX1"        ="#fee090", 
                              "High_hyperdiploid"="#ffffbf",
                              "Ph"               ="#e0f3f8",
                              "ZNF384"           ="#abd9e9", 
                              "KMT2A"            ="#74add1", 
                              "Low_hypodiploid"  ="#4575b4", 
                              "Target"           ="#a50026"))
}
ann_colors =list(clinical=c("ETV6-RUNX1-like" ="#a50026",
                            "DUX4"             ="#d73027",
                            "ETV6-RUNX1"       ="#f46d43",
                            "High_hyperdiploid"="#fdae61",
                            "HLF"              ="#fee090",
                            "iAMP21"           ="#ffffbf",
                            "IKZF1_N159Y"      ="#e0f3f8",
                            "KMT2A"            ="#abd9e9",
                            "Low_hyperdiploid" ="#74add1",
                            "ZNF384-like"      ="#4575b4",
                            "MEF2D"            ="#313695",
                            "Near_haploid"     ="#8e0152",
                            "NUTM1"            ="#c51b7d",
                            "PAX5_P80R"        ="#de77ae",
                            "PAX5alt"          ="#f1b6da",
                            "Ph"               ="#fde0ef",
                            "Ph-like"          ="#c7eae5",
                            "TCF3-PBX1"        ="#80cdc1",
                            "ZNF384"           ="#35978f",
                            "Target"           ="#01665e"))
annotation_col = data.frame(clinical=factor(group))
rownames(annotation_col) = colnames(hmMat)
pheatmap(hmMat,#filename=paste("sjdu_data_heatmap.pdf",sep="."),
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",show_colnames=FALSE,
         cluster_rows=TRUE,cluster_cols=TRUE,treeheight_row=0,
         treeheight_col=0,annotation_col = annotation_col,
         annotation_colors = ann_colors,fontsize_col=3,
         show_rownames=FALSE,
         main="Acute Lymphoblastic Leukemia mRNA Heatmap",
         scale = "row",color=colorRampPalette(rev(c("red","white","blue")))(102))


