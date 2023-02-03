load('01_Normalization.RData')
library(parallel)
library(Rtsne)

####  Unsupervised analysis
####  tSNE plot
tsne = Rtsne(data_train_norm,perplexity = 20,dims = 2);
df_plot = data.frame(sample = rownames(data_train_norm));
df_plot$SNE1 = tsne$Y[,1];
df_plot$SNE2 = tsne$Y[,2];
df_plot$class = group;
p = ggplot(df_plot,aes(x = SNE1,y = SNE2,color = class)) + 
  geom_point(size = 2) +
  labs(color = '') +
  theme_classic() + 
  theme(plot.title = element_text(size = 20,hjust = 0.5,face = 'bold'),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 10),
        legend.position = 'bottom',
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20))
print(p)
ggsave(filename = "tSNE_stju_19class.jpg",
       plot = p,width = 8,height = 10,dpi = 600)

tsne = Rtsne(data_test_norm,perplexity = 20,dims = 2);
df_plot = data.frame(sample = rownames(data_test_norm));
df_plot$SNE1 = tsne$Y[,1];
df_plot$SNE2 = tsne$Y[,2];
df_plot$class = sample_info$class[match(df_plot$sample,sample_info$sample)];
p = ggplot(df_plot[df_plot$sample%in%samp_noncap,],
           aes(x = SNE1,y = SNE2,color = class)) + 
           geom_point(size = 2) +
           labs(color = '') +
           theme_classic() + 
          theme(plot.title = element_text(size = 20,hjust = 0.5,face = 'bold'),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 10),
          legend.position = 'bottom',
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 20))
print(p)
ggsave(filename = "tSNE_righton_9class.jpg",
       plot = p,width = 7,height = 8,dpi = 600)

####  Matrix distance predict
library(class)
pred_train = knn(data_train_norm,data_train_norm,group,
                k = 5,prob = T);
pred_test = knn(data_train_norm,data_test_norm,group,
                k = 5,prob = T);
mean(group==pred_train)
group_test = sample_info$class[match(rownames(data_test_norm),sample_info$sample)];
index_noncap = rownames(data_test_norm)%in%samp_noncap;
mean(group_test==pred_test);
mean(group_test[index_noncap]==pred_test[index_noncap]);
res_out_train_knn = data.frame(sample = rownames(data_train_norm),
                               class = group,pred = pred_train);
write.table(res_out_train_knn,'02_Model_train_knn.txt',
            col.names = T,row.names = T,sep = '\t',quote = F);
res_out_test_knn = data.frame(sample = rownames(data_test_norm),
                               class = group_test,pred = pred_test);
write.table(res_out_test_knn,'02_Model_test_knn.txt',
            col.names = T,row.names = T,sep = '\t',quote = F);

####  Supervised analysis
##  randomForest model training
library(randomForest)
randf = randomForest(x = data_train_norm,y = factor(group),
                     ntree = 1000);
mean(randf$predicted==group)
res_out_train = data.frame(randf$votes);
res_out_train$class = group;
res_out_train$pred_primary = randf$predicted;
res_out_train$pred_second = apply(randf$votes,1,function(x){
  res = randf$classes[order(x,decreasing = T)[2]];
  return(res)
});
mean(res_out_train$class==res_out_train$pred_primary);
mean(res_out_train$class==res_out_train$pred_primary|
       res_out_train$class==res_out_train$pred_second);
write.table(res_out_train,'02_Model_train_cv.txt',
            col.names = T,row.names = T,sep = '\t',quote = F);
res_train_table = table(res_out_train$class,res_out_train$pred_primary);
write.table(res_train_table,'02_Model_train_cv_table.txt',
            col.names = T,row.names = T,sep = '\t',quote = F);

##  Single-blind test
Mat_gene800_database = fread('~/Desktop/RNA_seq_cluster/gene800_database.csv');
sample_info = data.frame(sample = colnames(Mat_gene800_database)[-1],
                         class = as.character(as.matrix(Mat_gene800_database[1,-1])));
sample_info$class = gsub('BCR-ABL1','Ph',sample_info$class);
sample_info$class = gsub('_rearrange','',sample_info$class);
sample_info$class = gsub('PAX5.P80R','PAX5_P80R',sample_info$class);
sample_info$class = gsub('PAX5ALT','PAX5alt',sample_info$class);
sample_info$class = gsub('TCF3-HLF','HLF',sample_info$class);

vote_test_randf = predict(randf,data_test_norm,type = 'vote');
pred_test_randf = predict(randf,data_test_norm);
res_out_test = data.frame(vote_test_randf);
res_out_test$class = sample_info$class[match(rownames(data_test_norm),sample_info$sample)];
res_out_test$pred_primary = apply(vote_test_randf,1,function(x){
  res = randf$classes[order(x,decreasing = T)[1]];
  return(res)
});
res_out_test$pred_second = apply(vote_test_randf,1,function(x){
  res = randf$classes[order(x,decreasing = T)[2]];
  return(res)
});
samp_noncap = sample_info$sample[1:230];
index_noncap = rownames(data_test_norm)%in%samp_noncap;
res_out_test$type = ifelse(rownames(res_out_test)%in%samp_noncap,'Non-captured','Captured');
res_out_test_nc = res_out_test[res_out_test$type=='Non-captured',];
mean(res_out_test_nc$class==res_out_test_nc$pred_primary);
mean(res_out_test_nc$class==res_out_test_nc$pred_primary|
       res_out_test_nc$class==res_out_test_nc$pred_second);

write.table(res_out_test,'02_Model_test.txt',
            col.names = T,row.names = T,sep = '\t',quote = F);
res_test_table = table(res_out_test$class,res_out_test$pred_primary);
write.table(res_test_table,'02_Model_test_table_all.txt',
            col.names = T,row.names = T,sep = '\t',quote = F);
res_test_table_nc = table(res_out_test$class[res_out_test$type=='Non-captured'],
                          res_out_test$pred_primary[res_out_test$type=='Non-captured']);
write.table(res_test_table_nc,'02_Model_test_table_nc.txt',
            col.names = T,row.names = T,sep = '\t',quote = F);
save(data_train_norm,data_test_norm,group,Mat_gene800_database,randf,
     file = '02_Model.RData')
