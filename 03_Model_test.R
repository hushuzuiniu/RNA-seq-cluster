load('01_Normalization.RData')

data_test = read.csv('../分析数据及代码——胡书/RA_add_database.csv');

test_df = data.frame(data_test[-1,-1]);
rownames(test_df) = data_test$Gene[-1];
test_mat = matrix(as.numeric(as.matrix(test_df)),nrow=nrow(test_df));
rownames(test_mat) = data_test$Gene[-1];
##  筛选read count>1的基因
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
group_test = as.character(as.matrix(data_test[1,-1]));
data_test_norm_add = test_mat_total*scale_factor_norm_test;

load('02_Model.RData')
vote_test_randf = predict(randf,data_test_norm_add,type = 'vote');
pred_test_randf = predict(randf,data_test_norm_add);
res_out_test = data.frame(vote_test_randf);
res_out_test$class = group_test
res_out_test$pred_primary = apply(vote_test_randf,1,function(x){
  res = randf$classes[order(x,decreasing = T)[1]];
  return(res)
});
res_out_test$pred_second = apply(vote_test_randf,1,function(x){
  res = randf$classes[order(x,decreasing = T)[2]];
  return(res)
});
res_out_test$type = 'Non-captured';
res_out_test_nc = res_out_test[res_out_test$type=='Non-captured',];
mean(res_out_test_nc$class==res_out_test_nc$pred_primary);
mean(res_out_test_nc$class==res_out_test_nc$pred_primary|
       res_out_test_nc$class==res_out_test_nc$pred_second);
write.table(res_out_test,'02_Model_test_add.txt',
            col.names = T,row.names = T,sep = '\t',quote = F);

data_test_all = rbind(data_test_norm,data_test_norm_add);
tsne = Rtsne(data_test_all,perplexity = 20,dims = 2);
df_plot = data.frame(sample = rownames(data_test_all));
df_plot$SNE1 = tsne$Y[,1];
df_plot$SNE2 = tsne$Y[,2];
df_plot$class = c(sample_info$class[match(rownames(data_test_norm),sample_info$sample)],group_test);
p = ggplot(df_plot[df_plot$sample%in%samp_noncap|
                     df_plot$sample%in%rownames(res_out_test_nc),],
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
ggsave(filename = "tSNE_righton_12class.jpg",
       plot = p,width = 7,height = 8,dpi = 600)


