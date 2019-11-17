rm(list = ls())

# library(ggplot2)
library(dplyr)
library(lme4)

setwd("./")


ROI = "EVC"

# load and arrange all data:
df_FB =  read.csv('EVC_FB.csv')
df_FO =  read.csv('EVC_FO.csv')

df_FB = mutate(df_FB, voxel_uid = subj*10^6+voxel_ind)
df_FO = mutate(df_FO, voxel_uid = subj*10^6+voxel_ind)

df_all=inner_join(df_FB, df_FO, by = "voxel_uid")

df_all = filter(df_all,  df_all$sum_FO>0.1  &  df_all$sum_FB>0.1)

# since EVC is much larger than the category-selective cortex, randomly choose spheres from EVC such that it would have
# the same size as category-selective cortex (2023 spheres)
category_selective_length = 2023
set.seed(1)
df_all = df_all[sample(1:nrow(df_all), category_selective_length, replace=FALSE),]


df_all$subj.x = as.factor(df_all$subj.x)

N_subj = nlevels(df_all$subj.x)

subj_list = unique(df_all$subj.x)


cor_mat = matrix(NA,N_subj,12)

for (subj_itr in 1:N_subj){
  
  curr_data= filter(df_all, subj.x==subj_list[subj_itr])
  
  cor_mat[subj_itr,1] = cor(curr_data$Wf_FB, curr_data$face_others_t_loc.x) 
  cor_mat[subj_itr,2] = cor(curr_data$Wf_FB, curr_data$body_others_t_loc.x)
  cor_mat[subj_itr,3] = cor(curr_data$Wf_FB, curr_data$object_others_t_loc.x)
  cor_mat[subj_itr,4] = cor(curr_data$Wb_FB, curr_data$face_others_t_loc.x) 
  cor_mat[subj_itr,5] = cor(curr_data$Wb_FB, curr_data$body_others_t_loc.x)
  cor_mat[subj_itr,6] = cor(curr_data$Wb_FB, curr_data$object_others_t_loc.x)
  
  cor_mat[subj_itr,7] = cor(curr_data$Wf_FO, curr_data$face_others_t_loc.x) 
  cor_mat[subj_itr,8] = cor(curr_data$Wf_FO, curr_data$body_others_t_loc.x)
  cor_mat[subj_itr,9] = cor(curr_data$Wf_FO, curr_data$object_others_t_loc.x)
  cor_mat[subj_itr,10] = cor(curr_data$Wo_FO, curr_data$face_others_t_loc.x) 
  cor_mat[subj_itr,11] = cor(curr_data$Wo_FO, curr_data$body_others_t_loc.x)
  cor_mat[subj_itr,12] = cor(curr_data$Wo_FO, curr_data$object_others_t_loc.x)
  
}

z_cor_mat = atanh(cor_mat)

# mean fisher z for correlations
z_mean = colMeans(z_cor_mat)
r_mean = tanh(z_mean)

p_values = vector()

for (cor_itr in 1:12){
  
  t=t.test(z_cor_mat[,cor_itr])
  p_values[cor_itr] = t$p.value
}

# Bonferroni correction for multiple comparisons (24 beacause it includes same set of comparisons for EVC)
p_values = p_values*24

p_vals_for_plot = matrix(data = c(p_values), nrow = 4,ncol =3,byrow = TRUE)
rownames(p_vals_for_plot) = c("Wf","Wb","Wf","Wo")
colnames(p_vals_for_plot) = c("face_t","body_t","object_t")


z_mean_for_plot = matrix(data = c(z_mean), nrow = 4,ncol =3,byrow = TRUE)
rownames(z_mean_for_plot) = c("Wf","Wb","Wf","Wo")
colnames(z_mean_for_plot) = c("face_t","body_t","object_t")

## plot figure 8b:

png(filename = paste("both",ROI,  "corr_weights_selectivity_ellipse", '.png'), height = 500, width = 400)
corrplot(z_mean_for_plot,method ="ellipse",
         tl.pos='n',
         cl.length = 5,
         cl.cex=1.5,
         cl.ratio = 0.2,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation)
         number.cex = 2,
         cl.align.text = 'l',
         cl.offset = 0.5)

dev.off()

df_fig8 = select(df_all, subj.x, Wf_FB, Wb_FB, 
                 Wf_FO, Wo_FO, 
                 face_others_t_loc.x, body_others_t_loc.x,object_others_t_loc.x)
write.csv(df_fig8, "Fig8_data_2_evc.csv")

