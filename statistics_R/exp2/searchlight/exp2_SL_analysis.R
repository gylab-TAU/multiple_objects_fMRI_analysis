rm(list = ls())

# library(ggplot2)
library(dplyr)
library(lme4)


setwd("./")

ROI = "category_right"

# load and arrange all data:

df_FB =  read.csv('Category_Right_FB.csv')
df_FO =  read.csv('Category_Right_FO.csv')


# combine the two datasets to the same df
df_FB = mutate(df_FB, voxel_uid = subj*10^6+voxel_ind)
df_FO = mutate(df_FO, voxel_uid = subj*10^6+voxel_ind)

df_all=inner_join(df_FB, df_FO, by = "voxel_uid")

df_all$subj.x = as.factor(df_all$subj.x)
N_subj = nlevels(df_all$subj.x)
subj_list = unique(df_all$subj.x)

df_fig6 = select(df_all, subj.x, Wf_FB, Wb_FB, sum_FB, diff_FB, Rsq_FB, 
                 Wf_FO, Wo_FO, sum_FO, diff_FO, Rsq_FO,
                 face_body_t_loc.x, face_object_t_loc.x,voxel_num.x)
write.csv(df_fig6, "Fig6_data.csv")


### calculate correlation between difference between beta coefficients and relative selectivity for both pair types

cor_mat = matrix(NA,N_subj,2)

for (subj_itr in 1:N_subj){
  
  curr_data= filter(df_all, subj.x==subj_list[subj_itr])
  # curr_data = select(df_FO, -subj)
  

  cor_mat[subj_itr,1] = cor(curr_data$diff_FO, curr_data$face_object_t_loc.x)
  cor_mat[subj_itr,2] = cor(curr_data$diff_FB, curr_data$face_body_t_loc.x)

}

z_cor_mat = atanh(cor_mat)

# mean fisher z for correlations
z_mean = colMeans(z_cor_mat)
r_mean = tanh(z_mean)

## t test for fisher z correlations:
# FO:
t.test(z_cor_mat[,1])

# FB:
t.test(z_cor_mat[,2])


## calculate mean sum
# mean sum of coefficients

mean_sum_FB = tapply(df_all$sum_FB, df_all$subj.x, mean)
t.test(mean_sum_FB, mu=1)

mean_sum_FO = tapply(df_all$sum_FO, df_all$subj.x, mean)
t.test(mean_sum_FO, mu=1)



#####
## plot figure 6a

data_line = data.frame(Wf_l=c(-0.72, 1.72), Wb_l=c(1.72, -0.72))
data_line2 = data.frame(Wf_l=c(-0.72, 1.72), Wb_l=c(-0.72, 1.72))
data_x_axis = data.frame(x_x_axis= c(-0.72,1.72),y_x_axis= c(0,0))
data_y_axis = data.frame(x_y_axis= c(0,0),y_y_axis=c(-0.72,1.72)) 


gg = ggplot(data = df_all) +
  geom_point(aes(x=Wf_FB, y=Wb_FB, color = face_body_t_loc.x))

gg = gg+scale_colour_gradient2(low = "red", mid = "grey80",
                               high = "blue", midpoint = 0,
                               limits = c(-15,15),
                               name=" ") 

gg = gg+
  geom_line(data=data_line, aes(x=Wf_l, y=Wb_l), size = 1, linetype="dashed") +
  geom_line(data=data_line2, aes(x=Wf_l, y=Wb_l), size = 1, linetype="dashed") +
  
  geom_line(data = data_x_axis, aes(x = x_x_axis, y = y_x_axis)) +
  geom_line(data = data_y_axis, aes(x = x_y_axis, y = y_y_axis)) +
  
  scale_x_continuous(name = "Wf",
                     limits = c(-0.72, 1.72),
                     breaks = seq(-0.5, 1.5, by = 0.5)) +
  scale_y_continuous(name = "Wb",
                     limits = c(-0.72, 1.72),
                     breaks = seq(-0.5, 1.5, by = 0.5)) +
  
  theme_minimal()  +
  theme(legend.position="right") +
  theme(legend.text = element_text(face = "italic",size=10),
        legend.key.size = unit(1, "cm"),
        legend.key.width = unit(0.5,"cm") )+
  
  theme(axis.text = element_text(face = "italic", size = 20),
        axis.title.x = element_text(vjust = -3, size = 20), # move title away from axis
        axis.title.y = element_text(vjust = 3, size = 20)) # move away for axis

gg  = gg+ theme(legend.position="bottom",
                legend.key.size = unit(0.4, "cm"),
                legend.key.width = unit(1,"cm"))

print(gg)

ggsave(gg, file=paste("FB", ROI, "Wb over Wf color by selectivity paper",'.png'), width = 4.5, height = 5)



#####
## plot figure 6b

data_line = data.frame(Wf_l=c(-0.72, 1.72), Wb_l=c(1.72, -0.72))
data_line2 = data.frame(Wf_l=c(-0.72, 1.72), Wb_l=c(-0.72, 1.72))
data_x_axis = data.frame(x_x_axis= c(-0.72,1.72),y_x_axis= c(0,0))
data_y_axis = data.frame(x_y_axis= c(0,0),y_y_axis=c(-0.72,1.72)) 


gg = ggplot(data = df_all) +
  geom_point(aes(x=Wf_FO, y=Wo_FO, color = face_object_t_loc.x))

gg = gg+scale_colour_gradient2(low = "green2", mid = "grey80",
                               high = "blue", midpoint = 0,
                               limits = c(-15,15),
                               name=" ") 

gg = gg+
  geom_line(data=data_line, aes(x=Wf_l, y=Wb_l), size = 1, linetype="dashed") +
  geom_line(data=data_line2, aes(x=Wf_l, y=Wb_l), size = 1, linetype="dashed") +
  
  geom_line(data = data_x_axis, aes(x = x_x_axis, y = y_x_axis)) +
  geom_line(data = data_y_axis, aes(x = x_y_axis, y = y_y_axis)) +
  
  scale_x_continuous(name = "Wf",
                     limits = c(-0.72, 1.72),
                     breaks = seq(-0.5, 1.5, by = 0.5)) +
  scale_y_continuous(name = "Wo",
                     limits = c(-0.72, 1.72),
                     breaks = seq(-0.5, 1.5, by = 0.5)) +
  
  theme_minimal()  +
  theme(legend.position="right") +
  theme(legend.text = element_text(face = "italic",size=10),
        legend.key.size = unit(1, "cm"),
        legend.key.width = unit(0.5,"cm") )+
  
  theme(axis.text = element_text(face = "italic", size = 20),
        axis.title.x = element_text(vjust = -3, size = 20), # move title away from axis
        axis.title.y = element_text(vjust = 3, size = 20)) # move away for axis

gg  = gg+ theme(legend.position="bottom",
                legend.key.size = unit(0.4, "cm"),
                legend.key.width = unit(1,"cm"))

print(gg)

ggsave(gg, file=paste("FO", ROI, "Wo over Wf color by selectivity paper",'.png'), width = 4.5, height = 5)



#################


### calculate correlations between single coefficients and category selectivity


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


## plot figure 8a:

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

M_empty = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0),4,3)


png(filename = paste("both", 'corr_weights_selectivity_ellipse', 'empty.png'), height = 500, width = 400)
corrplot(M_empty,method ="ellipse" ,
         tl.pos='n',
         tl.col = "black", tl.srt = 90,
         cl.length = 5,
         cl.cex=1.5,
         cl.ratio = 0.2,
         cl.align.text = 'l',
         cl.offset = 0.5) # Text label color and rotation

dev.off()

df_fig8 = select(df_all, subj.x, Wf_FB, Wb_FB, 
                 Wf_FO, Wo_FO, 
                 face_others_t_loc.x, body_others_t_loc.x,object_others_t_loc.x)
write.csv(df_fig8, "Fig8_data_1_cat.csv")

