rm(list = ls())

# library(ggplot2)
library(dplyr)
library(lme4)

setwd("./")

# load and arrange all data:
df_cat =  read.csv('Category_Right_exclude_neighbors.csv')
df_cat$subj = as.factor(df_cat$subj)

N_subj = nlevels(df_cat$subj)

# save source data for Figure 5
df_fig4 = select(df_cat, subj, Wf,Wb,sum,diff,Rsq,face_body_s,voxel_num)
write.csv(df_fig4, "Fig4_data.csv")

### calaculate statistics

# mean sum of coefficients
mean_sum = tapply(df_cat$sum, df_cat$subj, mean)
t.test(mean_sum, mu=1)


subj_list = unique(df_cat$subj)

#### Correlation between difference between beta coefficients (curr_data$diff) and selctivity to faces relative to 
# bodies (t value for Face>Body from localizer data curr_data$face_body_s)
cor_mat = vector(length = N_subj)
  
for (subj_itr in 1:N_subj){
 
  curr_data= filter(df_cat, subj==subj_list[subj_itr])
  
  cor_mat[subj_itr] = cor(curr_data$diff, curr_data$face_body_s)

}
    

z_cor_mat = atanh(cor_mat)
t.test(z_cor_mat)

Z_mean = mean(z_cor_mat)
r_mean_value = tanh(z_mean)


# median R-squared
mean(df_cat$Rsq)
median(df_cat$Rsq)
median(df_evc$Rsq)




##### plot all data - color by selectivity (Fig. 4a)

data_line = data.frame(Wf_l=c(-0.7, 1.7), Wb_l=c(1.7, -0.7))
data_line2 = data.frame(Wf_l=c(-0.7, 1.7), Wb_l=c(-0.7, 1.7))
data_x_axis = data.frame(x_x_axis= c(-0.7,1.7),y_x_axis= c(0,0))
data_y_axis = data.frame(x_y_axis= c(0,0),y_y_axis=c(-0.7,1.7)) 


gg = ggplot(data = df_cat) +
  geom_point(aes(x=Wf, y=Wb, color = face_body_s))

gg = gg+scale_colour_gradient2(low = "red", mid = "grey80",
                               high = "blue", midpoint = 0,
                               limits = c(-15,15),
                               name=" ")
gg =  gg+ 
  

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
print(gg)

ggsave(gg, file=paste("FB_Right", "Wb over Wf color by selectivity",'.png'), width = 6, height = 5)



#### plot R-squared figure (Fig. 4b)

gg = ggplot(data=df_cat, aes(Rsq)) + 
  #geom_freqpoly (bins=50) +
  geom_histogram(binwidth = 0.01 ) +
  theme_minimal()+
  theme(axis.text = element_text(face = "italic", size = 14),
        axis.title.x = element_text(vjust = -3, size = 14), # move title away from axis
        axis.title.y = element_text(vjust = 3, size = 14)) # move away for axis

gg = gg+ scale_x_continuous(name = "Rsq",
                            limits = c(-0.5, 1),
                            breaks = seq(-0.5, 1, by = 0.5)) 

gg  = gg+ theme(legend.position="none",
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                axis.text = element_text(face = "italic", size = 20))

ggsave(gg, file=paste('fig', "Rsq hist",'.png'), width = 4, height = 5)


