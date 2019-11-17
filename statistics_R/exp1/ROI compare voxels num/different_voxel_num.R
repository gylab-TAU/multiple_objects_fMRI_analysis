rm(list = ls())

library(ggplot2)
library (plyr)
library(dplyr)
library(ez)
library(tidyr)


setwd("./")

# load and arrange all data:
df_all =  read.csv('table_no_intercept.csv')

df_all = filter (df_all, ROI == "FFA_right" | ROI == "FBA_right" | ROI== "Fus_overlap_right")
df_all$ROI = droplevels(df_all$ROI)
df_all$diff = df_all$A_coef - df_all$B_coef
df_all$sum = df_all$A_coef + df_all$B_coef

mean_df = aggregate(x=df_all[,c(7,8,13,14)], by = list(ROI=df_all$ROI,voxels_num=df_all$voxels_num), FUN = mean)
count_df = aggregate(x=df_all[,1], by = list(ROI=df_all$ROI,voxels_num=df_all$voxels_num), FUN = length)
sem_df = aggregate(x=df_all[,c(7,8,13,14)], by = list(ROI=df_all$ROI,voxels_num=df_all$voxels_num), FUN =  function(x) sd(x)/sqrt(length(x)))


color_vec_all = c("FBA_right" = "#e74c3c",
                  "FFA_right" = "#2980b9",
                  "Fus_overlap_right" = "#9932CC")

# plot num of ss

gg = ggplot(data = count_df, aes (x=voxels_num, y=x, color=ROI)) +
  geom_point (size=3)+
  geom_line (size=1) 
  
# gg = gg+geom_vline(xintercept=30, color = "black", linetype = "dashed", size = 1)
gg = gg+ scale_colour_manual(values = color_vec_all, drop = FALSE)  

gg = gg + scale_x_continuous(name = "# voxels",
                            limits = c(0, 125),
                            breaks = seq(0, 120, by = 10)) +
  
  scale_y_continuous(name = "# subjects",
                     limits = c(0, 16),
                     breaks = seq(0, 16, by = 2)) 
  
gg = gg + theme_minimal()  +
    theme(legend.position="none") +
    theme(axis.text = element_text(face = "italic", size = 16),
          axis.title.x = element_text(vjust = -3, size = 16), # move title away from axis
          axis.title.y = element_text(vjust = 3, size = 16))+ # move away for axis
    theme(legend.text = element_text(face = "italic",size=16))
gg

ggsave(gg, file=paste('fig', "subj num ", '.png'), width = 6.5, height = 5)


# plot Wf (A_coef)

gg = ggplot(data = mean_df, aes (x=voxels_num, y=A_coef, color=ROI)) +
  geom_point (size=3)+
  geom_line (size=1) 

# gg = gg+geom_vline(xintercept=30, color = "black", linetype = "dashed", size = 1)


gg = gg+ geom_errorbar(data = mean_df,
                       aes(x=voxels_num, 
                           min=A_coef - sem_df$A_coef,
                           ymax=A_coef + sem_df$A_coef),
                       width=1) 

gg = gg+ scale_colour_manual(values = color_vec_all, drop = FALSE)  

gg = gg + scale_x_continuous(name = "# voxels",
                             limits = c(0, 125),
                             breaks = seq(0, 120, by = 10)) +
  
  scale_y_continuous(name = "Wf",
                     limits = c(0, 1),
                     breaks = seq(-0.5, 1.5, by = 0.25)) 

gg = gg + theme_minimal()  +
  theme(legend.position="none") +
  theme(axis.text = element_text(face = "italic", size = 16),
        axis.title.x = element_text(vjust = -3, size = 16), # move title away from axis
        axis.title.y = element_text(vjust = 3, size = 16))+ # move away for axis
  theme(legend.text = element_text(face = "italic",size=16))
gg

ggsave(gg, file=paste('fig', "Wf ", '.png'), width = 6.5, height = 5)



# plot Wb (B_coef)

gg = ggplot(data = mean_df, aes (x=voxels_num, y=B_coef, color=ROI)) +
  geom_point (size=3)+
  geom_line (size=1) 

# gg = gg+geom_vline(xintercept=30, color = "black", linetype = "dashed", size = 1)


gg = gg+ geom_errorbar(data = mean_df,
                       aes(x=voxels_num, 
                           min=B_coef - sem_df$B_coef,
                           ymax=B_coef + sem_df$B_coef),
                       width=1) 

gg = gg+ scale_colour_manual(values = color_vec_all, drop = FALSE)  

gg = gg + scale_x_continuous(name = "# voxels",
                             limits = c(0, 125),
                             breaks = seq(0, 120, by = 10)) +
  
  scale_y_continuous(name = "Wb",
                     limits = c(0, 1),
                     breaks = seq(-0.5, 1.5, by = 0.25)) 

gg = gg + theme_minimal()  +
  theme(legend.position="none") +
  theme(axis.text = element_text(face = "italic", size = 16),
        axis.title.x = element_text(vjust = -3, size = 16), # move title away from axis
        axis.title.y = element_text(vjust = 3, size = 16))+ # move away for axis
  theme(legend.text = element_text(face = "italic",size=16))
gg

ggsave(gg, file=paste('fig', "Wb ", '.png'), width = 6.5, height = 5)


# plot diff 

gg = ggplot(data = mean_df, aes (x=voxels_num, y=diff, color=ROI)) +
  geom_point (size=3)+
  geom_line (size=1) 

# gg = gg+geom_vline(xintercept=30, color = "black", linetype = "dashed", size = 1)


gg = gg+ geom_errorbar(data = mean_df,
                       aes(x=voxels_num, 
                           min=diff - sem_df$diff,
                           ymax=diff + sem_df$diff),
                       width=1) 

gg = gg+ scale_colour_manual(values = color_vec_all, drop = FALSE)  

gg = gg + scale_x_continuous(name = "# voxels",
                             limits = c(0, 125),
                             breaks = seq(0, 120, by = 10)) +
  
  scale_y_continuous(name = "diff",
                     limits = c(-0.75, 0.75),
                     breaks = seq(-0.75, 0.75, by = 0.25)) 

gg = gg + theme_minimal()  +
  theme(legend.position="none") +
  theme(axis.text = element_text(face = "italic", size = 16),
        axis.title.x = element_text(vjust = -3, size = 16), # move title away from axis
        axis.title.y = element_text(vjust = 3, size = 16))+ # move away for axis
  theme(legend.text = element_text(face = "italic",size=16))
gg

ggsave(gg, file=paste('fig', "diff ", '.png'), width = 6.5, height = 5)



gg = ggplot(data = count_df, aes (x=voxels_num, y=x, color=ROI)) +
  geom_point (size=3)+
  geom_line (size=1) 

# gg = gg+geom_vline(xintercept=30, color = "black", linetype = "dashed", size = 1)
gg = gg+ scale_colour_manual(values = color_vec_all, drop = FALSE)  

gg = gg + scale_x_continuous(name = "# voxels",
                             limits = c(0, 125),
                             breaks = seq(0, 120, by = 10)) +
  
  scale_y_continuous(name = "# subjects",
                     limits = c(0, 16),
                     breaks = seq(0, 16, by = 2)) 

gg = gg + theme_minimal()  +
  theme(legend.position="bottom") +
  theme(axis.text = element_text(face = "italic", size = 16),
        axis.title.x = element_text(vjust = -3, size = 16), # move title away from axis
        axis.title.y = element_text(vjust = 3, size = 16))+ # move away for axis
  theme(legend.text = element_text(face = "italic",size=16))
gg

ggsave(gg, file=paste('fig', "legend ", '.png'), width = 6.5, height = 5)




