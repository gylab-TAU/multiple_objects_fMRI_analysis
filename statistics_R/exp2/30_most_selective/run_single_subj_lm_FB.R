rm(list = ls())

library(ggplot2)
library(dplyr)
library(ez)

setwd("./")

# load and arrange all data:
df_all =  read.csv('study2_ROI_data_30_voxels_FB.csv')
df_all$subj = as.factor(df_all$subj)


#### calculate beta coeficients for each subject and each ROI

ROI_list = levels(df_all$ROI)
hemi_list = levels(df_all$hemi)
curr_hemi = "Right"

betas_df = data.frame(ROI = factor(,levels(df_all$ROI)),
                      hemi = factor(,levels(df_all$hemi)),
                      subj = factor(,levels(df_all$subj)),
                      beta_face = double(),
                      beta_body = double())


for (hemi_itr in seq_along(hemi_list)){
  curr_hemi = hemi_list[hemi_itr]
  
  for (roi_itr in seq_along(ROI_list)) {
    
    # filter data to curr ROI
    curr_ROI = ROI_list[roi_itr]
    curr_data = filter(df_all, ROI==curr_ROI, hemi==curr_hemi)
    curr_data$subj = droplevels(curr_data$subj)
    
    subj_list = unique(curr_data$subj)
    curr_N_subj = length(subj_list)
    
    betas = matrix(NA,curr_N_subj, 2)
    
    for (subj_itr in 1:curr_N_subj){
      
      curr_data_subj = filter(curr_data, subj==subj_list[subj_itr])
      
      fit = lm(curr_data_subj, formula = "person_PSC ~ 0+face_PSC + body_PSC")
      betas[subj_itr,] = fit$coefficients
      
      betas_df[nrow(betas_df) + 1,] = list(curr_ROI,curr_hemi,subj_list[subj_itr],fit$coefficients[1],fit$coefficients[2])
      
    }
  }
}


# calculate diff and sum of betas
betas_df$diff = betas_df$beta_face - betas_df$beta_body
betas_df$sum = betas_df$beta_face + betas_df$beta_body


# calculate summary statistics over the betas
betas_mean = aggregate(x = betas_df[,4:6], by = list(ROI=betas_df$ROI,hemi=betas_df$hemi), FUN = mean)
betas_sem = aggregate(x = betas_df[,4:6], by = list(ROI=betas_df$ROI,hemi=betas_df$hemi), FUN = function(x) sd(x)/sqrt(length(x)))
betas_Ns = aggregate(x = betas_df[,4:6], by = list(ROI=betas_df$ROI,hemi=betas_df$hemi), FUN = length)


# calculate statistics

ROI_list = c("FFA", "FBA")

for (roi_itr in seq_along(ROI_list)) {
  # perform t-test for each ROI: compare the betas to 0 or 1, compare diff to 0 , compare sum to 2
  curr_betas = filter(betas_df, hemi=="Right", ROI==ROI_list[roi_itr])
  print(ROI_list[roi_itr])
  print(t.test (curr_betas$beta_face, mu=0))
  print(t.test (curr_betas$beta_body, mu=0))
  print(t.test (curr_betas$beta_face, mu=1))
  print(t.test (curr_betas$beta_body, mu=1))
  print(t.test (curr_betas$diff, mu=0))
  print(t.test(curr_betas$sum, mu=2))
}




##### Plot Fig. 5a  

# define function - plot froup man and sem
plot_means = function(curr_betas_mean,  curr_betas_sem, color_vec_all, fig_title){
  
  data_points = data.frame(Wf_p=c(1,0, 0.5,1), Wb_p= c(0,1,0.5,1))
  data_line = data.frame(Wf_l=c(-0.25, 1.25), Wb_l=c(1.25, -0.25))
  data_line2 = data.frame(Wf_l=c(-0.25, 1.25), Wb_l=c(-0.25, 1.25))
  
  data_x_axis = data.frame(x_x_axis= c(-0.25,1.25),y_x_axis= c(0,0))
  data_y_axis = data.frame(x_y_axis= c(0,0),y_y_axis=c(-0.25,1.25)) 
  
  gg = ggplot()
  gg = gg + geom_point(data = curr_betas_mean,
                       aes(x = beta_face, y = beta_body, color = ROI),
                       size = 8,
                       shape = 18) +
    
    scale_colour_manual(values = color_vec_all, drop = FALSE)  +
    #       
    geom_errorbar(data = curr_betas_mean,
                  aes(x=curr_betas_mean$beta_face, 
                      min=beta_body - curr_betas_sem$beta_body,
                      ymax=beta_body + curr_betas_sem$beta_body),
                  width=.05) +
    
    geom_errorbarh(data = curr_betas_mean,
                   aes(y= curr_betas_mean$beta_body, 
                       xmin=beta_face - curr_betas_sem$beta_face,
                       xmax=beta_face + curr_betas_sem$beta_face), 
                   height=.05) 
  
  
  gg =  gg+ 
    
    #geom_point(data=data_points, aes(x=Wf_p, y=Wb_p), size=5, alpha=1)  +
    
    geom_line(data=data_line, aes(x=Wf_l, y=Wb_l), size = 1, linetype="dashed") +
    geom_line(data = data_line2, aes(x=Wf_l, y=Wb_l), size = 1, linetype="dashed") +
    geom_line(data = data_x_axis, aes(x = x_x_axis, y = y_x_axis)) +
    geom_line(data = data_y_axis, aes(x = x_y_axis, y = y_y_axis)) +
    
    scale_x_continuous(name = "Wf",
                       limits = c(-0.25, 1.25),
                       breaks = seq(0, 1, by = 0.5)) +
    scale_y_continuous(name = "Wb",
                       limits = c(-0.25, 1.25),
                       breaks = seq(0, 1, by = 0.5)) +
    
    theme_minimal()  +
    theme(legend.position="none") +
    theme(axis.text = element_text(face = "italic", size = 20),
          axis.title.x = element_text(vjust = -3, size = 16), # move title away from axis
          axis.title.y = element_text(vjust = 3, size = 16)) # move away for axis
  
  print(gg)
  
  ggsave(gg, file=paste('fig', "means ",fig_title ,'.png'), width = 6, height = 5.5)
  
  return (gg)
} 


# define function - add single subject data to mean figure
add_subjects = function(curr_betas_df, gg, fig_title){
  
  gg = gg+ geom_point(data = curr_betas_df,
                      aes(x = beta_face, y = beta_body,color = ROI),
                      size = 2,
                      shape = 16)
  
  
  print(gg)
  
  #ggsave(gg, file=paste('fig', "means with subjects ",fig_title ,'.png'), width = 6, height = 5.5)
  ggsave(gg, file=paste('fig', "means with subjects ",fig_title ,'.png'), width = 4.5, height = 4)
  
  
  return (gg)
}




# filter data and plot

color_vec_all = c("FBA" = "#e74c3c",
                  "FFA" = "#2980b9")


curr_betas_mean = filter(betas_mean, ((ROI=="FFA" | ROI=="FBA" ) & hemi=="Right"))
curr_betas_sem = filter(betas_sem, ((ROI=="FFA" | ROI=="FBA" )& hemi=="Right"))
curr_betas_df = filter(betas_df, ((ROI=="FFA" | ROI=="FBA" )& hemi=="Right"))

curr_betas_mean = droplevels(curr_betas_mean)
curr_betas_sem = droplevels(curr_betas_sem)
curr_betas_df = droplevels(curr_betas_df)

gg = plot_means(curr_betas_mean, curr_betas_sem, color_vec_all, "FB_Right Ventral FFA FBA lm")
add_subjects(curr_betas_df, gg, "FB_Right Ventral FFA FBA lm")

write.csv(curr_betas_df, "Fig5_data_1_FB.csv")



