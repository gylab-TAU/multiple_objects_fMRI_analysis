rm(list = ls())

library(ggplot2)
library(dplyr)
library(ez)

setwd("./")


#### load and arrange all data:
df_all =  read.csv('study1_roi_data_30_voxels.csv')
df_all$subj = as.factor(df_all$subj)


attach(df_all)
df_all$region[ROI=="FFA" | ROI=="FBA" | ROI=="Overlap_ventral" ] <- "Ventral"
df_all$region[ROI=="OFA" | ROI=="EBA" | ROI=="Overlap_lat_occ"] <- "Lateral"
detach(df_all)

attach(df_all)
df_all$region_type[ROI=="FFA" | ROI=="OFA" ] <- "Face"
df_all$region_type[ROI=="FBA" | ROI=="EBA"] <- "Body"
df_all$region_type[ROI=="Overlap_ventral" | ROI=="Overlap_lat_occ"] <- "Overlap"
detach(df_all)


####  Calculate beta coeficients for each subject and each ROI:
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

betas_df$diff = betas_df$beta_face - betas_df$beta_body
betas_df$sum = betas_df$beta_face + betas_df$beta_body


#### calculate mean and s.e.m across subjects of beta coefficients
betas_mean = aggregate(x = betas_df[,4:6], by = list(ROI=betas_df$ROI,hemi=betas_df$hemi), FUN = mean)
betas_sem = aggregate(x = betas_df[,4:6], by = list(ROI=betas_df$ROI,hemi=betas_df$hemi), FUN = function(x) sd(x)/sqrt(length(x)))
betas_Ns = aggregate(x = betas_df[,4:6], by = list(ROI=betas_df$ROI,hemi=betas_df$hemi), FUN = length)



ROI_list = c("FFA", "FBA", "Overlap_ventral")

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


#### calculate mean and s.e.m across subjects of beta coefficients
betas_mean = aggregate(x = betas_df[,4:6], by = list(ROI=betas_df$ROI,hemi=betas_df$hemi), FUN = mean)
betas_sem = aggregate(x = betas_df[,4:6], by = list(ROI=betas_df$ROI,hemi=betas_df$hemi), FUN = function(x) sd(x)/sqrt(length(x)))
betas_Ns = aggregate(x = betas_df[,4:6], by = list(ROI=betas_df$ROI,hemi=betas_df$hemi), FUN = length)



### compare response to single and multiple objects and check how many voxels have 
### lower response to multiple objects than one of the single objects

df_all$max_face_body = pmax(df_all$face_PSC, df_all$body_PSC) 

df_all$max_diff = df_all$max_face_body- df_all$person_PSC

curr_df = filter(df_all, hemi=="Right", (ROI=="Overlap_ventral" | ROI=="FFA" | ROI=="FBA"))
curr_df$ROI = droplevels(curr_df$ROI)

curr_df$ROI = relevel(curr_df$ROI, "FFA")

median(curr_df$max_diff)
which(sort(curr_df$max_diff)>0)[[1]]


median((df_all$max_face_body- df_all$person_PSC)/df_all$max_face_body)




#### define functions for plotting beta coefficients

## plot means - a function that plots only the mean and s.e.m of given data 
plot_means = function(curr_betas_mean,  curr_betas_sem, color_vec_all, fig_title){
  
  #data_points = data.frame(Wf_p=c(1,0, 0.5,1), Wb_p= c(0,1,0.5,1))
  data_points = data.frame(Wf_p=c(1), Wb_p= c(1))
  data_line = data.frame(Wf_l=c(-0.5, 1.5), Wb_l=c(1.5, -0.5))
  data_line2 = data.frame(Wf_l=c(-0.5, 1.5), Wb_l=c(-0.5, 1.5))
  data_x_axis = data.frame(x_x_axis= c(-0.5,1.5),y_x_axis= c(0,0))
  data_y_axis = data.frame(x_y_axis= c(0,0),y_y_axis=c(-0.5,1.5)) 
  
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
    # geom_point(data=data_points, aes(x=Wf_p, y=Wb_p), size=5, alpha=1)  +
    
    geom_line(data=data_line, aes(x=Wf_l, y=Wb_l), size = 1, linetype="dashed") +
    geom_line(data = data_line2, aes(x=Wf_l, y=Wb_l), size = 1, linetype="dashed") +
    geom_line(data = data_x_axis, aes(x = x_x_axis, y = y_x_axis)) +
    geom_line(data = data_y_axis, aes(x = x_y_axis, y = y_y_axis)) +
    
    scale_x_continuous(name = "Wf",
                       limits = c(-0.5, 1.5),
                       breaks = seq(0, 1, by = 0.5)) +
    scale_y_continuous(name = "Wb",
                       limits = c(-0.5, 1.5),
                       breaks = seq(0, 1, by = 0.5)) +
    
    theme_minimal()  +
    theme(legend.position="right") +
    theme(axis.text = element_text(face = "italic", size = 20),
          axis.title.x = element_text(vjust = -3, size = 20), # move title away from axis
          axis.title.y = element_text(vjust = 3, size = 20))+ # move away for axis
    theme(legend.text = element_text(face = "italic",size=8))
  
  print(gg)
  
  ggsave(gg, file=paste('fig', "means ",fig_title ,'2.png'), width = 6.5, height = 5)
  
  return (gg)
} 


## add subjects - a function that adds a layer of single subject data over the mean plot from the previous figure
add_subjects = function(curr_betas_df, gg, fig_title){
  
  gg = gg+ geom_point(data = curr_betas_df,
                      aes(x = beta_face, y = beta_body,color = ROI),
                      size = 2,
                      shape = 16)
  
  
  print(gg)
  
  ggsave(gg, file=paste('fig', "means with subjects  ",fig_title ,'2.png'), width = 6.5, height = 5)
  
  return (gg)
}



#### filter data to according to ROIs and hemispheres and draw figures

## Right ventral ROIs:

color_vec_all = c("FBA" = "#e74c3c",
                  "FFA" = "#2980b9",
                  "Overlap_ventral" = "#9932CC")



curr_betas_mean = filter(betas_mean, ((ROI=="FFA" | ROI=="FBA" |ROI=="Overlap_ventral") & hemi=="Right"))
curr_betas_sem = filter(betas_sem, ((ROI=="FFA" | ROI=="FBA" |ROI=="Overlap_ventral")& hemi=="Right"))
curr_betas_df = filter(betas_df, ((ROI=="FFA" | ROI=="FBA" |ROI=="Overlap_ventral")& hemi=="Right"))

curr_betas_mean = droplevels(curr_betas_mean)
curr_betas_sem = droplevels(curr_betas_sem)
curr_betas_df = droplevels(curr_betas_df)

write.csv(curr_betas_df,'Fig3_data.csv')


gg = plot_means(curr_betas_mean, curr_betas_sem, color_vec_all, "Right Ventral lm")
add_subjects(curr_betas_df, gg, "Right Ventral lm")


#######
# Right lateral ROIs

curr_betas_mean = filter(betas_mean, ((ROI=="OFA" | ROI=="EBA") & hemi=="Right"))
curr_betas_sem = filter(betas_sem, ((ROI=="OFA" | ROI=="EBA") & hemi=="Right"))
curr_betas_df = filter(betas_df, ((ROI=="OFA" | ROI=="EBA") & hemi=="Right"))

curr_betas_mean = droplevels(curr_betas_mean)
curr_betas_sem = droplevels(curr_betas_sem)
curr_betas_df = droplevels(curr_betas_df)


color_vec_all = c("EBA" = "#e74c3c",
                  "OFA" = "#2980b9")

rm(gg)
gg = plot_means(curr_betas_mean, curr_betas_sem, color_vec_all, "Right lateral lm")
add_subjects(curr_betas_df, gg, "Right lateral lm")




