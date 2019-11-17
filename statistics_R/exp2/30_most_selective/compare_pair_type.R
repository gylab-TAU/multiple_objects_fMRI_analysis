rm(list = ls())

library(ggplot2)
library(dplyr)
library(ez)
library(lme4)


setwd("./")


## FO = Face-Object data
## FB = Face-Body data

# load and arrange all data:
# FO data:
df_all =  read.csv('study2_ROI_data_30_voxels_FO.csv')
df_all$subj = as.factor(df_all$subj)

ROI_list = levels(df_all$ROI)
hemi_list = levels(df_all$hemi)
curr_hemi = "Right"

betas_df = data.frame(design = factor(,list("FB","FO")),
                      ROI = factor(,levels(df_all$ROI)),
                      hemi = factor(,levels(df_all$hemi)),
                      subj = factor(,levels(df_all$subj)),
                      beta_face = double(),
                      beta_other = double())


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
      
      fit = lm(curr_data_subj, formula = "integratedFO_PSC ~ 0+face_PSC + object_PSC")
      betas[subj_itr,] = fit$coefficients
      
      betas_df[nrow(betas_df) + 1,] = list("FO", curr_ROI,curr_hemi,subj_list[subj_itr],fit$coefficients[1],fit$coefficients[2])
      
    }
  }
}


# FB data:
df_all =  read.csv('study2_ROI_data_30_voxels_FB.csv')
df_all$subj = as.factor(df_all$subj)

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
      
      betas_df[nrow(betas_df) + 1,] = list("FB", curr_ROI,curr_hemi,subj_list[subj_itr],fit$coefficients[1],fit$coefficients[2])
      
    }
  }
}



# calculate statistics:

betas_df$diff = betas_df$beta_face - betas_df$beta_other
betas_df$sum = betas_df$beta_face + betas_df$beta_other

curr_betas = filter(betas_df, 
                    hemi=="Right",
                    ((ROI=="FFA" | ROI=="FBA") & design=="FB") |( (ROI=="FFA" | ROI=="medial_obj")& design=="FO"))

curr_betas$new_ROI[curr_betas$ROI=="FFA"]= "upper"
curr_betas$new_ROI[curr_betas$ROI!="FFA"]= "lower"

curr_betas$new_ROI = as.factor(curr_betas$new_ROI)

# remove subjects with incomplete data (missing ROIs) from analysis
curr_betas2 = filter(curr_betas, !(subj=="3" | subj=="6"| subj=="12"))
curr_betas2$subj = droplevels(curr_betas2$subj)

anova_results <- ezANOVA(curr_betas2, # specify data frame
                     dv = diff, # specify dependent variable 
                     wid = subj, # specify the subject variable
                     within = .(new_ROI, design), # specify within-subject variables
                     detailed = TRUE # get a detailed table that includes SS
)

print (anova_results)
