rm(list = ls())

# library(ggplot2)
library(dplyr)
# library(ez)
library(rstan) # observe startup messages
library(bridgesampling)

# show warnings when they occur:
options(warn=1)

# for stan:
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

setwd("./")

# load and arrange all data:
df_all =  read.csv('study1_roi_data_30_voxels.csv')
# df_right = filter(df_all, hemi=="Right")


attach(df_all)
df_all$region[ROI=="FFA" | ROI=="FBA" | ROI=="Overlap_ventral" ] <- "Ventral"
df_all$region[ROI=="OFA" | ROI=="EBA" | ROI=="Overlap_lat_occ"] <- "Lateral"
detach(df_all)

attach(df_all)
df_all$region_type[ROI=="FFA" | ROI=="OFA" ] <- "Face"
df_all$region_type[ROI=="FBA" | ROI=="EBA"] <- "Body"
df_all$region_type[ROI=="Overlap_ventral" | ROI=="Overlap_lat_occ"] <- "Overlap"
detach(df_all)

bf_stats = data.frame(bf_fit_over_intercept = double(),
                      bf_fit_over_interaction = double(),
                      ROI =  character())


# arrange data for a single ROI:
 # curr_ROI="Overlap_ventral"
curr_hemi="Right"

ROI_list = c("FFA","FBA","Overlap_ventral")

for (roi_itr in seq_along(ROI_list)) {

  curr_ROI = ROI_list[roi_itr]
  curr_data = filter(df_all, ROI==curr_ROI, hemi==curr_hemi)
    
    
  stan_data = list( N = nrow(curr_data),
                    L = nlevels(droplevels(curr_data$subj)),
                    y = curr_data$person_PSC,
                    ll = as.numeric(droplevels(curr_data$subj)),
                    x1 = curr_data$face_PSC,
                    x2 = curr_data$body_PSC
                   )
  
  # run stan:
  
  cat ("*** ",curr_ROI, curr_hemi,": ***\n") 
  fit = stan(file='subj_model.stan', data=stan_data,
             iter=3500, warmup=1000, chains=4, seed=1,
             refresh=3500,
             control = list( adapt_delta = 0.9,max_treedepth = 10))
  
  fit_intercept = stan(file='subj_model_intercept.stan', data=stan_data,
             iter=3500, warmup=1000, chains=4, seed=1,
             refresh=3500,
             control = list( adapt_delta = 0.9,max_treedepth = 10))
  
  fit_interaction = stan(file='subj_model_only_interaction.stan', data=stan_data,
                       iter=3500, warmup=1000, chains=4, seed=1,
                       refresh=3500,
                       control = list( adapt_delta = 0.9,max_treedepth = 10))
  
  b_fit = bridge_sampler(fit)
  b_fit_intercept = bridge_sampler(fit_intercept)
  b_fit_interaction = bridge_sampler(fit_interaction)
  
  bf_fit_over_intercept = bf(b_fit, b_fit_intercept)
  bf_fit_over_interaction = bf(b_fit, b_fit_interaction)
  
  curr_bf_stats = data.frame(bf_fit_over_intercept = bf_fit_over_intercept$bf[1],
                               bf_fit_over_interaction = bf_fit_over_interaction$bf[1],
                               ROI = curr_ROI)
  
  bf_stats = union(bf_stats,curr_bf_stats)
  
}
 

