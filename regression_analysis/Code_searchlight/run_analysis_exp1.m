clear; clc; close all;


subj = dir ('/data/libi/MRI_data/Face_Body_Integration/data_mat_files_m/*.mat');

load_settings_file_list = {'./load_settings_params_searchlight_exp1.m'



MAX_SUBJ_NUM = 20;


for subj_itr = 1:min(length(subj), MAX_SUBJ_NUM)
    
    for settings_file_itr = 1: length(load_settings_file_list)
        
        settings.subj_file = subj(subj_itr).name;
        params = [];
        [data,settings, params] = main_searchlight (settings, params, load_settings_file_list{settings_file_itr}); 
        fprintf('\n');

    end
end
    