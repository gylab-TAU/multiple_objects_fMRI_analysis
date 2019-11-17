clear; clc; close all;


subj = dir ('/data/libi/MRI_data/Face_Body_Integration/data_mat_files_m/*.mat');

load_settings_file_list = { './load_settings_params_exp1.m'};



MAX_SUBJ_NUM = 20;


table_no_intercept = [];
table_with_intercept = [];
table_compare_models = [];


for subj_itr = 1:min(length(subj), MAX_SUBJ_NUM)
    
    for settings_file_itr = 1: length(load_settings_file_list)
        
        settings.subj_file = subj(subj_itr).name;
        params = [];
        [data,settings, params, results, result_tables] = main (settings, params, load_settings_file_list{settings_file_itr}); 
        
         table_no_intercept = [table_no_intercept; result_tables.table_no_intercept];
         table_with_intercept = [table_with_intercept; result_tables.table_with_intercept];
         table_compare_models = [table_compare_models; result_tables.table_compare_models];
         
        fprintf('\n');

    end
end


t = table_no_intercept(table_no_intercept.voxels_num>0, :);
writetable(t, [settings.path_Results filesep 'table_no_intercept.csv'],'WriteRowNames',true);
% writetable(table_with_intercept, [settings.path_Results filesep 'table_with_intercept.csv'],'WriteRowNames',true);
% writetable(table_compare_models, [settings.path_Results filesep 'table_compare_models.csv'],'WriteRowNames',true);
