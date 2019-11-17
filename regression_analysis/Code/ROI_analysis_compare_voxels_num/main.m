function [data,settings, params, results, result_tables] = main(settings, params, load_setting_file_name)

addpath ('/data/libi/fMRI_analyses/gyeuclidean/Code');

if isempty(load_setting_file_name)
    [settings, params] = load_settings_params(settings, params);
    
else
    currDir = pwd;
    [pathstr,name,ext] = fileparts(load_setting_file_name);
    cd (pathstr);
    [settings, params] = feval(name, settings, params);
    cd (currDir);
        
end


   
rng(params.seed);
% addpath (settings.path_spm8);

% load mask and relevant pairs of conditions (data from all runs)
data = load_data(settings);

% prepare data for analysis
data = divide_data_to_conditions(settings, data);
data = create_labels_for_train_and_cv(settings, data);

result_tables.table_no_intercept = [];
result_tables.table_with_intercept = [];
result_tables.table_compare_models = [];

orig_data = data;

for voxels_num = [5,10,15,20,25,30,35,40,50,60,70,80,100,120]
    
    settings.min_voxel_num = voxels_num;
    settings.exact_voxel_num = voxels_num;
    
    data = calculate_roi_masks(orig_data, settings);

    % perform analysis

    results = perform_classification_and_cv_with_intercept(data, settings, params);
    
    result_tables.table_no_intercept = [result_tables.table_no_intercept; results.table_no_intercept];
    result_tables.table_with_intercept = [result_tables.table_with_intercept; results.table_with_intercept];
    result_tables.table_compare_models = [result_tables.table_compare_models; results.table_compare_models];
    
end
end


function [settings, params] = load_settings_params()
    [FileName,PathName,FilterIndex] = uigetfile('CreateSettingsParams.m');
    
    currDir = pwd;
    cd (PathName);
    funcName = FileName(1:(end-2));
    [settings, params] = feval(funcName);
    cd (currDir);
end
