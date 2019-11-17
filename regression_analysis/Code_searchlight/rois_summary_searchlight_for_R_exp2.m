
results_dir_FB = './Results/searchlight_FB';
results_dir_FO = './Results/searchlight_FO';


files_FB = dir ([results_dir_FB filesep '*.mat']);
files_FO = dir ([results_dir_FO filesep '*.mat']);

load([results_dir_FB filesep files_FB(1).name]);

roi_num = length(results.roi.roi_size);

subj_num = length(files_FB);
voxel_num = length(results.alpha);

roi_names = {'FBO_Left',  'FBO_Right', 'EVC','CSF'};


roi_loc = {'FBO', 'FBO','EVC','CSF'};

rois_hemis = {'Left','Right','both','both'};

results_cell = {'Wf_FB', 'Wb_FB', 'sum_FB', 'diff_FB', 'Rsq_FB', 'face_PSC_FB', 'body_PSC_FB', 'person_PSC_FB',...
                'Wf_FO', 'Wo_FO', 'sum_FO', 'diff_FO', 'Rsq_FO', 'face_PSC_FO', 'object_PSC_FO', 'integrated_PSC_FO',...    
                'face_others_t_loc', 'body_others_t_loc', 'object_others_t_loc',...
                'face_body_t_loc', 'face_object_t_loc', ...
                'X','Y','Z','num_of_voxels','ROI', 'subj', 'voxel_ind'};

for subj_itr = 1:length(files_FB)
    
    load([results_dir_FB filesep files_FB(subj_itr).name]);    
    results_FB = results;
    
    results_FO_file = dir ([ results_dir_FO filesep '*' results_FB.subj_name '.mat']);
    load([results_dir_FO filesep results_FO_file(1).name]);    
    results_FO = results;
    
    if ~strcmp(results_FB.subj_name, results_FO.subj_name)
       fprintf ('different subjects!!') 
       error = 1;
    end
        
        
    for roi_itr = 1:3 %length(results.roi_masks_headers)
        curr_roi_size = results.roi.roi_size(roi_itr);
        data_cell_FB = num2cell(results_FB.roi.corr_data{roi_itr});
        data_cell_FO = num2cell(results_FO.roi.corr_data{roi_itr});
        subj_name_cells = repmat ({results_FB.subj_name},curr_roi_size,1);
        roi_name_cells = repmat (roi_names(roi_itr),curr_roi_size,1);
        
        curr_roi_cell = [data_cell_FB(:,1:5),...    % Wf,Wb,summ,diff,Rsq- FB
                        data_cell_FB(:,22:24),...   % PSC FB
                        data_cell_FO(:,1:5),...     % Wf,Wo,summ,diff,Rsq - FO
                        data_cell_FO(:,22:24),...   % PSC FO
                        data_cell_FB(:,6:21),...    % t values
                        data_cell_FB(:,24:27),...   % coordinates and num of voxels
                        roi_name_cells,...
                        subj_name_cells];
                    
        results_cell = [results_cell; curr_roi_cell];
    end
    
    
end



file_name = [results_dir_FB filesep 'searchlite_rois_summary_for_R_FB_FO.csv'];
    cell2csv(file_name, results_cell);




