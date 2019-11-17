
results_dir = './Results/overlap_30_most_selective';

files = dir ([results_dir filesep '*.mat']);

load([results_dir filesep files(1).name]);
results_cell = {'face_PSC', 'body_PSC', 'person_PSC',...
                'subj', 'ROI', 'hemi'};
line_ind = 2;
roi_num = length(results.roi_name);
rois = {'FFA','FBA','Overlap_ventral', 'FFA','FBA','Overlap_ventral', 'OFA','EBA','Overlap_lat_occ','OFA','EBA','Overlap_lat_occ'}';
hemis = {'Left', 'Left', 'Left', 'Right', 'Right', 'Right','Left', 'Left', 'Left', 'Right', 'Right', 'Right'}';

for file_itr = 1:length(files)
    
    load([results_dir filesep files(file_itr).name]);
    for roi_itr = 1:length(results.roi_size)
        if (results.roi_size(roi_itr) > 0)
            
            face_PSC = num2cell(results.vals_psc_A_cond{roi_itr});
            body_PSC = num2cell(results.vals_psc_B_cond{roi_itr});
            person_PSC = num2cell(results.vals_psc_C_cond{roi_itr});
            
            roi_name_cells = repmat (rois(roi_itr),results.roi_size(roi_itr),1);
            roi_hemi_cells = repmat (hemis(roi_itr),results.roi_size(roi_itr),1);
            subj_name_cells = repmat({results.subj_name},results.roi_size(roi_itr),1);
            
            curr_cell = [face_PSC, body_PSC, person_PSC,...
                        subj_name_cells, roi_name_cells, roi_hemi_cells];
            
            results_cell = [results_cell ; curr_cell];
            
        end
        
    end
end    


file_name = [results_dir filesep 'rois_summary_data_for_R_all_together_single_voxels_m_30.csv'];
cell2csv(file_name, results_cell);