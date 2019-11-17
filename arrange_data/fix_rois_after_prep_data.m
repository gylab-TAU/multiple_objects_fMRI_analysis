function [] = fix_rois_after_prep_data(PARAMS)



                

for subj_itr = 1:length(PARAMS.subjects_list)
    
    
    subj_file_name = [PARAMS.data_dir filesep PARAMS.subjects_list{subj_itr}];
    fprintf('Subject %i - %s \n',subj_itr, PARAMS.subjects_list{subj_itr}); 
    
    load(subj_file_name);
%     headers = subj_data.header;
    
    curr_subj_dir = [PARAMS.SpmDir PARAMS.subjects_list{subj_itr} filesep];
    curr_roi_dir = [curr_subj_dir PARAMS.unite_rois.new_rois_dir];
    
    voxels_num = size(subj_data.data,1);
    
    for roi_itr = 1:length(PARAMS.all_rois)
       
        curr_roi_name = PARAMS.all_rois{roi_itr};
        roi_files = dir ([curr_roi_dir filesep curr_roi_name(1:end-4) '*.mat']);
        
        % handle midding rois
        if isempty(roi_files)
            
            new_roi_ind = length(subj_data.header)+1;
            subj_data.data(:,new_roi_ind) = zeros(voxels_num,1);        
            subj_data.header{new_roi_ind} = curr_roi_name;

            
        % handle rois with multiple clusters 
        elseif length(roi_files) > 1
           
            new_roi_ind = length(subj_data.header)+1;
            subj_data.header{new_roi_ind} = curr_roi_name;

            subj_data.data(:,new_roi_ind) = zeros(voxels_num,1);

            for files_itr = 1:length(roi_files)
                subj_data.data(:,new_roi_ind) = subj_data.data(:,new_roi_ind) | subj_data.data(:,find(strcmp(roi_files(files_itr).name(1:end-4),subj_data.header)));
            end
                        
        end
        
    end
    
    % unite m_FFA and p_FFA:
    
    % left
    new_roi_ind = length(subj_data.header)+1;
    subj_data.header{new_roi_ind} = 'FFA_left_roi';
    subj_data.data(:,new_roi_ind) = subj_data.data(:,find(strcmp('p_FFA_left_roi',subj_data.header))) | subj_data.data(:,find(strcmp('m_FFA_left_roi',subj_data.header)));    
    
    % right
    new_roi_ind = length(subj_data.header)+1;
    subj_data.header{new_roi_ind} = 'FFA_right_roi';
    subj_data.data(:,new_roi_ind) = subj_data.data(:,find(strcmp('p_FFA_right_roi',subj_data.header))) | subj_data.data(:,find(strcmp('m_FFA_right_roi',subj_data.header)));
    
        
    save(subj_file_name, 'subj_data');
   
end
    
%     % missing_rois
%     missing_rois = PARAMS.fix_rois.missing_rois{subj_itr};
%     if ~isempty(missing_rois)
%         missing_rois_num = length(missing_rois);
%         new_roi_ind = length(subj_data.header)+1;
%         
%         subj_data.data(:,(new_roi_ind:(new_roi_ind+missing_rois_num-1))) = zeros(size(subj_data.data,1),missing_rois_num);
%         for roi_itr = 1:missing_rois_num
%             subj_data.header{new_roi_ind} = missing_rois{roi_itr};
%             new_roi_ind = new_roi_ind+1;
%         end
%     end

    
    %unite_rois
%     unite_rois_cells = PARAMS.fix_rois.unite_rois{subj_itr};
%     if ~isempty(unite_rois_cells)
%         for unite_cell_itr = 1:length(unite_rois_cells)
%             unite_rois = unite_rois_cells{unite_cell_itr};
%             new_roi_ind = length(subj_data.header)+1;
%             subj_data.header{new_roi_ind} = unite_rois{1};
% 
%             subj_data.data(:,new_roi_ind) = subj_data.data(:,find(strcmp(unite_rois{2},subj_data.header)));
% 
%             for rois_itr = 3:length(unite_rois)
%                 subj_data.data(:,new_roi_ind) = subj_data.data(:,new_roi_ind) | subj_data.data(:,find(strcmp(unite_rois{rois_itr},subj_data.header)));
%             end
%         end
%         
%     end
    


        
        
