function [results] = perform_searchlight(data, settings, params)

subj_name = settings.subj_file(1:4);
% sub_folder_path = [settings.method];
results_path = [settings.path_Results];
results.subj_name = subj_name;

this_time = fix(clock);
this_time_str = [num2str(this_time(3)) '_' num2str(this_time(2)) '_' num2str(this_time(1)) '_' num2str(this_time(4)) '_' num2str(this_time(5))];

results.file_name_prefix = ['searchlight_weights_' settings.data.data_design '_' this_time_str];

results.settings = settings;
results.params = params;

if ~exist(results_path, 'dir')
    mkdir (results_path);
end


voxels_num = length(data.inds.voxel_ind);
% min_X_ind = min(data.inds.X_ind);
% max_X_ind = max(data.inds.X_ind);
% min_Y_ind = min(data.inds.Y_ind);
% max_Y_ind = max(data.inds.Y_ind);
% min_Z_ind = min(data.inds.Z_ind);
% max_Z_ind = max(data.inds.Z_ind);

volume_size = [max(data.inds.X_ind), max(data.inds.Y_ind), max(data.inds.Z_ind)];

results.alpha = nan(voxels_num,1);
results.beta = nan(voxels_num,1);
results.Rsq = nan(voxels_num,1);
results.voxel_num_for_calc = zeros(voxels_num,1);
results.mean_loc_t_vals = zeros(voxels_num,length(data.loc_t_names));
results.loc_t_vals_names = data.loc_t_names;
results.conj_mask = data.conj_mask;


A_vector_all = mean(data.A_cond_data,2);
B_vector_all = mean(data.B_cond_data,2);
C_vector_all = mean(data.C_cond_data,2);

masked_voxels = find(data.conj_mask);

% go through all possible center voxels that is included in the mask
% for voxel_itr = 1:voxels_num
% 
%     % if current center voxel is outside the mask then skip it
%     if ~data.conj_mask(voxel_itr)
%         continue
%     end    
    


for masked_voxel_itr = 1:length(masked_voxels)
    
    voxel_itr = masked_voxels(masked_voxel_itr);
    
    % curr center voxel in XYZ
    [center_X, center_Y, center_Z] = ind2sub(volume_size, voxel_itr);
    
    voxel_list_X = zeros(1,27);
    voxel_list_Y = zeros(1,27);
    voxel_list_Z = zeros(1,27);
    voxel_list_itr = 0;
    
    for X_itr = 1:3
        
        curr_X = center_X-2+X_itr;
        if curr_X < 1 || curr_X > volume_size(1)
            continue
        end
        
        for Y_itr = 1:3
            
            curr_Y = center_Y-2+Y_itr;
            if curr_Y < 1 || curr_Y > volume_size(2)
                continue
            end
            
            for Z_itr = 1:3
                
                curr_Z = center_Z-2+Z_itr;
                if curr_Z < 1 || curr_Z > volume_size(3)
                    continue
                end
                
                % update voxels list
                voxel_list_itr = voxel_list_itr + 1;
                voxel_list_X(voxel_list_itr) = curr_X;
                voxel_list_Y(voxel_list_itr) = curr_Y;
                voxel_list_Z(voxel_list_itr) = curr_Z;
                
            end
        end
    end
    
    % get voxel inds
    voxel_list_inds = sub2ind   (volume_size,...
                                voxel_list_X(1:voxel_list_itr),...
                                voxel_list_Y(1:voxel_list_itr),...
                                voxel_list_Z(1:voxel_list_itr)); 
    
    % remove voxels from list that are out of the mask
%     included_inds_of_list = find(data.conj_mask(voxel_list_inds));
%     voxel_list_inds = voxel_list_inds(included_inds_of_list);    
    voxel_list_inds = voxel_list_inds(data.conj_mask(voxel_list_inds));
    
    
    % run weight calculator
%     % get A, B and C vectors for current iteration:
    A_vector = A_vector_all(voxel_list_inds);
    B_vector = B_vector_all(voxel_list_inds);
    C_vector = C_vector_all(voxel_list_inds);

%     [b,~,~,~,stats] = regress(C_vector,[A_vector, B_vector]);
%     
%     results.alpha(voxel_itr) = b(1);
%     results.beta(voxel_itr) = b(2);
%     results.Rsq(voxel_itr) = stats(1);
    
    % calculate dot products for later calculations
    AA = dot(A_vector,A_vector);
    BB = dot(B_vector,B_vector);
    AB = dot(A_vector,B_vector);
    CA = dot(C_vector,A_vector);
    CB = dot(C_vector,B_vector);

    % calculate alpha and beta (alpha is (W-face) and beta is (W-body))
    results.alpha(voxel_itr) = (CB*AB - CA*BB) / (AB*AB - AA*BB);
    results.beta(voxel_itr) = (CA - results.alpha(voxel_itr)*AA) / AB;
    
    SS_total = sum((C_vector - mean(C_vector)).^2);
    SS_resid = sum((C_vector - (A_vector*results.alpha(voxel_itr) + B_vector*results.beta(voxel_itr))).^2);
   
    results.Rsq(voxel_itr) = (SS_total - SS_resid)/SS_total;
    
    
    % update results vectors with weights and voxel num
%     results.alpha(voxel_itr) = alpha;
%     results.beta(voxel_itr) = beta;
    
    results.voxel_num_for_calc(voxel_itr) = length(voxel_list_inds);
    
    % update results vector with the mean localizer t-vals and model % signal change 
%     results.mean_loc_t_vals(voxel_itr,:) = mean(data.loc_t_vals(voxel_list_inds,:),1);
%     results.mean_model_psc_vals(voxel_itr,:) = mean([A_vector,B_vector,C_vector],1);
  
end

%% calculate sum and diff of weights for all voxels

results.sum = results.alpha + results.beta;
results.diff_normed = (results.alpha - results.beta)./(results.alpha + results.beta);
results.center_loc_t_vals = data.loc_t_vals;
results.coords = [data.coords.X,data.coords.Y,data.coords.Z];
results.center_model_psc_vals = [A_vector_all,B_vector_all,C_vector_all];
results.roi_masks = data.roi_masks;
results.roi_masks_headers = data.roi_masks_headers;
results.conj_mask = data.conj_mask;


%% get all results in a table format

results.data_mat = [results.alpha, results.beta, results.sum, results.diff_normed,...
                    results.Rsq, results.center_loc_t_vals, results.center_model_psc_vals,...
                    results.coords, results.voxel_num_for_calc];
                
results.data_mat = results.data_mat(results.conj_mask,:);

results.data_mat_headers = {'alpha','beta','sum','diff','Rsq', 'loc_t','model_PSC','coords','voxel_num'};
                
%% filter neighboring voxels

% neighbors_mask = zeros(volume_size);
% neighbors_mask(1:2:end,1:2:end,1:2:end)=1;
% neighbors_mask = reshape(neighbors_mask, voxels_num,[]);

%% calculate results for roi

% num_of_roi = length(data.roi_masks_headers);
% 
% results.roi.roi_size = zeros(1, num_of_roi);
% results.roi.corr_r = cell(1, num_of_roi);
% results.roi.corr_p = cell(1, num_of_roi);
% results.roi.mean_alpha = zeros(1, num_of_roi);
% results.roi.mean_beta = zeros(1, num_of_roi);
% results.roi.mean_t_vals = zeros(length(data.loc_t_names), num_of_roi);
% results.roi.corr_data = cell(1, num_of_roi);
% 
% results.roi.roi_masks = data.roi_masks;

% results.roi.corr_headers = {'Wf', 'Wb', 'weights sum', 'weights diff (norm)',...
%                             results.loc_t_vals_names{:},...
% %                             'loc_face_t', 'loc_body_t', 'loc_person_t',...
% %                             'loc_face_body_t', 'loc_person_face_t', 'loc_person_bosy_t',...
%                             'model_face_psc', 'model_body_psc', 'model_person_psc',...
%                             'X','Y','Z'};

% for roi_itr = 1: num_of_roi
%     
%     results.roi.roi_size(roi_itr) = length(find(data.roi_masks(:,roi_itr)));
%     
%     alpha_roi = results.alpha(data.roi_masks(:,roi_itr));
%     beta_roi = results.beta(data.roi_masks(:,roi_itr));
%     sum_roi = alpha_roi + beta_roi;
% %     diff_roi = alpha_roi - beta_roi;
%     diff_normalized_roi = (alpha_roi - beta_roi)./(alpha_roi + beta_roi);
% %     mean_t_vals_roi = results.mean_loc_t_vals(data.roi_masks(:,roi_itr),:);
%     Rsq_roi = results.Rsq(data.roi_masks(:,roi_itr));
%     voxels_num = results.voxel_num_for_calc(data.roi_masks(:,roi_itr));
%     t_vals_roi = results.center_loc_t_vals(data.roi_masks(:,roi_itr),:);
%     model_psc_vals_roi = results.center_model_psc_vals(data.roi_masks(:,roi_itr),:);
%     coords = results.coords(data.roi_masks(:,roi_itr),:);
% 
%     % calculate correlations
%     results.roi.corr_data{roi_itr} = [alpha_roi, beta_roi, sum_roi, diff_normalized_roi, Rsq_roi, t_vals_roi, model_psc_vals_roi,coords, voxels_num];
%     [results.roi.corr_r{roi_itr}, results.roi.corr_p{roi_itr}] = corrcoef(results.roi.corr_data{roi_itr});
%     
    
    
%     % calculate means and sem
%     results.roi.mean_alpha(roi_itr) = mean(alpha_roi);
%     results.roi.mean_beta(roi_itr) = mean(beta_roi);
%     results.roi.mean_t_vals(roi_itr,:) = mean(t_vals_roi,1);
%     
%     results.roi.sem_alpha(roi_itr) = std(alpha_roi)./sqrt(results.roi.roi_size(roi_itr));
%     results.roi.sem_beta(roi_itr) = std(beta_roi)./sqrt(results.roi.roi_size(roi_itr));
%     results.roi.sem_t_vals(roi_itr,:) = std(t_vals_roi,0,1)./sqrt(results.roi.roi_size(roi_itr));
    
% end

%% save img files to plot on freesurfer

curr_vol = data.spm_vol;
% curr_vol.dt = [16 0];
% 
% curr_vol.fname = [results_path filesep results.file_name_prefix 'alpha_' subj_name '.img'];
% V = spm_write_vol(curr_vol,reshape(results.alpha, volume_size));
% 
% curr_vol.fname = [results_path filesep results.file_name_prefix 'beta_' subj_name '.img'];
% V = spm_write_vol(curr_vol,reshape(results.beta, volume_size));
% 
% curr_vol.fname = [results_path filesep results.file_name_prefix 'sum_' subj_name '.img'];
% V = spm_write_vol(curr_vol,reshape(results.sum, volume_size));
% 
% curr_vol.fname = [results_path filesep results.file_name_prefix 'diff_normed_' subj_name '.img'];
% V = spm_write_vol(curr_vol,reshape(results.diff_normed, volume_size));

% curr_vol.fname = [results_path filesep results.file_name_prefix 'loc_face_t_' subj_name '.img'];
% V = spm_write_vol(curr_vol,reshape(results.mean_loc_t_vals(:,1), volume_size));
% 
% curr_vol.fname = [results_path filesep results.file_name_prefix 'loc_body_t_' subj_name '.img'];
% V = spm_write_vol(curr_vol,reshape(results.mean_loc_t_vals(:,2), volume_size));
% 
% curr_vol.fname = [results_path filesep results.file_name_prefix 'loc_person_t_' subj_name '.img'];
% V = spm_write_vol(curr_vol,reshape(results.mean_loc_t_vals(:,3), volume_size));

results.curr_vol = curr_vol;
%% save results
save ([results_path filesep results.file_name_prefix '_' subj_name '.mat'], 'results');

