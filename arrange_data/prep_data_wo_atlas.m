function []=prep_data_wo_atlas(PARAMS)



%%
addpath (PARAMS.path_spm);
addpath (PARAMS.path_marsbar);
mkdir (PARAMS.data_dir);

%%
subj_all_dirs_names = dir(PARAMS.SpmDir);
if ~isempty(PARAMS.subjects_list)
    num_subj_dirs_total = length(PARAMS.subjects_list);
else
    num_subj_dirs_total = length(subj_all_dirs_names);
end

fprintf('Function ''prepData'' - extracting data for subjects in folder %s\n\n', PARAMS.SpmDir);


% Go over all folders/subjects in directory and extract the data for them.
% Assuming that the directory includes only folders for subjects with the
% appropriate data and structure inside.
subj_count = 0; % Counting subjects (and not directories)

for num_subj_ind = 1 : num_subj_dirs_total
    if ~isempty(PARAMS.subjects_list) && subj_count<length(PARAMS.subjects_list)        
        subj_name = PARAMS.subjects_list{num_subj_ind};
        subj_count = subj_count + 1;
    else
        this_subj_dir = subj_all_dirs_names(num_subj_ind);
        if subj_count+1<=PARAMS.max_subj_num && length(this_subj_dir.name) > 2 && strcmp(this_subj_dir.name(1:4), PARAMS.roi_anal_prefix_dir_name) == 0 % not '.' or '..' or 'ROI_Analysis' directories
            subj_count = subj_count + 1;
            subj_name = this_subj_dir.name;
        else 
            continue
        end
    end
    

            
        % Extract data for subjects
        fprintf('Extracting voxels data for subject %i - %s \n',subj_count, subj_name);     
        
        subj_spm_path = [PARAMS.SpmDir subj_name filesep];
        beta_headers = cell(1,length(PARAMS.design_dirs));
        percent_signal_change_headers = cell(1,length(PARAMS.design_dirs));
        mean_across_sess_beta_headers = cell(1,length(PARAMS.design_dirs));
        mean_across_sess_percent_signal_change_headers = cell(1,length(PARAMS.design_dirs));
        brain_activation_masks_headers = cell(1,length(PARAMS.design_dirs));
        contrast_headers = cell(1,length(PARAMS.design_dirs));
        brain_mask_values = cell(1,length(PARAMS.design_dirs));
        beta_values = cell(1,length(PARAMS.design_dirs));
        percent_signal_change_values = cell(1,length(PARAMS.design_dirs));
        contrasts_values = cell(1,length(PARAMS.design_dirs));
        mean_across_sess_beta_values = cell(1,length(PARAMS.design_dirs));
        mean_across_sess_percent_signal_change_values = cell(1,length(PARAMS.design_dirs));
        
        for design_ind = 1:length(PARAMS.design_dirs)
        
            fprintf('Extracting %s data\n',PARAMS.design_dirs{design_ind});
            
            % load SPM.mat
            subj_design_path = [subj_spm_path PARAMS.design_dirs{design_ind}];
            subj_spm_file = [subj_design_path filesep 'SPM.mat'];         
            load(subj_spm_file);            
            
            % get num of sessions, num of conditions and condition names
            % from SPM.matcontrasts_values
            num_sess = length(SPM.Sess);     
            num_conds = length(SPM.Sess(1).Fc);
            num_motion_regs = zeros(1,num_sess);
            for sess_itr = 1:num_sess
                num_motion_regs(sess_itr) = length(SPM.Sess(sess_itr).C.name);
            end
            [cond_names{1:num_conds}] = deal(SPM.Sess(1).Fc(:).name);
            
            % get contrasts names and weights from SPM.mat
            num_contrasts = length(SPM.xCon);
            if num_contrasts
                [contrasts_names{1:num_contrasts}] = deal(SPM.xCon(:).name);
                [contrasts_weights{1:num_contrasts}] = deal(SPM.xCon(:).c);
            end
            
            fprintf('\t Found %i conditions in %i sessions and %i contrasts \n',num_conds, num_sess, num_contrasts);     
            
            % for later mat file saving:
            subj_data.designs{design_ind}.path = PARAMS.design_dirs{design_ind};
            subj_data.designs{design_ind}.header_prefix = PARAMS.header_prefix{design_ind};
            subj_data.designs{design_ind}.num_sess = num_sess;
            subj_data.designs{design_ind}.num_conds = num_conds;
            subj_data.designs{design_ind}.conds_names = cond_names;
            subj_data.designs{design_ind}.num_contrasts = num_contrasts;
            if num_contrasts               
                subj_data.designs{design_ind}.contrasts_names = contrasts_names;
                subj_data.designs{design_ind}.contrasts_weights = contrasts_weights;
            end
            subj_data.designs{design_ind}.SPM = SPM;
            
            % beta and percent signal change headers
            beta_headers{design_ind} = cell(1,num_sess * num_conds);
            percent_signal_change_headers{design_ind} = cell(1,num_sess * num_conds);
            header_ind = 1;
            for sess_ind = 1:num_sess
                for cond_ind = 1:num_conds
                    beta_headers{design_ind}{header_ind} = ['Beta_' PARAMS.header_prefix{design_ind} '_'  cond_names{cond_ind} '_Sn' num2str(sess_ind)];
                    percent_signal_change_headers{design_ind}{header_ind} = ['PercSigCh_' PARAMS.header_prefix{design_ind} '_'  cond_names{cond_ind} '_Sn' num2str(sess_ind)];
                    header_ind = header_ind + 1;
                    
                end
            end
            
            
            % mean across sessions beta and percent signal change headers
            mean_across_sess_beta_headers{design_ind} = cell(1,num_conds);
            mean_across_sess_percent_signal_change_headers{design_ind} = cell(1,num_conds);
            header_ind = 1;
            for cond_ind = 1:num_conds
                mean_across_sess_beta_headers{design_ind}{header_ind} = ['Beta_' PARAMS.header_prefix{design_ind} '_'  cond_names{cond_ind} '_mean_across_all_sess'];
                mean_across_sess_percent_signal_change_headers{design_ind}{header_ind} = ['PercSigCh_' PARAMS.header_prefix{design_ind} '_'  cond_names{cond_ind} '_mean_across_all_sess'];
                header_ind = header_ind + 1;

            end
           
            % brain activation masks headers
            brain_activation_masks_headers{design_ind} = ['Brain_activation_mask_' PARAMS.header_prefix{design_ind}];
            
            
            % contrasts headers
            if num_contrasts
                contrast_headers{design_ind} = cell(1,num_contrasts);
                for contrast_ind = 1:num_contrasts
                    contrast_headers{design_ind}{contrast_ind} = ['Cont_' PARAMS.header_prefix{design_ind} '_'  contrasts_names{contrast_ind}];           
                end
                
            else
                contrast_headers{design_ind} = [];
            end
  
            
            % Extract betas
            fprintf('\t Extracting %s brain activation mask \n', PARAMS.design_dirs{design_ind}); 

            % Extract mask
            mask_file_name = [subj_design_path filesep 'mask.nii'];
            mask_file_name = strvcat(mask_file_name);
            
            mask_vols = spm_vol(mask_file_name);
            [brain_mask_values{design_ind}, mask_XYZ] = spm_read_vols(mask_vols);

            % save the properties of the spm volumes for later writing new volumes
            % with the results
            if design_ind == 1
                subj_data.spm_vol = struct('fname',   [],...
                                           'dim',     mask_vols.dim,...
                                           'dt',      mask_vols.dt,...
                                           'pinfo',   mask_vols.pinfo,...
                                           'mat',     mask_vols.mat,...
                                           'n',       mask_vols.n,...
                                           'descrip', [],...
                                           'private', []);
            end

            % Get segmentation maps
            fprintf('Get segmntation maps\n');
            curr_anatomy_path = [subj_spm_path filesep PARAMS.anatomy_path{num_subj_ind}];
            segmentation_files = dir ([curr_anatomy_path filesep 'rc*.nii']);
            
            segmentation_files_names = cell(1,3);
%             segmentation_headers = cell(1,length(segmentation_files));
            
            for seg_file_itr = 1:3
                segmentation_files_names{seg_file_itr} = [curr_anatomy_path filesep segmentation_files(seg_file_itr).name];
%                 segmentation_roi_headers{seg_file_itr} = segmentation_files(seg_file_itr).name(3:end-4);
            end
            
            segmentation_headers = {'GM','WM','CSF'};
            segmentation_files_names = strvcat(segmentation_files_names);
            
            segmentation_vols = spm_vol(segmentation_files_names);
            segmentation_values = spm_read_vols(segmentation_vols);
            
            
            fprintf('\t Extracting %s betas \n', PARAMS.design_dirs{design_ind}); 
            conds_betas_num = num_sess * num_conds;
            total_betas_num = conds_betas_num + sum(num_motion_regs); % without the last betas for session mean
            beta_files_names = cell(1,conds_betas_num);
            cond_beta_count =1;
            spm_beta_counter = 0;
            for session_itr = 1:num_sess
                for cond_itr = spm_beta_counter + (1:num_conds)
                  
                    beta_files_names{cond_beta_count} = [subj_design_path filesep 'beta_', num2str(cond_itr,'%.4d'),'.nii'];
                    cond_beta_count = cond_beta_count+1;
                    
                end
                
                spm_beta_counter = spm_beta_counter +num_conds + num_motion_regs(session_itr);
            end
            


            beta_files_names = strvcat(beta_files_names);
            beta_vols = spm_vol(beta_files_names);
            [beta_values{design_ind}, XYZ] = spm_read_vols(beta_vols);
            
            
            % Calculate percent signal change
            fprintf('\t Calculating %s percent signal change \n', PARAMS.design_dirs{design_ind}); 
                       
            % Extract betas for session means (the last betas of the design matrix)                
            sessions_betas_files_names = cell(1,num_sess);
            for beta_ind = 1:num_sess  % The betas for each sessions are the last betas at the design matrix.
                sessions_betas_files_names{beta_ind} = [subj_design_path filesep 'beta_', num2str(beta_ind + total_betas_num,'%.4d'),'.nii']; % The betas for each sessions are the last betas at the design matrix.
            end

            sessions_betas_files_names = strvcat(sessions_betas_files_names);
            beta_sessions_vols = spm_vol(sessions_betas_files_names);
            beta_sessions_values = spm_read_vols(beta_sessions_vols);


            dupl_beta_sessions_values = zeros(size(beta_values{design_ind}));
            for sess_ind = 1:num_sess
                dupl_beta_sessions_values(:,:,:,(1:num_conds)+(num_conds*(sess_ind-1))) = repmat(beta_sessions_values(:,:,:,sess_ind), 1,1,1,num_conds);
            end

            % dupl_beta_sessicontrasts_valuesons_values{design_ind} = repelem(beta_sessions_values{design_ind},1,1,1,num_conds);
            % the above line replaces the for loop but works only in matlabr2015a and above.

            max_regressor_height = max(SPM.xX.X(:,1)); % we are taking the max value from the first regressor (first col) in the regression matrix

            percent_signal_change_values{design_ind} = 100 * beta_values{design_ind} * max_regressor_height ./ dupl_beta_sessions_values;
            
            
            % Extract contrasts values
            if num_contrasts
                fprintf('\t Extracting %s contrasts \n', PARAMS.design_dirs{design_ind}); 
                contrasts_files_names = cell(1,num_contrasts);
                for contrast_ind = 1:num_contrasts
                    contrasts_files_names{contrast_ind} = [subj_design_path filesep 'spmT_', num2str(contrast_ind,'%.4d'),'.nii'];
                end

                contrasts_files_names = strvcat(contrasts_files_names);
                contrasts_vols = spm_vol(contrasts_files_names);
                contrasts_values{design_ind} = spm_read_vols(contrasts_vols);
            end
 
            
            % Calculate mean values for conditions across all sessions
            % for betas and percent signal change
            fprintf('\t Calculating %s mean values for beta and percent signal change across sessions\n', PARAMS.design_dirs{design_ind});
            [size_x, size_y, size_z, ~] = size(beta_values{design_ind});
            mean_across_sess_beta_values{design_ind} = zeros(size_x, size_y, size_z, num_conds);
            mean_across_sess_percent_signal_change_values{design_ind} = zeros(size_x, size_y, size_z, num_conds);
            
            for cond_ind = 1:num_conds
               curr_cond_inds = (cond_ind : num_conds : num_sess*num_conds);
               
               mean_across_sess_beta_values{design_ind}(:,:,:,cond_ind) = mean(beta_values{design_ind}(:,:,:,curr_cond_inds),4);
               mean_across_sess_percent_signal_change_values{design_ind}(:,:,:,cond_ind) = mean(percent_signal_change_values{design_ind}(:,:,:,curr_cond_inds),4);
            end
            
        end
       
        % Get ROIs (marsbar)
        fprintf('Get ROIs from marsbar\n') 
        curr_roi_path = [subj_spm_path filesep PARAMS.rois_dir];
        roi_files = dir ([curr_roi_path filesep '*.mat']);
        roi_data = cell(1, length(roi_files));
        roi_headers = cell(1, length(roi_files));
        
        for roi_ind = 1:length(roi_files)
        
            roi_headers{roi_ind} = roi_files(roi_ind).name(1:(end-4));        
            roi_obj = maroi([curr_roi_path filesep roi_files(roi_ind).name]);
            [~, ~, roi_XYZ, ~] = getdata(roi_obj,mask_vols);
            
            % temp - ver. 1
%             roi_mask = zeros(size(brain_mask_values{1}));
%             for vox_num = 1: size (roi_XYZ,2)
%                 roi_mask(uint8(roi_XYZ(1,vox_num)),uint8(roi_XYZ(2,vox_num)),uint8(roi_XYZ(3,vox_num))) = 1;
%             end
%             
            % temp - ver. 2
            roi_XYZ = uint8(roi_XYZ);
            roi_XYZ_ind = sub2ind(size(brain_mask_values{1}), roi_XYZ(1,:), roi_XYZ(2,:), roi_XYZ(3,:));
            roi_data{roi_ind} = zeros(numel(brain_mask_values{1}),1);
            roi_data{roi_ind}(uint64(roi_XYZ_ind)) = 1;
            
        end 
            
        % Get atlas masks
%         fprintf('Get atlas parcellation masks\n');
%         atlas_path = [subj_spm_path filesep 'Atlas_rois'];
%         atlas_files = dir ([atlas_path filesep 'rwROI*.nii']);
%         
%         atlas_files_names = cell(1,length(atlas_files));
%         atlas_roi_headers = cell(1,length(atlas_files));
%         
%         for atlas_file_itr = 1:length(atlas_files)
%             atlas_files_names{atlas_file_itr} = [atlas_path filesep atlas_files(atlas_file_itr).name];
%             atlas_roi_headers{atlas_file_itr} = atlas_files(atlas_file_itr).name(3:end-4);
%         end
%         
%         atlas_files_names = strvcat(atlas_files_names);
% 
%         atlas_vols = spm_vol(atlas_files_names);
%         atlas_values = spm_read_vols(atlas_vols);
%         
        
        % Get FreeSurfer parcellation masks                
        if PARAMS.use_fressurfer
            fprintf('Get FreeSurfer parcellation masks\n') 
            subj_freesurfer_path = [PARAMS.FreesurferDir subj_name filesep];

            aparc_a2009_file_names = { [subj_freesurfer_path 'aparc_a2009_mask_lh.img']
                                        [subj_freesurfer_path 'aparc_a2009_mask_rh.img']};
            aparc_file_names = { [subj_freesurfer_path 'aparc_mask_lh.img']
                                        [subj_freesurfer_path 'aparc_mask_rh.img']};

            aparc_a2009_file_names = strvcat(aparc_a2009_file_names);
            aparc_a2009_vols = spm_vol(aparc_a2009_file_names);
            aparc_a2009_values = spm_read_vols(aparc_a2009_vols);

            aparc_file_names = strvcat(aparc_file_names);
            aparc_vols = spm_vol(aparc_file_names);
            aparc_values = spm_read_vols(aparc_vols);        

            aparc_a2009_ctab_file_names = [subj_freesurfer_path 'label' filesep 'aparc.annot.a2009s.ctab'];
            aparc_a2009_ctab_fileId = fopen(aparc_a2009_ctab_file_names);
            aparc_a2009_ctab_fileData = textscan(aparc_a2009_ctab_fileId, '%d %s %d %d %d %d'); 
            fclose(aparc_a2009_ctab_fileId);
            subj_data.parcellation.aparc_a2009_labels = deal(aparc_a2009_ctab_fileData{2}); 
            subj_data.parcellation.aparc_a2009_numbers = deal(aparc_a2009_ctab_fileData{1});


            aparc_ctab_file_names = [subj_freesurfer_path 'label' filesep 'aparc.annot.ctab'];
            aparc_ctab_fileId = fopen(aparc_ctab_file_names);
            aparc_ctab_fileData = textscan(aparc_ctab_fileId, '%d %s %d %d %d %d'); 
            fclose(aparc_ctab_fileId);
            subj_data.parcellation.aparc_labels = deal(aparc_ctab_fileData{2}); 
            subj_data.parcellation.aparc_numbers = deal(aparc_ctab_fileData{1});
        end
        
        % Write to csv array
        fprintf('Saving to files \n'); 
        
        % get voxel inds
        vol_size = size(beta_values{1});
        vol_size = vol_size(1:3);
        num_voxels = size(XYZ,2);
        
        [X_ind, Y_ind, Z_ind] = ind2sub(vol_size, 1:num_voxels);
                
        

        % reshape voxels values from a 4D mat (with XTZ location inds) to
        % indexed listsbetas
%         atlas_values_list = reshape(atlas_values, num_voxels, []);
        
        beta_values_lists = cell(1, length(PARAMS.design_dirs));
        percent_signal_change_values_lists = cell(1, length(PARAMS.design_dirs));
        mean_across_sess_beta_values_lists = cell(1, length(PARAMS.design_dirs));
        mean_across_sess_percent_signal_change_values_lists = cell(1, length(PARAMS.design_dirs));
        contrasts_values_lists = cell(1, length(PARAMS.design_dirs));
        brain_activation_mask_lists = cell(1, length(PARAMS.design_dirs));
                
        
        if PARAMS.use_fressurfer
            aparc_a2009_values_lh_list = reshape(aparc_a2009_values(:,:,:,1),num_voxels, []);
            aparc_a2009_values_rh_list = reshape(aparc_a2009_values(:,:,:,2),num_voxels, []);
            aparc_values_lh_list = reshape(aparc_values(:,:,:,1),num_voxels, []);
            aparc_values_rh_list = reshape(aparc_values(:,:,:,2),num_voxels, []);
        end
        
        for design_ind = 1:length(PARAMS.design_dirs)
            beta_values_lists{design_ind} = reshape(beta_values{design_ind},num_voxels, []);
            percent_signal_change_values_lists{design_ind} = reshape(percent_signal_change_values{design_ind},num_voxels, []);
            mean_across_sess_beta_values_lists{design_ind} = reshape(mean_across_sess_beta_values{design_ind},num_voxels, []);
            mean_across_sess_percent_signal_change_values_lists{design_ind} = reshape(mean_across_sess_percent_signal_change_values{design_ind},num_voxels, []);
            brain_activation_mask_lists{design_ind} = reshape(brain_mask_values{design_ind}, num_voxels,[]);
            
            if subj_data.designs{design_ind}.num_contrasts
                contrasts_values_lists{design_ind} = reshape( contrasts_values{design_ind},num_voxels,[]);
            else
                contrasts_values_lists{design_ind} = [];
            end
        end
        
        segmentation_values_list = reshape( segmentation_values,num_voxels,[]);
        
        % empty cells from contrasts_values_list (they are empty because
        % some designs doesn't have contrasts
        contrasts_values_lists = contrasts_values_lists(~cellfun(@isempty,contrasts_values_lists));
        
        
        % perform conjunction between model and loclizer brain masks       
        brain_activation_conj_mask_list = brain_activation_mask_lists{1};        
        for design_ind = 2:length(brain_activation_mask_lists)
            brain_activation_conj_mask_list = brain_activation_conj_mask_list & brain_activation_mask_lists{design_ind};
        end
        
        
        
        % put all data at a single mat files with rows=voxels and
        % columns=features
        if PARAMS.use_fressurfer            
            subj_data.data = [(1:num_voxels)', X_ind', Y_ind', Z_ind', XYZ',beta_values_lists{:}, percent_signal_change_values_lists{:},...
                                contrasts_values_lists{:}, mean_across_sess_beta_values_lists{:}, mean_across_sess_percent_signal_change_values_lists{:},...
                                brain_activation_mask_lists{:}, brain_activation_conj_mask_list, segmentation_values_list, roi_data{:},...
                                aparc_a2009_values_lh_list, aparc_a2009_values_rh_list, aparc_values_lh_list, aparc_values_rh_list];
            
        else 
            subj_data.data = [(1:num_voxels)', X_ind', Y_ind', Z_ind', XYZ',beta_values_lists{:}, percent_signal_change_values_lists{:},...
                                contrasts_values_lists{:}, mean_across_sess_beta_values_lists{:}, mean_across_sess_percent_signal_change_values_lists{:},...
                                brain_activation_mask_lists{:}, brain_activation_conj_mask_list, segmentation_values_list, roi_data{:}];
    %                             aparc_a2009_values_lh_list, aparc_a2009_values_rh_list, aparc_values_lh_list, aparc_values_rh_list,...
        end
                        
%         subj_data.data = subj_data.data(brain_activation_mask_inds,:);
        
       
        if PARAMS.save_csv    
            fprintf('Writing to csv file \n');
            subj_csv_file_name = [PARAMS.data_dir subj_name '.csv'];
            csvwrite(subj_csv_file_name, subj_data.data);
        end
        
        if PARAMS.save_mat
            fprintf('Writing to mat file \n');
            beta_headers_list = {};
            percent_signal_change_headers_list = {};
            mean_across_sess_beta_headers_list = {};
            mean_across_sess_percent_signal_change_headers_list = {};
            contrast_headers_list = {};
            brain_activation_masks_headers_list = {};
            
            for design_ind = 1: length(PARAMS.design_dirs)
                beta_headers_list = [ beta_headers_list, beta_headers{design_ind}];
                percent_signal_change_headers_list = [percent_signal_change_headers_list, percent_signal_change_headers{design_ind}];
                mean_across_sess_beta_headers_list = [mean_across_sess_beta_headers_list, mean_across_sess_beta_headers{design_ind}];
                mean_across_sess_percent_signal_change_headers_list = [mean_across_sess_percent_signal_change_headers_list, mean_across_sess_percent_signal_change_headers{design_ind}];
                brain_activation_masks_headers_list = [brain_activation_masks_headers_list, brain_activation_masks_headers{design_ind}];
                if subj_data.designs{design_ind}.num_contrasts
                    contrast_headers_list = [contrast_headers_list, contrast_headers{design_ind}];
                end
            end
            
            if PARAMS.use_fressurfer 
                subj_data.header = {'voxel_ind', 'X_ind', 'Y_ind', 'Z_ind', 'X_coord', 'Y_coord', 'Z_coord', ...
                                beta_headers_list{:}, percent_signal_change_headers_list{:}, contrast_headers_list{:}, ...
                                mean_across_sess_beta_headers_list{:}, mean_across_sess_percent_signal_change_headers_list{:}, ...
                                brain_activation_masks_headers_list{:}, 'Brain_activation_conj_mask', segmentation_headers{:},roi_headers{:}, ...
                                'aparc_a2009_lh', 'aparc_a2009_rh', 'aparc_lh', 'aparc_rh'} ;
            else
                subj_data.header = {'voxel_ind', 'X_ind', 'Y_ind', 'Z_ind', 'X_coord', 'Y_coord', 'Z_coord', ...
                                beta_headers_list{:}, percent_signal_change_headers_list{:}, contrast_headers_list{:}, ...
                                mean_across_sess_beta_headers_list{:}, mean_across_sess_percent_signal_change_headers_list{:}, ...
                                brain_activation_masks_headers_list{:}, 'Brain_activation_conj_mask', segmentation_headers{:},roi_headers{:}} ;
%                                 'aparc_a2009_lh', 'aparc_a2009_rh', 'aparc_lh', 'aparc_rh'} ;
            end
            
            subj_data.header = strrep(subj_data.header,'+','_');
            subj_mat_file_name = [PARAMS.data_dir filesep subj_name '.mat'];
            save (subj_mat_file_name, 'subj_data');
        end
                    
                    
        fprintf('Extracting voxels data for subject %i - %s is done.\n\n',subj_count, subj_name);  
        


       
end
     
% % Write headers to csv file
% if PARAMS.save_csv
%     fprintf('Extracting headers to csv\n\n');  
% 
%     beta_headers_list = beta_headers{1};
%     percent_signal_change_headers_list = percent_signal_change_headers{1};
%     contrast_headers_list = contrast_headers{1};
% 
%     for design_ind = 2: length(PARAMS.design_dirs)
%        beta_headers_list = [ beta_headers_list, beta_headers{design_ind}];
%        percent_signal_change_headers_list = [percent_signal_change_headers_list, percent_signal_change_headers{design_ind}];
%        contrast_headers_list = [contrast_headers_list, contrast_headers{design_ind}];
%     end
%             
%     data_header = {'"X_ind"', '"Y_ind"', '"X_ind"', '"X_coord"', '"Y_coord"', '"Z_coord"', ...
%                     beta_headers_list{:}, percent_signal_change_headers_list{:} , contrast_headers_list{:}};
% 
%     cell2csv ([PARAMS.data_dir 'headers.csv'], data_header);
% 
% end
fprintf('Done!!\n');  


       
          
