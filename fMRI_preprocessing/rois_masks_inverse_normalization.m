function [] = rois_masks_inverse_normalization(Settings)
% convert MNI masks to subject's native space using inverse normalization
% here we use Brodman EVC mask


Settings.roi_masks_dir = './masks';


for subj_itr = 1:length(Settings.Sessions)      
    curr_subject_session_info = Settings.Sessions{subj_itr};
    curr_subject_spm_path = [Settings.SpmDir filesep curr_subject_session_info{2}];
    means_file_path = [curr_subject_spm_path filesep curr_subject_session_info{3}];
    anatomy_dir_path = [curr_subject_spm_path filesep curr_subject_session_info{end}];

    disp(['ROI MASKS INVERSE NORMALIZATION - subj: ' curr_subject_session_info{2}]);

    roi_mask_files = dir ([Settings.roi_masks_dir filesep '*.nii']);

    jobs{1}.spm.spatial.normalise.write.woptions.bb = [-200 -200 -200 ; 200 200 200];
    jobs{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1]; % this should be the voxel size of the episo that itwould be easier to comapreit to the results
    jobs{1}.spm.spatial.normalise.write.woptions.interp = 4;
    jobs{1}.spm.spatial.normalise.write.woptions.prefix = 'w';


    inv_deformation_file = dir ([anatomy_dir_path filesep 'iy_*.nii']);
    jobs{1}.spm.spatial.normalise.write.subj.def = {[anatomy_dir_path filesep inv_deformation_file(1).name]};

    subj_atlas_path = [curr_subject_spm_path filesep 'Atlas_rois'];
    resample_files = [];
    mkdir (subj_atlas_path);
    for file_itr = 1:length (roi_mask_files)
        copyfile ([Settings.roi_masks_dir filesep roi_mask_files(file_itr).name], [subj_atlas_path filesep roi_mask_files(file_itr).name]);
  
        resample_files{end+1,1}=[subj_atlas_path filesep roi_mask_files(file_itr).name ',1'];
    end

    jobs{1}.spm.spatial.normalise.write.subj.resample = resample_files;

    % run atlas inverse normalization
    spm_jobman('run' , jobs);

    jobs = [];


    % now reslice the normalized atlas mask so it would be with the same bounding box dimensions as
    means_file = dir([means_file_path filesep 'means*.img']);
    reslice_files = {[means_file_path filesep means_file(1).name ',1']};

    roi_mask_files = dir ([subj_atlas_path filesep 'w*.nii']);

    for mask_itr = 1:length(roi_mask_files)
        reslice_files{end+1,1}=[subj_atlas_path filesep roi_mask_files(mask_itr).name ',1'];
    end


    jobs{1}.spm.spatial.realign.write.data = reslice_files;
    jobs{1}.spm.spatial.realign.write.roptions.which = [1 0]; % don't reslice the first image
    jobs{1}.spm.spatial.realign.write.roptions.interp = 4;
    jobs{1}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
    jobs{1}.spm.spatial.realign.write.roptions.mask = 1;
    jobs{1}.spm.spatial.realign.write.roptions.prefix = 'r';

    % run atlas reslice
    spm_jobman('run' , jobs);

    jobs = [];

end


