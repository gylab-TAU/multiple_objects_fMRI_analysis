function [] = create_marsbar_anatomy_rois(Settings)
% create_mask_rois  create marsbar ROI files that are based on anatomy
% (Brodman EVC, CSF data)


Settings.roi_masks_dir = '/data/libi/fMRI_analyses/gyfmri_libi_12/masks';


for subj_itr = 1:length(Settings.Sessions)      
    
    curr_subject_session_info = Settings.Sessions{subj_itr};
    curr_subject_spm_path = [Settings.SpmDir filesep curr_subject_session_info{2}];
%     means_file_path = [curr_subject_spm_path filesep curr_subject_session_info{3}];
    anatomy_dir_path = [curr_subject_spm_path filesep curr_subject_session_info{end}];
    subj_atlas_path = [curr_subject_spm_path filesep 'Atlas_rois'];
    ROI_Analysis_path = [curr_subject_spm_path filesep '/ResultsLocalizer_m/ROI_Analysis/'];
%     rmdir([ROI_Analysis_path ']'],'s');
    
    % EVC roi
    curr_roi_dir = [ROI_Analysis_path filesep 'EVC'];
    mkdir (curr_roi_dir);
    img_file = [subj_atlas_path filesep 'rwBrodmann_EVC_mask.nii'];
    o = maroi_image(struct('vol', spm_vol(img_file), 'binarize',0, 'func', 'img'));
    saveroi(o, [curr_roi_dir filesep 'EVC_roi.mat']);
    
    % CSF roi
    curr_roi_dir = [ROI_Analysis_path filesep 'CSF'];
    mkdir (curr_roi_dir);
    img_file = dir ([anatomy_dir_path filesep 'rc3*.nii']);
    img_file = [anatomy_dir_path filesep img_file(1).name];
    o = maroi_image(struct('vol', spm_vol(img_file), 'binarize',0, 'func', 'img>0.99'));
    saveroi(o, [curr_roi_dir filesep 'CSF_th_0.99_roi.mat']);
end
