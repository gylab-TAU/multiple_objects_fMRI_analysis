
function [] = resample_img(Settings)
% resample segmentation files so that they would be in the same resolution
% as the EPI data. (used to later crate maskes from CSF data)



for subj_itr = 1:length(Settings.Sessions)
    
    curr_subj_data = Settings.Sessions{subj_itr};
    subj_name = curr_subj_data{2};
    
    subj_spm_path = [Settings.SpmDir filesep subj_name];
    
    mean_file_path = [subj_spm_path filesep curr_subj_data{3}];
    files = dir ([mean_file_path filesep 'means*.img']);
    mean_file_full_path = [mean_file_path filesep files(1).name];
    V_mean_file = spm_vol(mean_file_full_path);
           
    anatomy_path = [subj_spm_path filesep curr_subj_data{end}];
    segmentation_files = dir ([anatomy_path filesep 'c*.nii']);
    
    for file_itr = 1:length(segmentation_files)
       
        V = spm_vol([anatomy_path filesep segmentation_files(file_itr).name]);
        
        VV = [V_mean_file,V];
        spm_reslice(VV,struct('mean',false,'which',1,'interp',4,'prefix','r')); % 1 for linear
        
    end
    
end