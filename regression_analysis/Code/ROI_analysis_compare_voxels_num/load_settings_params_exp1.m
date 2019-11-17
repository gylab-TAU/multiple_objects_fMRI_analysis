function [ settings, params ] = load_settings_params_exp1(settings, params )
    
%% paths

settings.path_Data = '/data/libi/MRI_data/Face_Body_Integration/data_mat_files_m';
settings.path_Results = './Results';
settings.path_Figures = './Figures';


%%
settings.method = 'Euclidean'; % options: 'SVM' 'Correlationa'

settings.file_header_info = 'FaceBody_ROIs';

settings.data.data_type = 'PercSigCh'; % PercSigCh, Beta
settings.data.cv = 'ALL_TOGETHER'; % LORO = Leave-One-Run-Out cross validation
settings.data.data_design = 'model';
settings.data.sample_type = 'OVERLAP_ROIS'; %'UNITE_ROIS'; 'OVERLAP_ROIS' 'CONTRASTS', 'ROIs'
settings.data.normalize = 0; % if 1 than normalization (zscore of data) is performed




settings.min_voxel_num = 30;
settings.exact_voxel_num = 30;

settings.cond_names = {'Face' 'Body' 'Person'};
settings.A_cond_name = 'Face'; % first single, weight is alpha
settings.B_cond_name = 'Body'; % second cond, weight is beta
settings.C_cond_name = 'Person'; % combined cond


%% masks 
% if there are no masks to apply for the data, set an empty value for the masks:
% settings.masks = [];
settings.masks(1).type = 'mask_values';
settings.masks(1).header = 'Brain_activation_mask_loc';

settings.masks(2).type = 'mask_values';
settings.masks(2).header = 'Brain_activation_mask_model';


%% rois overlaps and exclusions
% relevant only if settings.data.sample_type='ROIs'
settings.rois.names = { 'FBA_left_roi'    % 1
                        'FBA_right_roi'   % 2
                        'FFA_left_roi'    % 3
                        'FFA_right_roi'   % 4
                        'EBA_left_roi'    % 5
                        'EBA_right_roi'   % 6
                        'OFA_left_roi'    % 7
                        'OFA_right_roi'}; % 8
                        
settings.roi.cont_t_vals = {'Cont_loc_Bodies>Objects'
                            'Cont_loc_Bodies>Objects'
                            'Cont_loc_Faces>Objects'
                            'Cont_loc_Faces>Objects'
                            'Cont_loc_Bodies>Objects'
                            'Cont_loc_Bodies>Objects'
                            'Cont_loc_Faces>Objects'
                            'Cont_loc_Faces>Objects'};
%                         
settings.t_vals.cont_t_vals = {'Cont_loc_Faces>Objects'
                            'Cont_loc_Bodies>Objects'
                            'Cont_loc_Persons>Objects'
                            'Cont_loc_Faces>Bodies'
                            'Cont_loc_Persons>Faces'
                            'Cont_loc_Persons>Bodies'}; 

settings.t_vals.names = {'Face > Object'
                         'Body > Object'
                         'Person > Object'
                         'Face > Body'
                         'Person > Face'
                         'Person > Body'}; 
%                         
settings.rois.pairs_for_exclusion_overlap = {[3,1] [4,2] [7,5] [8,6]};

settings.rois.exclusion_roi_headers = {'FFA_left', 'FBA_left'
                                        'FFA_right', 'FBA_right'
                                        'OFA_left', 'EBA_left'
                                        'OFA_right', 'EBA_right'};                                    

settings.rois.overlap_roi_headers = {'Fus_overlap_left'
                                    'Fus_overlap_right'
                                    'Loc_overlap_left'
                                    'Loc_overlap_right'};
                                

%% fixed params

params.seed = 1; % the seed for the randomize

