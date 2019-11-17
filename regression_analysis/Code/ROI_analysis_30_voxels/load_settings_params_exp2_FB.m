function [ settings, params ] = load_settings_params_exp2_FB(settings, params )
    
%% paths

settings.path_Data = '/data/libi/MRI_data/Face_Object_Integration/data_mat_files_m/';
settings.path_Results = './Results/exp2_FB';
settings.path_Figures = './Figures/exp2_FB';



%%
settings.method = 'Euclidean'; % options: 'SVM' 'Correlationa'

settings.file_header_info = 'FB';

settings.data.data_type = 'PercSigCh'; % PercSigCh, Beta
settings.data.cv = 'ALL_TOGETHER'; % LORO = Leave-One-Run-Out cross validation
settings.data.data_design = 'FB_model';
settings.data.sample_type = 'ROIs_all_else_excluded'; %'UNITE_ROIS'; 'OVERLAP_ROIS' 'CONTRASTS', 'ROIs'
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
settings.masks(2).header = 'Brain_activation_mask_FB_model';


%% rois overlaps and exclusions
% relevant only if settings.data.sample_type='ROIs'

settings.rois.names = { 'FFA_left_roi'    % 1
                        'FFA_right_roi'   % 2
                        'FBA_left_roi'    % 3
                        'FBA_right_roi'   % 4
                        'medial_obj_left_roi'   % 5
                        'medial_obj_right_roi'  % 6
                        'OFA_left_roi'    % 7
                        'OFA_right_roi'   % 8
                        'EBA_left_roi'    % 9
                        'EBA_right_roi'   % 10
                        'LO_left_roi'     % 11
                        'LO_right_roi'    % 12
                        }; % 


settings.roi.cont_t_vals = {'Cont_loc_Faces>all_others'
                            'Cont_loc_Faces>all_others'
                            'Cont_loc_Bodies>all_others'
                            'Cont_loc_Bodies>all_others'
                            'Cont_loc_Objects>all_others'
                            'Cont_loc_Objects>all_others'
                            'Cont_loc_Faces>all_others'
                            'Cont_loc_Faces>all_others'
                            'Cont_loc_Bodies>all_others'
                            'Cont_loc_Bodies>all_others'
                            'Cont_loc_Objects>all_others'
                            'Cont_loc_Objects>all_others'
                            };

                        
settings.rois.pairs_for_exclusion_overlap = {[1,3] [2,4] [1,5] [2,6] [1,7] [2,8] [9,11] [10,12] [9,13], [10,14]};

settings.rois.exclusion_roi_headers = { 'FFA_left_fba', 'FBA_left'
                                        'FFA_right_fba', 'FBA_right'
                                        'FFA_left_pfs', 'pFs_left'
                                        'FFA_right_pfs', 'pFs_right'
                                        'FFA_left_medObj','medObj_left'
                                        'FFA_right_medObj','medObj_right'
                                        'OFA_left_fba', 'EBA_left'
                                        'OFA_right_fba', 'EBA_right'
                                        'OFA_left_lo', 'LO_left'
                                        'OFA_right_lo', 'LO_right'};                                    

settings.rois.overlap_roi_headers = {   'Overlap_ffa_fba_left'
                                        'Overlap_ffa_fba_right'
                                        'Overlap_ffa_pfs_left'
                                        'Overlap_ffa_pfs_right'
                                        'Overlap_ffa_medObj_left'
                                        'Overlap_ffa_medObj_right'
                                        'Overlap_ofa_eba_left'
                                        'Overlap_ofa_eba_right'
                                        'Overlap_ofa_lo_left'
                                        'Overlap_ofa_lo_right'
                                    };

                                
settings.t_vals.cont_t_vals = { 'Cont_loc_Faces>all_others'
                                'Cont_loc_Bodies>all_others'
                                'Cont_loc_Objects>all_others'
                                }; 

settings.t_vals.names = {'Face > All others'
                         'Body > All others'
                         'Object > All others'
                         };                                 
                                

%% fixed params

params.seed = 1; % the seed for the randomize

