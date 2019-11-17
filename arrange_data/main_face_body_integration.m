
PARAMS.exp_path = '/data/libi/MRI_data/Face_Body_Integration';
PARAMS.SpmDir = [PARAMS.exp_path filesep 'spm' filesep];
PARAMS.FreesurferDir = '';

PARAMS.path_spm = '/usr/local/spm12/';
PARAMS.path_marsbar = '/usr/local/spm12/marsbar';

PARAMS.unite_rois.orig_rois_dir = ['ResultsLocalizer_m' filesep 'ROI_Analysis'];
PARAMS.unite_rois.new_rois_dir = ['ResultsLocalizer_m' filesep 'ROIs_for_prepData'];
                    

% PARAMS for prepData

PARAMS.use_fressurfer = 0;

PARAMS.design_dirs = {'ResultsLocalizer_m','ResultsModel_per_condition_new_m','ResultsLocalizer_unsmoothed_m'};
PARAMS.header_prefix = {'loc', 'model','locUnsmoothed'};
PARAMS.rois_dir = PARAMS.unite_rois.new_rois_dir;

PARAMS.data_dir = '/data/libi/MRI_data/Face_Body_Integration/data_mat_files_m/';
PARAMS.roi_anal_prefix_dir_name = 'ROI_';

PARAMS.max_subj_num = 15;
PARAMS.save_csv = false;
PARAMS.save_mat = true;


                    
PARAMS.subjects_list = {
                        'subj_01' % 1
                        'subj_02' % 2
                        'subj_03'  % 3
                        'subj_04'  % 4
                        'subj_05'  % 5
                        'subj_06'  % 6
                        'subj_07'  % 7
                        'subj_08'  % 8
                        'subj_09'  % 9
                        'subj_10' % 10
                        'subj_11' % 11
                        'subj_12' % 12
                        'subj_13' % 13
                        'subj_14' % 14
                        'subj_15' % 15
                        };
                    
                    
PARAMS.anatomy_path = { '10' %1
                        '12' %2
                        '12' %4
                        '12' %3
                        '14' %5
                        '12' %6
                        '12' %7   
                        '12' %8 
                        '12' %9 
                        '13' %10 
                        '14' %11 
                        '16' %12 
                        '12' %13 
                        '12' %14 
                        '22' %15 
                        };
 
PARAMS.all_rois = { 'EBA_left_roi'
                    'EBA_right_roi'
                    'FBA_left_roi'
                    'FBA_right_roi'
                    'p_FFA_left_roi'
                    'p_FFA_right_roi'
                    'm_FFA_left_roi'
                    'm_FFA_right_roi'
                    'OFA_left_roi'
                    'OFA_right_roi'
                    'LOC_FaceBOdy_left_roi' 
                    'LOC_FaceBOdy_right_roi'
                    'FG_FaceBody_left_roi'
                    'FG_FaceBody_right_roi'
                    'CSF_roi'
                    'EVC_roi'};
                    

%% Run code

unite_ROIs_for_prepData(PARAMS);

prep_data_wo_atlas(PARAMS);
 
fix_rois_after_prep_data(PARAMS);


