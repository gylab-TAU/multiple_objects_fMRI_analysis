

addpath(pwd);

Settings = loadSettings();
  
    
    
%% Phase 1
% prepare data for analysis

convertDicoms(Settings);

removeDummyScans(Settings);

copy_design_files(Settings);

     
%% Phase 2
% preprocessing

Phase2c_preprocessing(Settings);

create_motion_regressors(Settings)

resample_img(Settings)

rois_masks_inverse_normalization(Settings)


%% Stop for CheckReg!!

%% phase 3
% run GLM model and calculate contrasts

Phase3_GLM(Settings)


%% phase 4
% extract ROIs
Phase4a_new(Settings)
create_marsbar_anatomy_rois(Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%
% 4b. Find ROIs with MarsBar (based on localizer data
%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [Settings] = loadSettings()
    [FileName,PathName,FilterIndex] = uigetfile('CreateSettings.m');
    
    currDir = pwd;
    cd (PathName);
    funcName = FileName(1:(end-2));
    Settings = feval(funcName);
    cd (currDir);
end

