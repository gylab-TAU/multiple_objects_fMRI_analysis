function Phase3_GLM(Settings)
% run GLM and calculate contrasts

addpath(Settings.spmpath);
spm fmri;

    
for DesignNum=2:length(Settings.ExpDesign)
     
    
          
    filePrefix = Settings.ExpDesign{DesignNum}.FilePrefix;
    
    IsContrasted = Settings.ExpDesign{DesignNum}.contrastsAcrossRuns.HasContrasts || Settings.ExpDesign{DesignNum}.contrastsSingleRun.HasContrasts; % needs to add contrasts

    for subj_itr=  1:length(Settings.Sessions)
        
        curr_subject_spm_path = [Settings.SpmDir Settings.Sessions{subj_itr}{2}];

        curr_subj_design_stats_path = [curr_subject_spm_path filesep 'Results' Settings.ExpDesign{DesignNum}.Name '_m'];
        
        if ~exist(curr_subj_design_stats_path,'dir')
            system(['mkdir ' curr_subj_design_stats_path]);
        end
        
        
        dataFolders = Settings.ExpDesign{DesignNum}.runDirs{subj_itr};
        num_of_sess = length(dataFolders);

        jobs{1}.spm.stats.fmri_spec.timing.units = 'secs';
        jobs{1}.spm.stats.fmri_spec.timing.RT = Settings.TR;
        jobs{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        jobs{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
        
        jobs{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        jobs{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        jobs{1}.spm.stats.fmri_spec.volt = 1;
        jobs{1}.spm.stats.fmri_spec.global = 'None';
        jobs{1}.spm.stats.fmri_spec.mthresh = 0.8;
        jobs{1}.spm.stats.fmri_spec.mask = {''};
        jobs{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
        jobs{1}.spm.stats.fmri_spec.dir{1} = curr_subj_design_stats_path;
        
        jobs{1}.spm.stats.fmri_spec.sess = struct(  'scans', {''},...
                                                'cond', struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {}),...
                                                'multi', cell(1,num_of_sess),...
                                                'regress', struct('name', {}, 'val', {}),...
                                                'multi_reg',{{''}},...
                                                'hpf', 128);
               

        for run_itr = 1:num_of_sess
            
            curr_run_path = [curr_subject_spm_path filesep dataFolders{run_itr}];

            d=dir([curr_run_path filesep char(filePrefix) '*.img']);
            files={d.name}';
            
            jobs{1}.spm.stats.fmri_spec.sess(1,run_itr).scans=cellstr(strcat([curr_run_path],filesep ,files,',1'));
            
            curr_design_file_path = [curr_subject_spm_path filesep 'design_files'];
            design_file = dir ([curr_design_file_path filesep Settings.ExpDesign{DesignNum}.design_file_prefix '*design_' num2str(run_itr) '*.mat']); 
            
            jobs{1}.spm.stats.fmri_spec.sess(1,run_itr).multi{1}= [curr_design_file_path filesep design_file(1).name];

%             motion_file = dir([curr_run_path filesep '*.txt']);
            motion_file = dir([curr_run_path filesep 'motion_regressors_for_GLM.mat']);
            
            jobs{1}.spm.stats.fmri_spec.sess(1,run_itr).multi_reg = {[curr_run_path filesep motion_file(1).name]};
            
            
        end
                
        % run  first level specification
        spm_jobman('run',jobs);
        
        jobs = [];
        
        SPM_file_full_path = [curr_subj_design_stats_path filesep 'SPM.mat'];
        
        jobs{1}.spm.stats.fmri_est.spmmat = {SPM_file_full_path};
        jobs{1}.spm.stats.fmri_est.write_residuals = 0;
        jobs{1}.spm.stats.fmri_est.method.Classical = 1;

        spm_jobman('run',jobs);

        jobs = [];
        
        
        %% adding contrasts
        if (IsContrasted)
            AddContrasts(SPM_file_full_path, curr_subj_design_stats_path, num_of_sess, DesignNum, Settings)
        end
                                                
            
    end
end



helpdlg(['Phase 3 complete!']);

function AddContrasts(SPM_file_full_path, curr_subj_design_stats_path,num_of_sess, DesignNum, Settings)

disp(['Adding contrasts to ' SPM_file_full_path]);
load(SPM_file_full_path);


% delete old contrasts
SPM.xCon = [];
save (SPM_file_full_path,'SPM');
delete([curr_subj_design_stats_path filesep 'spmT*.nii']);
delete([curr_subj_design_stats_path filesep 'con*.nii']);

contrastsAcrossRunsCounter = 0;
if Settings.ExpDesign{DesignNum}.contrastsAcrossRuns.HasContrasts % =1

    for constrats_id = 1:length(Settings.ExpDesign{DesignNum}.contrastsAcrossRuns.Contrasts)
        currentContrast = Settings.ExpDesign{DesignNum}.contrastsAcrossRuns.Contrasts{constrats_id};
        zeros_runsContrast = zeros(1,num_of_sess);
        
        beta_names = SPM.xX.name;
        current_full_contrast = zeros(1, length(beta_names));
        
        % bf(1) means that this regressor was convolved with a base
        % function (hrf). only condition regressors are convolved with bf
        cond_inds = find(~cellfun(@isempty,regexp(beta_names(1,:), ['*.*bf(1)*.'] )));
        current_full_contrast(cond_inds) = repmat(currentContrast,1, num_of_sess);
       
        
        currentName = Settings.ExpDesign{DesignNum}.contrastsAcrossRuns.ContrastsNames{constrats_id};
        jobs{1}.spm.stats.con.consess{constrats_id}.tcon.name = currentName;
        jobs{1}.spm.stats.con.consess{constrats_id}.tcon.weights = current_full_contrast;

    end
    
%     for constrats_id = 1:length(Settings.ExpDesign{DesignNum}.contrastsAcrossRuns.Contrasts)
%         currentContrast = Settings.ExpDesign{DesignNum}.contrastsAcrossRuns.Contrasts{constrats_id};
%         zeros_runsContrast = zeros(1,sessions_number_effective);
%         zeros_motionContrast = zeros(1,6);
%         currentContrast = [repmat([currentContrast,zeros_motionContrast],1,sessions_number_effective),zeros_runsContrast];
%         currentValue = mat2str(currentContrast);
%         currentName = Settings.ExpDesign{DesignNum}.contrastsAcrossRuns.ContrastsNames{constrats_id};
%         jobs{1}.spm.stats.con.consess{constrats_id}.tcon.name = currentName;
%         jobs{1}.spm.stats.con.consess{constrats_id}.tcon.weights = currentContrast;
% 
%     end
    
    contrastsAcrossRunsCounter = length(Settings.ExpDesign{DesignNum}.contrastsAcrossRuns.Contrasts);
end    

if Settings.ExpDesign{DesignNum}.contrastsSingleRun.HasContrasts
    N_sessions = num_of_sess;
    N_conds = length(Settings.ExpDesign{DesignNum}.ConditionsNames);
    contrasts_count = contrastsAcrossRunsCounter +1;
    templateContrast = zeros(1,(N_sessions*(N_conds+1)));
    for run_ind = 1:N_sessions
        current_session_inds = (1:N_conds) + (run_ind-1)*N_conds;
        for constrats_ind = 1:length(Settings.ExpDesign{DesignNum}.contrastsSingleRun.Contrasts)

            currentContrast = templateContrast;
            currentContrast(current_session_inds) = Settings.ExpDesign{DesignNum}.contrastsSingleRun.Contrasts{constrats_ind};
            currentValue = mat2str(currentContrast);
            if N_sessions < 10
                currentName = ['Sn_0' , num2str(run_ind), '_' Settings.ExpDesign{DesignNum}.contrastsSingleRun.ContrastsNames{constrats_ind}];
            else
                if run_ind <10
                    currentName = ['Sn_00' , num2str(run_ind), '_' Settings.ExpDesign{DesignNum}.contrastsSingleRun.ContrastsNames{constrats_ind}];
                else
                    currentName = ['Sn_0' , num2str(run_ind), '_' Settings.ExpDesign{DesignNum}.contrastsSingleRun.ContrastsNames{constrats_ind}];
                end
            end


            jobs{1}.spm.stats.con.consess{contrasts_count}.tcon.name = currentName;
            jobs{1}.spm.stats.con.consess{contrasts_count}.tcon.weights = currentContrast;

            contrasts_count = contrasts_count+1;

        end
    end
        
end    
    


jobs{1}.spm.stats.con.spmmat = {SPM_file_full_path};


spm_jobman('run',jobs);





