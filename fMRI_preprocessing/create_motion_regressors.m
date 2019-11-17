function []=create_motion_regressors(Settings) 
% check_motion  Calculate motion regressors and print moition summary and
% figures for each subject and each run.


Settings.FD_th = 0.9;

for subj_itr = 1:length(Settings.Sessions)
    
    curr_subject_session_info = Settings.Sessions{subj_itr};
    curr_subject_spm_path = [Settings.SpmDir filesep curr_subject_session_info{2}];
    dataFolders = curr_subject_session_info(3:(end-1)); % all EPI folders, not including the MPRANGE which is the last folder
    
    num_of_runs = length(dataFolders);
    
    rms_relative_to_first_single_dim = zeros(num_of_runs,6);
    rms_relative_to_first_all = zeros(num_of_runs,1);
    rms_sequential_diff_single_dim = zeros(num_of_runs,6);
    rms_sequential_diff_all = zeros(num_of_runs,1);
    max_abs_diff = zeros(num_of_runs,1);
    high_FD_num = zeros(num_of_runs,1);
    high_FD_percent = zeros(num_of_runs,1);
    
    for runs_itr = 1:length(dataFolders)

        currRunPath = [curr_subject_spm_path filesep dataFolders{runs_itr}];
                
        d = dir([currRunPath filesep 'rp_*.txt']);    %load file names of files for preprocessing
        
        movements = dlmread([currRunPath filesep d(1).name]);
        m = zeros(size(movements));
        m(:,1:3) = movements(:,1:3);
        % the rotation estimation is beeing transformed to mm on a 50 mm radius sphere
        m(:,4:6) = 2*pi*50 * movements(:,4:6) /360;
        rms_relative_to_first_single_dim(runs_itr,:) = sqrt(mean(m.^2));
        rms_relative_to_first_all(runs_itr) = sqrt(mean(m(:).^2));
        
        m_plus_1 = m(2:end,:);
        diff = [zeros(1,6);m_plus_1 - m(1:end-1,:)];
        abs_diff = abs(diff);
        FD = sum(abs_diff,2);
                
        high_FD_inds = find(FD>Settings.FD_th);
        high_FD_vol_num = length(high_FD_inds);
        scrubbing_regressors = zeros(length(FD),high_FD_vol_num);
        for FD_itr = 1:high_FD_vol_num
            scrubbing_regressors(high_FD_inds(FD_itr),FD_itr)=1;
        end
        
        m_sq = m.^2;
        m_sq_plus_1 = m_sq(2:end,:);
        diff_m_sq = [zeros(1,6);m_sq_plus_1 - m_sq(1:end-1,:)];
        
        R = [m, diff, m_sq, diff_m_sq, scrubbing_regressors];
        
        save ([currRunPath filesep 'motion_regressors_for_GLM.mat'] ,'R');
        
        high_FD_num(runs_itr) = high_FD_vol_num;
        high_FD_percent(runs_itr) = high_FD_vol_num/length(FD);
        
        rms_sequential_diff_single_dim(runs_itr,:) = sqrt(mean(diff.^2));
        rms_sequential_diff_all(runs_itr) = sqrt(mean(diff(:).^2));
        max_abs_diff(runs_itr) = max(diff(:));
        
        figure
        subplot(2,1,1)
        plot (m)
        legend({'x', 'y', 'z', 'pitch', 'roll', 'yaw'}, 'Location', 'northwest', 'NumColumns', 2); 
        title('Displacement relative to first volume')
        
        subplot(2,1,2)
        plot (diff)
        legend({'x', 'y', 'z', 'pitch', 'roll', 'yaw'}, 'Location', 'northwest', 'NumColumns', 2); 
        title('Displacement relative previous volume')
        saveas(gcf, [curr_subject_spm_path filesep 'motion_estmators_run_' dataFolders{runs_itr} '.jpg']);
        close(gcf)
    end
    
    save ([curr_subject_spm_path filesep 'motion_estimators.mat'] ,...
            'rms_relative_to_first_single_dim',...
            'rms_relative_to_first_all',...
            'rms_sequential_diff_single_dim',...
            'rms_sequential_diff_all',...
            'max_abs_diff',...
            'high_FD_num',...
            'high_FD_percent'...
            );
        
end