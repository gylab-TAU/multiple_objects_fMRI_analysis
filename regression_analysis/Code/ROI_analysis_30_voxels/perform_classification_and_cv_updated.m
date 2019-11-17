function [results] = perform_classification_and_cv_updated(data, settings, params)


subj_name = settings.subj_file(1:4);
% sub_folder_path = [settings.method];
results_path = [settings.path_Results];

this_time = fix(clock);
this_time_str = [num2str(this_time(3)) '_' num2str(this_time(2)) '_' num2str(this_time(1)) '_' num2str(this_time(4)) '_' num2str(this_time(5))];

results.file_name_prefix = ['class_ROI_' settings.data.data_design '_' this_time_str];

results.settings = settings;
results.params = params;

if ~exist(results_path, 'dir')
    mkdir (results_path);
end

num_of_masks = length(data.roi_masks_headers);

results.subj_name = subj_name;
results.roi_name = data.roi_masks_headers;

results.roi_size = zeros(1, num_of_masks);

results.roi_masks = data.roi_masks;
results.spm_vol = data.spm_vol;
results.conj_mask = data.conj_mask;


results.vals_psc_A_cond_train = cell(1,num_of_masks);
results.vals_psc_B_cond_train = cell(1,num_of_masks);
results.vals_psc_C_cond_train = cell(1,num_of_masks);
results.vals_psc_A_cond_test = cell(1,num_of_masks);
results.vals_psc_B_cond_test = cell(1,num_of_masks);
results.vals_psc_C_cond_test = cell(1,num_of_masks);

results.voxel_vals_psc_A_cond_train = cell(1,num_of_masks);
results.voxel_vals_psc_B_cond_train = cell(1,num_of_masks);
results.voxel_vals_psc_C_cond_train = cell(1,num_of_masks);
results.voxel_vals_psc_A_cond_test = cell(1,num_of_masks);
results.voxel_vals_psc_B_cond_test = cell(1,num_of_masks);
results.voxel_vals_psc_C_cond_test = cell(1,num_of_masks);


results.vals_psc_A_cond = cell(1,num_of_masks);
results.vals_psc_B_cond = cell(1,num_of_masks);
results.vals_psc_C_cond = cell(1,num_of_masks);
results.mean_psc_A_cond = zeros(1,num_of_masks);
results.mean_psc_B_cond = zeros(1,num_of_masks);
results.mean_psc_C_cond = zeros(1,num_of_masks);
results.vals_alpha = cell(1,num_of_masks);
results.vals_beta = cell(1,num_of_masks);
results.mean_alpha = zeros(1, num_of_masks);
results.mean_beta = zeros(1, num_of_masks);
results.mean_weights_sum = zeros(1, num_of_masks);
results.mean_weights_ratio = zeros(1, num_of_masks);
results.mean_weights_log_ratio = zeros(1, num_of_masks);
results.mean_normalized_delta_C_prediction_test = zeros(1, num_of_masks);
results.mean_normalized_delta_A_train_test = zeros(1, num_of_masks);
results.mean_normalized_delta_B_train_test = zeros(1, num_of_masks);
results.mean_normalized_delta_C_train_test = zeros(1, num_of_masks);

results.sem_psc_A_cond = zeros(1,num_of_masks);
results.sem_psc_B_cond = zeros(1,num_of_masks);
results.sem_psc_C_cond = zeros(1,num_of_masks);

results.std_alpha = zeros(1, num_of_masks);
results.std_beta = zeros(1, num_of_masks);
results.std_weights_sum = zeros(1, num_of_masks);
results.std_weights_ratio = zeros(1, num_of_masks);
results.std_weights_log_ratio = zeros(1, num_of_masks);
results.std_normalized_delta_C_prediction_test = zeros(1, num_of_masks);
results.std_normalized_delta_A_train_test = zeros(1, num_of_masks);
results.std_normalized_delta_B_train_test = zeros(1, num_of_masks);
results.std_normalized_delta_C_train_test = zeros(1, num_of_masks);


for roi_itr = 1:num_of_masks
    
    % check if there are enough voxels in this mask
    % for now we do the classificatiom based on the maximum available
    % voxels for each roi
    num_of_mask_voxels = length(find(data.roi_masks(:,roi_itr)));
    if num_of_mask_voxels < settings.min_voxel_num
        continue
    end
    
    results.roi_size(roi_itr) = num_of_mask_voxels;
    
    num_of_cv_itr = length(data.train_inds);
    
    alpha = zeros(1,num_of_cv_itr);
    beta = zeros(1,num_of_cv_itr);
    weights_sum = zeros(1,num_of_cv_itr);
    weights_ratio = zeros(1,num_of_cv_itr);
    weights_log_ratio = zeros(1,num_of_cv_itr);
    normalized_delta_C_prediction_test = zeros(1,num_of_cv_itr);
    normalized_delta_A_train_test = zeros(1,num_of_cv_itr);
    normalized_delta_B_train_test = zeros(1,num_of_cv_itr);
    normalized_delta_C_train_test = zeros(1,num_of_cv_itr);
   
    results.vals_psc_A_cond_train{roi_itr} = zeros(1,num_of_cv_itr);
    results.vals_psc_B_cond_train{roi_itr} = zeros(1,num_of_cv_itr);
    results.vals_psc_C_cond_train{roi_itr} = zeros(1,num_of_cv_itr);
    results.vals_psc_A_cond_test{roi_itr} = zeros(1,num_of_cv_itr);
    results.vals_psc_B_cond_test{roi_itr} = zeros(1,num_of_cv_itr);
    results.vals_psc_C_cond_test{roi_itr} = zeros(1,num_of_cv_itr);

    results.voxel_vals_psc_A_cond_train{roi_itr} = zeros(num_of_mask_voxels,num_of_cv_itr);
    results.voxel_vals_psc_B_cond_train{roi_itr} = zeros(num_of_mask_voxels,num_of_cv_itr);
    results.voxel_vals_psc_C_cond_train{roi_itr} = zeros(num_of_mask_voxels,num_of_cv_itr);
    results.voxel_vals_psc_A_cond_test{roi_itr} = zeros(num_of_mask_voxels,num_of_cv_itr);
    results.voxel_vals_psc_B_cond_test{roi_itr} = zeros(num_of_mask_voxels,num_of_cv_itr);
    results.voxel_vals_psc_C_cond_test{roi_itr} = zeros(num_of_mask_voxels,num_of_cv_itr);

    for cv_itr= 1:num_of_cv_itr
        
        % get A, B and C vectors for current iteration:
        % A and B vectors are mean of train data vectors with ROI mask
        A_vector_train = mean(data.A_cond_data(data.roi_masks(:,roi_itr),data.train_inds{cv_itr}),2);
        B_vector_train = mean(data.B_cond_data(data.roi_masks(:,roi_itr),data.train_inds{cv_itr}),2);
        C_vector_train = mean(data.C_cond_data(data.roi_masks(:,roi_itr),data.train_inds{cv_itr}),2);
        
        % C vector is mean of test data vectors with ROI mask
        A_vector_test = mean(data.A_cond_data(data.roi_masks(:,roi_itr),data.test_inds{cv_itr}),2);
        B_vector_test = mean(data.B_cond_data(data.roi_masks(:,roi_itr),data.test_inds{cv_itr}),2);
        C_vector_test = mean(data.C_cond_data(data.roi_masks(:,roi_itr),data.test_inds{cv_itr}),2);
        
        
        % calculate dot products for later calculations
        AA = dot(A_vector_train,A_vector_train);
        BB = dot(B_vector_train,B_vector_train);
        AB = dot(A_vector_train,B_vector_train);
        CA = dot(C_vector_train,A_vector_train);
        CB = dot(C_vector_train,B_vector_train);
        
        % calculate alpha and beta
        alpha(cv_itr) = (CB*AB - CA*BB) / (AB*AB - AA*BB);
        beta(cv_itr) = (CA - alpha(cv_itr)*AA) / AB;
        
        % calculate delta_vector
        predicted_C_vector = alpha(cv_itr)*A_vector_train + beta(cv_itr)*B_vector_train;
        
        C_test_norm = sqrt(dot(C_vector_test,C_vector_train));
        normalized_delta_C_prediction_test(cv_itr) = sqrt(dot(C_vector_train - predicted_C_vector,C_vector_train - predicted_C_vector))/C_test_norm;
        normalized_delta_A_train_test(cv_itr) = sqrt(dot(A_vector_train - A_vector_test, A_vector_train - A_vector_test))/C_test_norm;
        normalized_delta_B_train_test(cv_itr) = sqrt(dot(B_vector_train - B_vector_test, B_vector_train - B_vector_test))/C_test_norm;
        normalized_delta_C_train_test(cv_itr) = sqrt(dot(C_vector_train - C_vector_test, C_vector_train - C_vector_test))/C_test_norm;
        
        weights_sum(cv_itr) = alpha(cv_itr) + beta(cv_itr);
        weights_ratio(cv_itr) = alpha(cv_itr) ./ beta(cv_itr);
        weights_log_ratio(cv_itr) = log(weights_ratio(cv_itr));
        
        results.vals_psc_A_cond_train{roi_itr}(cv_itr) = mean(A_vector_train);
        results.vals_psc_B_cond_train{roi_itr}(cv_itr) = mean(B_vector_train);
        results.vals_psc_C_cond_train{roi_itr}(cv_itr) = mean(C_vector_train);
        results.vals_psc_A_cond_test{roi_itr}(cv_itr) = mean(A_vector_test);
        results.vals_psc_B_cond_test{roi_itr}(cv_itr) = mean(B_vector_test);
        results.vals_psc_C_cond_test{roi_itr}(cv_itr) = mean(C_vector_test);
        
        results.voxel_vals_psc_A_cond_train{roi_itr}(:,cv_itr) = A_vector_train;
        results.voxel_vals_psc_B_cond_train{roi_itr}(:,cv_itr) = B_vector_train;
        results.voxel_vals_psc_C_cond_train{roi_itr}(:,cv_itr) = C_vector_train;
        results.voxel_vals_psc_A_cond_test{roi_itr}(:,cv_itr) = A_vector_test;
        results.voxel_vals_psc_B_cond_test{roi_itr}(:,cv_itr) = B_vector_test;
        results.voxel_vals_psc_C_cond_test{roi_itr}(:,cv_itr) = C_vector_test;
        
        
        
    end
    
    results.vals_psc_A_cond{roi_itr} = mean(data.A_cond_data(data.roi_masks(:,roi_itr),:),2);
    results.vals_psc_B_cond{roi_itr} = mean(data.B_cond_data(data.roi_masks(:,roi_itr),:),2);
    results.vals_psc_C_cond{roi_itr} = mean(data.C_cond_data(data.roi_masks(:,roi_itr),:),2);
    results.mean_psc_A_cond(roi_itr) = mean(results.vals_psc_A_cond{roi_itr});
    results.mean_psc_B_cond(roi_itr) = mean(results.vals_psc_B_cond{roi_itr});
    results.mean_psc_C_cond(roi_itr) = mean(results.vals_psc_C_cond{roi_itr});
    results.sem_psc_A_cond(roi_itr) = std(results.vals_psc_A_cond{roi_itr})/sqrt(length(results.vals_psc_A_cond{roi_itr}));
    results.sem_psc_B_cond(roi_itr) = std(results.vals_psc_B_cond{roi_itr})/sqrt(length(results.vals_psc_B_cond{roi_itr}));
    results.sem_psc_C_cond(roi_itr) = std(results.vals_psc_C_cond{roi_itr})/sqrt(length(results.vals_psc_C_cond{roi_itr}));
    
    results.vals_alpha{roi_itr} = alpha;
    results.vals_beta{roi_itr} = beta;
    results.mean_alpha(roi_itr) = mean(alpha);
    results.mean_beta(roi_itr) = mean(beta);
    results.mean_weights_sum(roi_itr) = mean(weights_sum);
    results.mean_weights_ratio(roi_itr) = mean(weights_ratio);
    results.mean_weights_log_ratio(roi_itr) = mean(weights_log_ratio);
    results.mean_normalized_delta_C_prediction_test(roi_itr) = mean(normalized_delta_C_prediction_test);
    results.mean_normalized_delta_A_train_test(roi_itr) = mean(normalized_delta_A_train_test);
    results.mean_normalized_delta_B_train_test(roi_itr) = mean(normalized_delta_B_train_test);
    results.mean_normalized_delta_C_train_test(roi_itr) = mean(normalized_delta_C_train_test);
    
    results.std_alpha(roi_itr) = std(alpha);
    results.std_beta(roi_itr) = std(beta);
    results.std_weights_sum(roi_itr) = std(weights_sum);
    results.std_weights_ratio(roi_itr) = std(weights_ratio);
    results.std_weights_log_ratio(roi_itr) = std(weights_log_ratio);
    results.std_normalized_delta_C_prediction_test(roi_itr) = std(normalized_delta_C_prediction_test);
    results.std_normalized_delta_A_train_test(roi_itr) = std(normalized_delta_A_train_test);
    results.std_normalized_delta_B_train_test(roi_itr) = std(normalized_delta_B_train_test);
    results.std_normalized_delta_C_train_test(roi_itr) = std(normalized_delta_C_train_test);
    
end




save ([results_path filesep results.file_name_prefix '_' subj_name '.mat'], 'results');

% cell_header = fieldnames(results.data)';
% cell_data = squeeze(struct2cell(results.data))';
% cell_results = [cell_header; cell_data];
% cell2csv([results_path filesep results.file_name_prefix '_' subj_name '.csv'], cell_results);
% 
