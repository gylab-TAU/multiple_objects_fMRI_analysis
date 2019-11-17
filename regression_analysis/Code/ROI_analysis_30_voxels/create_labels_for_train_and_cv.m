function [data] = create_labels_for_train_and_cv(settings, data)




if strcmp(settings.data.cv, 'LORO') % the only option implemented at the moment

    % the number of session is the last character of the last header:
    num_of_sess = str2num(data.data_headers{end}(end));

    sess_vec = 1:num_of_sess;
    
    data.test_inds = cell(1, num_of_sess);
    data.train_inds = cell(1, num_of_sess);
    
    for cv_itr_num = 1:num_of_sess
        
        data.test_inds{cv_itr_num} = cv_itr_num;
        data.train_inds{cv_itr_num} = sess_vec(~ismember(sess_vec,cv_itr_num));
    
    end
end

if strcmp(settings.data.cv, 'ALL_TOGETHER') % the only option implemented at the moment

    % the number of session is the last character of the last header:
    num_of_sess = str2num(data.data_headers{end}(end));
    
    data.test_inds{1} = (1:num_of_sess);
    data.train_inds{1} = (1:num_of_sess);

  
end

if strcmp(settings.data.cv, 'HALVES') % the only option implemented at the moment

    % the number of session is the last character of the last header:
    num_of_sess = str2num(data.data_headers{end}(end));
    
    sess_vec = 1:num_of_sess;
    num_of_repetition = nchoosek(num_of_sess,num_of_sess/2);
    
    all_repetition_vec = nchoosek(sess_vec, num_of_sess/2);

    data.test_inds = cell(1, num_of_repetition);
    data.train_inds = cell(1, num_of_repetition);
    
    for cv_itr_num = 1:num_of_repetition
        
        data.test_inds{cv_itr_num} = all_repetition_vec(cv_itr_num,:);
        data.train_inds{cv_itr_num} = sess_vec(~ismember(sess_vec,all_repetition_vec(cv_itr_num,:)));
    
    end

    
%     data.test_inds = cell(1, num_of_repetition/2);
%     data.train_inds = cell(1, num_of_repetition/2);
%     
%     for cv_itr_num = 1:num_of_repetition/2
%         
%         data.test_inds{cv_itr_num} = all_repetition_vec(cv_itr_num,:);
%         data.train_inds{cv_itr_num} = sess_vec(~ismember(sess_vec,all_repetition_vec(cv_itr_num,:)));
%     
%     end


  
end