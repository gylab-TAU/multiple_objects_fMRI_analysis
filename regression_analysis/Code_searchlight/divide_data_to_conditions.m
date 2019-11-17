function [data] = divide_data_to_conditions(settings, data)


% cond A
filter_data_header = ['.*_' settings.A_cond_name '.*'];
filter_var_inds = find(~cellfun(@isempty,regexp(data.data_headers, filter_data_header)));
data.A_cond_data = data.orig_data(:,filter_var_inds);
data.A_cond_headers = data.data_headers(filter_var_inds);

% cond A
filter_data_header = ['.*_' settings.B_cond_name '.*'];
filter_var_inds = find(~cellfun(@isempty,regexp(data.data_headers, filter_data_header)));
data.B_cond_data = data.orig_data(:,filter_var_inds);
data.B_cond_headers = data.data_headers(filter_var_inds);

% cond A
filter_data_header = ['.*_' settings.C_cond_name '.*'];
filter_var_inds = find(~cellfun(@isempty,regexp(data.data_headers, filter_data_header)));
data.C_cond_data = data.orig_data(:,filter_var_inds);
data.C_cond_headers = data.data_headers(filter_var_inds);


