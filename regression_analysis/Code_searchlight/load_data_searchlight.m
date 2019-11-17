function [data] = load_data_searchlight(settings)

fprintf('%s: Load data from file %s\n', settings.subj_file(1:4), settings.subj_file);


% data = [];


%% load subj data file
load ([settings.path_Data filesep settings.subj_file]);

data.subj_name = settings.subj_file(1:4);
num_voxels = size(subj_data.data,1);
data.spm_vol = subj_data.spm_vol;



%% calculate masks

masks = cell(1,length(settings.masks));

data.conj_mask = logical(ones(num_voxels,1));

for mask_ind=1:length(settings.masks)

    switch settings.masks(mask_ind).type

        case 'mask_values'
            var_ind = find(strcmp(settings.masks(mask_ind).header, subj_data.header));
            masks{mask_ind} = subj_data.data(:,var_ind);
        
        case 'above_threshold'
            var_ind = find(strcmp(settings.masks(mask_ind).header, subj_data.header));
            masks{mask_ind} = subj_data.data(:,var_ind) > settings.masks(mask_ind).thrshold;
        
        case 'parcellation_values'
           
           var_ind = find(strcmp([settings.masks(mask_ind).header '_' settings.masks(mask_ind).hemisphere], subj_data.header));
            if strcmp(settings.masks(mask_ind).header, 'aparc_a2009')
                [~,~,label_ind] = intersect(settings.masks(mask_ind).labels, subj_data.parcellation.aparc_a2009_labels);
                masks{mask_ind} = subj_data.data(:,var_ind) == subj_data.parcellation.aparc_a2009_numbers(label_ind);
            elseif strcmp(settings.masks(mask_ind).header, 'aparc')
%                 label_ind = find(strcmp(settings.masks(mask_ind).labels, subj_data.parcellation.aparc_labels));
                
                [~,~,label_ind] = intersect(settings.masks(mask_ind).labels, subj_data.parcellation.aparc_labels);
                masks{mask_ind} = zeros(size(subj_data.data(:,var_ind)));
                for label_itr = 1:length(label_ind)
                    masks{mask_ind} = masks{mask_ind} | (subj_data.data(:,var_ind) == subj_data.parcellation.aparc_numbers(label_ind(label_itr)));
                end
            end
            
            
        case 'union'            
            curr_mask = settings.masks(mask_ind);
            submasks = cell(1,length(curr_mask.submask));

            masks{mask_ind} = zeros(num_voxels,1);
            
            for submask_ind = 1:length(curr_mask.submask)
                
                switch curr_mask.submask(submask_ind).type

                    case 'mask_values'
                        var_ind = find(strcmp(curr_mask.submask(submask_ind).header, subj_data.header));
                        submasks{submask_ind} = subj_data.data(:,var_ind);

                    case 'above_threshold'
                        var_ind = find(strcmp(curr_mask.submask(submask_ind).header, subj_data.header));
                        submasks{submask_ind} = subj_data.data(:,var_ind) > curr_mask.submask(submask_ind).thrshold;
                
                    case 'parcellation_values'
                        var_ind = find(strcmp([curr_mask.submask(submask_ind).header '_' curr_mask.submask(submask_ind).hemisphere], subj_data.header));
                        if strcmp(curr_mask.submask(submask_ind).header, 'aparc_a2009')
                            label_ind = find(strcmp(curr_mask.submask(submask_ind).labels, subj_data.parcellation.aparc_a2009_labels));
                            submasks{submask_ind} = subj_data.data(:,var_ind) == subj_data.parcellation.aparc_a2009_numbers(label_ind);
                        elseif strcmp(curr_mask.submask(submask_ind).header, 'aparc')
                            label_ind = find(strcmp(curr_mask.submask(submask_ind).labels, subj_data.parcellation.aparc_labels));
                            submasks{submask_ind} = subj_data.data(:,var_ind) == subj_data.parcellation.aparc_numbers(label_ind);
                        end
                end
               
                masks{mask_ind} = masks{mask_ind} | submasks{submask_ind};
            end
        
%             data.union_submasks = submasks;
%             test_data.union_submasks = submasks;
    end
    
    data.conj_mask = data.conj_mask & logical(masks{mask_ind});
end

%% save inds
var_ind = find(strcmp('voxel_ind', subj_data.header));
data.inds.voxel_ind = subj_data.data(:,var_ind);
% data.ind.voxel_inds = data.ind.voxel_inds(data.conj_mask,:);

var_ind = find(strcmp('X_ind', subj_data.header));
data.inds.X_ind = subj_data.data(:,var_ind);
% data.ind.X_ind = data.ind.X_ind(data.conj_mask,:);

var_ind = find(strcmp('Y_ind', subj_data.header));
data.inds.Y_ind = subj_data.data(:,var_ind);
% data.ind.Y_ind = data.ind.Y_ind(data.conj_mask,:);

var_ind = find(strcmp('Z_ind', subj_data.header));
data.inds.Z_ind = subj_data.data(:,var_ind);
% data.ind.Z_ind = data.ind.Z_ind(data.conj_mask,:);


var_ind = find(strcmp('X_coord', subj_data.header));
data.coords.X = subj_data.data(:,var_ind);
% data.ind.X_ind = data.ind.X_ind(data.conj_mask,:);

var_ind = find(strcmp('Y_coord', subj_data.header));
data.coords.Y = subj_data.data(:,var_ind);
% data.ind.Y_ind = data.ind.Y_ind(data.conj_mask,:);

var_ind = find(strcmp('Z_coord', subj_data.header));
data.coords.Z = subj_data.data(:,var_ind);
% data.ind.Z_ind = data.ind.Z_ind(data.conj_mask,:);



%% filter variables - model

% settings.data.data_type = 'Beta'; % PercSigCh, Beta
% settings.data.cv = 'LORO'; % LORO = Leave-One-Run-Out cross validation
% settings.data.data_design = {'model'};
% settings.data.normalize = 0; % if 1 than normalization (zscore of data) is performed

filter_data_header = [settings.data.data_type '_' settings.data.data_design '.*Sn.*' ];
filter_var_inds = find(~cellfun(@isempty,regexp(subj_data.header, filter_data_header)));


data.orig_data = subj_data.data(:,filter_var_inds);
% data.orig_data = data.orig_data(data.conj_mask,:);
data.data_headers = subj_data.header(filter_var_inds);

%% filter variables - localizer t-values

data.loc_t_vals = zeros(size(data.orig_data,1),length(settings.t_vals.cont_t_vals));
data.loc_t_names = settings.t_vals.names;


for cont_itr=1:length(settings.t_vals.cont_t_vals)
    
    filter_t_val_header = settings.t_vals.cont_t_vals{cont_itr};
    filter_t_val_ind = find(~cellfun(@isempty,regexp(subj_data.header, ['^' filter_t_val_header '$'])));

    data.loc_t_vals(:,cont_itr) = subj_data.data(:,filter_t_val_ind);
%     data.loc_t_vals(:,cont_itr) = data.loc_t_vals(:,cont_itr)(data.conj_mask,:);
end





%% ROI masks

switch settings.data.sample_type
    
    case 'AAL_ROIS'
        
        filter_var_inds = find(~cellfun(@isempty,regexp(subj_data.header, 'ROI_mask*')));

        data.roi_masks = logical(subj_data.data(:,filter_var_inds));
        data.roi_masks_headers = subj_data.header(filter_var_inds);
        
    case 'ROIs'
        data.roi_masks = zeros(size(data.orig_data,1),length(settings.rois.names));
        data.roi_masks_headers = cell(1,length(settings.rois.names));
        data.roi_t_vals = zeros(size(data.orig_data,1),length(settings.rois.names));
        data.roi_t_vals_headers = cell(1,length(settings.rois.names));

        for roi_itr = 1: length(settings.rois.names)

            filter_roi_mask_header = settings.rois.names{roi_itr};
            filter_roi_mask_ind = find(~cellfun(@isempty,regexp(subj_data.header, ['^' filter_roi_mask_header '$'])));

            data.roi_masks(:,roi_itr) = subj_data.data(:,filter_roi_mask_ind);
%             data.roi_masks(:,roi_itr) = roi_mask(data.conj_mask,:);
            data.roi_masks_headers(roi_itr) = subj_data.header(filter_roi_mask_ind);

            filter_roi_t_val_header = settings.roi.cont_t_vals{roi_itr};
            filter_roi_t_val_ind = find(~cellfun(@isempty,regexp(subj_data.header, ['^' filter_roi_t_val_header '$'])));

            data.roi_t_val = subj_data.data(:,filter_roi_t_val_ind);
%             data.roi_t_vals(:,roi_itr) = t_val(data.conj_mask,:);
            data.roi_t_vals_headers(roi_itr) = subj_data.header(filter_roi_t_val_ind);

        end
        
    case 'ROIs_no_voxel_exact_num'
        data.roi_masks = zeros(size(data.orig_data,1),length(settings.rois.names));
        data.roi_masks_headers = cell(1,length(settings.rois.names));
        data.roi_t_vals = zeros(size(data.orig_data,1),length(settings.rois.names));
        data.roi_t_vals_headers = cell(1,length(settings.rois.names));

        for roi_itr = 1: length(settings.rois.names)

            filter_roi_mask_header = settings.rois.names{roi_itr};
            filter_roi_mask_ind = find(~cellfun(@isempty,regexp(subj_data.header, ['^' filter_roi_mask_header '$'])));

            data.roi_masks(:,roi_itr) = subj_data.data(:,filter_roi_mask_ind);
%             data.roi_masks(:,roi_itr) = roi_mask(data.conj_mask,:);
            data.roi_masks_headers(roi_itr) = subj_data.header(filter_roi_mask_ind);

            filter_roi_t_val_header = settings.roi.cont_t_vals{roi_itr};
            filter_roi_t_val_ind = find(~cellfun(@isempty,regexp(subj_data.header, filter_roi_t_val_header)));

            data.roi_t_val = subj_data.data(:,filter_roi_t_val_ind);
%             data.roi_t_vals(:,roi_itr) = t_val(data.conj_mask,:);
            data.roi_t_vals_headers(roi_itr) = subj_data.header(filter_roi_t_val_ind);

        end
    case 'OVERLAP_ROIS'
        
        data.roi_masks = zeros(size(data.orig_data,1),length(settings.rois.names));
        data.roi_masks_headers = cell(1,length(settings.rois.names));
        data.roi_t_vals = zeros(size(data.orig_data,1),length(settings.rois.names));
        data.roi_t_vals_headers = cell(1,length(settings.rois.names));

        for roi_itr = 1: length(settings.rois.names)

            filter_roi_mask_header = settings.rois.names{roi_itr};
            filter_roi_mask_ind = find(~cellfun(@isempty,regexp(subj_data.header, ['^' filter_roi_mask_header '$'])));

            data.roi_masks(:,roi_itr) = subj_data.data(:,filter_roi_mask_ind);
%             data.roi_masks(:,roi_itr) = roi_mask(data.conj_mask,:);
            data.roi_masks_headers(roi_itr) = subj_data.header(filter_roi_mask_ind);

            filter_roi_t_val_header = settings.roi.cont_t_vals{roi_itr};
            filter_roi_t_val_ind = find(~cellfun(@isempty,regexp(subj_data.header, filter_roi_t_val_header)));

            data.roi_t_vals(:,roi_itr) = subj_data.data(:,filter_roi_t_val_ind);
%             data.roi_t_vals(:,roi_itr) = t_val(data.conj_mask,:);
            data.roi_t_vals_headers(roi_itr) = subj_data.header(filter_roi_t_val_ind);

        end 
        
        
    case '3_WAY_OVERLAP_ROIS'
        data.roi_masks = zeros(size(data.orig_data,1),length(settings.rois.names));
        data.roi_masks_headers = cell(1,length(settings.rois.names));
        data.roi_t_vals = zeros(size(data.orig_data,1),length(settings.rois.names));
        data.roi_t_vals_headers = cell(1,length(settings.rois.names));

        for roi_itr = 1: length(settings.rois.names)

            filter_roi_mask_header = settings.rois.names{roi_itr};
            filter_roi_mask_ind = find(~cellfun(@isempty,regexp(subj_data.header, ['^' filter_roi_mask_header '$'])));

            data.roi_masks(:,roi_itr) = subj_data.data(:,filter_roi_mask_ind);
%             data.roi_masks(:,roi_itr) = roi_mask(data.conj_mask,:);
            data.roi_masks_headers(roi_itr) = subj_data.header(filter_roi_mask_ind);

            filter_roi_t_val_header = settings.roi.cont_t_vals{roi_itr};
            filter_roi_t_val_ind = find(~cellfun(@isempty,regexp(subj_data.header, filter_roi_t_val_header)));

            data.roi_t_vals(:,roi_itr) = subj_data.data(:,filter_roi_t_val_ind);
%             data.roi_t_vals(:,roi_itr) = t_val(data.conj_mask,:);
            data.roi_t_vals_headers(roi_itr) = subj_data.header(filter_roi_t_val_ind);

        end     
    case 'UNITE_ROIS'
        data.roi_masks = zeros(size(data.orig_data,1),length(settings.rois.names));
        data.roi_masks_headers = cell(1,length(settings.rois.names));
        data.roi_t_vals = zeros(size(data.orig_data,1),length(settings.rois.names));
        data.roi_t_vals_headers = cell(1,length(settings.rois.names));

        for roi_itr = 1: length(settings.rois.names)

            filter_roi_mask_header = settings.rois.names{roi_itr};
            filter_roi_mask_ind = find(~cellfun(@isempty,regexp(subj_data.header, ['^' filter_roi_mask_header '$'])));

            data.roi_masks(:,roi_itr) = subj_data.data(:,filter_roi_mask_ind);
%             data.roi_masks(:,roi_itr) = roi_mask(data.conj_mask,:);
            data.roi_masks_headers(roi_itr) = subj_data.header(filter_roi_mask_ind);

            filter_roi_t_val_header = settings.roi.cont_t_vals{roi_itr};
            filter_roi_t_val_ind = find(~cellfun(@isempty,regexp(subj_data.header, ['^' filter_roi_t_val_header '$'])));

            data.roi_t_vals(:,roi_itr) = subj_data.data(:,filter_roi_t_val_ind);
%             data.roi_t_vals(:,roi_itr) = t_val(data.conj_mask,:);
            data.roi_t_vals_headers(roi_itr) = subj_data.header(filter_roi_t_val_ind);

        end 
        
        
    case 'CONTRASTS'
        
        data.contrasts_t_vals = zeros(size(data.orig_data,1),length(settings.contrasts_masks.original_masks));
        data.contrasts_t_vals_headers = cell(1,length(settings.contrasts_masks.original_masks));
        data.contrasts_inclusion_masks = zeros(size(data.orig_data,1),length(settings.contrasts_masks.original_masks));
        data.contrasts_exclusion_masks = zeros(size(data.orig_data,1),length(settings.contrasts_masks.original_masks));
        
        for cont_itr = 1: length(settings.contrasts_masks.original_masks)

            filter_cont_t_val_header = settings.contrasts_masks.original_masks{cont_itr};
            filter_cont_t_val_ind = find(~cellfun(@isempty,regexp(subj_data.header, [filter_cont_t_val_header '$'])));

            t_val = subj_data.data(:,filter_cont_t_val_ind);
            data.contrasts_t_vals(:,cont_itr) = t_val(data.conj_mask,:);
            data.contrasts_inclusion_masks(:,cont_itr) = data.contrasts_t_vals(:,cont_itr) > settings.contrasts_masks.inclusion_t_threshold;
            data.contrasts_exclusion_masks(:,cont_itr) = data.contrasts_t_vals(:,cont_itr) > settings.contrasts_masks.exclusion_t_threshold;
            data.contrasts_t_vals_headers(cont_itr) = subj_data.header(filter_cont_t_val_ind);
        end
        
    case 'NONE'
end



