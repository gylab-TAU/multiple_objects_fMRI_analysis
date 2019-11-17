function [data] = calculate_roi_masks(data, settings)

switch settings.data.sample_type
    
    case 'AAL_ROIS'
        
        data.roi_masks(~data.conj_mask,:) = 0;
        
        
    case 'ROIs'                              
 
%         num_of_voxels = size(data.roi_masks,1);        
        data.orig_roi_masks = data.roi_masks;
        data.roi_masks = zeros(size(data.orig_roi_masks));
        
        
        for roi_itr = 1:length(data.roi_masks_headers)
            if length(find(data.orig_roi_masks(:,roi_itr)))>= settings.exact_voxel_num % if there are less voxels than needed then the new mask is zeros
%                 curr_roi_t_val = data.roi_t_vals(:,roi_itr);
%                 curr_roi_t_val(~data.orig_roi_masks(:,roi_itr)) = -1000;
%                 [temp,sorted_inds] = sort(curr_roi_t_val,'descend');
%                 data.roi_masks(sorted_inds(1:settings.exact_voxel_num),roi_itr) = 1;

                  data.roi_masks(:,roi_itr) =  data.orig_roi_masks(:,roi_itr); 
            end
        end
        
    case 'ROIs_no_voxel_exact_num'
        
        
%         face_mask_wo_excl = data.roi_t_val_mask(:,1) > settings.roi.t_val_treshold;
%         body_mask_wo_excl = data.roi_t_val_mask(:,2) > settings.roi.t_val_treshold;
%         
%         overlap_mask = face_mask_wo_excl & body_mask_wo_excl;
%         face_mask = face_mask_wo_excl & ~body_mask_wo_excl;
%         body_mask = ~face_mask_wo_excl & body_mask_wo_excl;
%         none_mask = ~face_mask_wo_excl & ~body_mask_wo_excl;
%         
%         for orig_roi_itr = 1:length(data.orig_roi_masks)
%            
%             data.roi_mask(:,1 + (1-orig_roi_itr)) = data.orig_roi_masks & face_mask;
%             data.roi_mask(:,2 + (1-orig_roi_itr)) = data.orig_roi_masks & body_mask;
%             data.roi_mask(:,3 + (1-orig_roi_itr)) = data.orig_roi_masks & overlap_mask;
%             data.roi_mask(:,4 + (1-orig_roi_itr)) = data.orig_roi_masks & none_mask;
%             
%             data.roi_masks_headers{1 + (1-orig_roi_itr)} = 
%             
%         end
%         settings.rois.names = {
%                 'FG_person_left_roi'
%                 'FG_person_right_roi'
%                 'LOC_person_left_roi'
%                 'LOC_person_right_roi'};
%                     
% 
%                         
% settings.roi.t_val_for_areas = {'Cont_loc_Faces>Objects'
%                             'Cont_loc_Bodies>Objects'};
%                         
% settings.roi.t_names_for_areas = {'Face','Body'};
% settings.roi.t_val_treshold = 3.748312; % p = 0.0001

        
        
    case '3_WAY_OVERLAP_ROIS'
        num_of_new_rois = 7 * length(settings.rois.trios_for_exclusion_overlap);
        
        data.orig_roi_masks = data.roi_masks;
        data.roi_masks = zeros(size(data.orig_roi_masks,1), num_of_new_rois);
        data.roi_masks_headers = cell(1, num_of_new_rois);
        
        for trio_itr = 1:length(settings.rois.trios_for_exclusion_overlap)
            trio = settings.rois.trios_for_exclusion_overlap{trio_itr};
            orig_mask_1 = data.orig_roi_masks(:,trio(1));
            orig_mask_2 = data.orig_roi_masks(:,trio(2));
            orig_mask_3 = data.orig_roi_masks(:,trio(3));
            data.roi_masks(:,1 + 7*(trio_itr-1)) = orig_mask_1 & ~orig_mask_2 & ~orig_mask_3; % only 1
            data.roi_masks(:,2 + 7*(trio_itr-1)) = ~orig_mask_1 & orig_mask_2 & ~orig_mask_3; % only 2
            data.roi_masks(:,3 + 7*(trio_itr-1)) = ~orig_mask_1 & ~orig_mask_2 & orig_mask_3; % only 3
            data.roi_masks(:,4 + 7*(trio_itr-1)) = orig_mask_1 & orig_mask_2 & ~orig_mask_3; % 1 & 2
            data.roi_masks(:,5 + 7*(trio_itr-1)) = orig_mask_1 & ~orig_mask_2 & orig_mask_3; % 1 & 3
            data.roi_masks(:,6 + 7*(trio_itr-1)) = ~orig_mask_1 & orig_mask_2 & orig_mask_3; % 2 & 3
            data.roi_masks(:,7 + 7*(trio_itr-1)) = orig_mask_1 & orig_mask_2 & orig_mask_3; % 1 & 2 & 3
            
            data.roi_masks_headers((1:7) + 7*(trio_itr-1)) = settings.roi.trio_headers{trio_itr};
        end
        
        
    case 'OVERLAP_ROIS'
        num_of_new_rois = 3 * length(settings.rois.pairs_for_exclusion_overlap); % for every pair we have 2 exclusion rois and 1 overlap roi
        
        data.orig_roi_masks = data.roi_masks;
        data.roi_masks = zeros(size(data.orig_roi_masks,1), num_of_new_rois);
        data.roi_masks_headers = cell(1, num_of_new_rois);

        for roi_pairs_itr = 1:length(settings.rois.pairs_for_exclusion_overlap)

            curr_pair = settings.rois.pairs_for_exclusion_overlap{roi_pairs_itr};

            roi_mask_1 = data.orig_roi_masks(:,curr_pair(1));
            roi_mask_2 = data.orig_roi_masks(:,curr_pair(2));

            overlap_mask = roi_mask_1 & roi_mask_2;
            exclusion_mask_1 = roi_mask_1 & ~overlap_mask;
            exclusion_mask_2 = roi_mask_2 & ~overlap_mask;

            overlap_mask(~data.conj_mask) = 0;
            exclusion_mask_1(~data.conj_mask) = 0;
            exclusion_mask_2(~data.conj_mask) = 0;
            
            
%             
%             new_exclusion_mask_1 = zeros(size(exclusion_mask_1));
%             new_exclusion_mask_2 = zeros(size(exclusion_mask_2));
%             new_overlap_mask = zeros(size(overlap_mask));
%             
%             if length(find(exclusion_mask_1))>= settings.exact_voxel_num % if there are less voxels than needed then the new mask is zeros
%                 curr_roi_t_val = data.roi_t_vals(:,curr_pair(1));
%                 curr_roi_t_val(~exclusion_mask_1) = -1000;
%                 [temp,sorted_inds] = sort(curr_roi_t_val,'descend');
%                 new_exclusion_mask_1(sorted_inds(1:settings.exact_voxel_num)) = 1;
%             end
%             
%             if length(find(exclusion_mask_2))>= settings.exact_voxel_num % if there are less voxels than needed then the new mask is zeros
%                 curr_roi_t_val = data.roi_t_vals(:,curr_pair(2));
%                 curr_roi_t_val(~exclusion_mask_2) = -1000;
%                 [temp,sorted_inds] = sort(curr_roi_t_val,'descend');
%                 new_exclusion_mask_2(sorted_inds(1:settings.exact_voxel_num)) = 1;
%             end
%             
%             if length(find(overlap_mask))>= settings.exact_voxel_num % if there are less voxels than needed then the new mask is zeros
% %                 roi_t_val_1 = data.roi_t_vals(:,curr_pair(1));
% %                 roi_t_val_2 = data.roi_t_vals(:,curr_pair(2));
%                 roi_t_vals = [data.roi_t_vals(:,curr_pair(1)), data.roi_t_vals(:,curr_pair(2))];
%                 max_t_vals = repmat(max(roi_t_vals(overlap_mask,:)), size(roi_t_vals,1),1);
%                 
%                 dist_vec = sum((max_t_vals-roi_t_vals).^2,2).^0.5;
%                 dist_vec(~overlap_mask) = 1000;
%                 
%                 [temp,sorted_inds] = sort(dist_vec,'ascend');
%                 new_overlap_mask(sorted_inds(1:settings.exact_voxel_num)) = 1;
%             end
%             
            data.roi_masks(:,(1:3)+3*(roi_pairs_itr-1)) = [exclusion_mask_1, exclusion_mask_2, overlap_mask];
            data.roi_masks_headers((1:3)+3*(roi_pairs_itr-1)) = {settings.rois.exclusion_roi_headers{roi_pairs_itr,:}, settings.rois.overlap_roi_headers{roi_pairs_itr}};

        end

        
    case 'UNITE_ROIS'
        
        num_of_new_rois =length(settings.rois.union_roi_headers); 
        
        data.orig_roi_masks = data.roi_masks;
        data.roi_masks = zeros(size(data.orig_roi_masks,1), num_of_new_rois);
        data.roi_masks_headers = settings.rois.union_roi_headers;

        for roi_comb_itr = 1:length(settings.rois.orig_rois_for_unite)

            curr_comb = settings.rois.orig_rois_for_unite{roi_comb_itr};            

            for orig_roi_itr = 1:length(curr_comb)
               
                data.roi_masks(:,roi_comb_itr) = data.roi_masks(:,roi_comb_itr) | (data.orig_roi_masks(:,curr_comb(orig_roi_itr)));
            end  
            
            data.roi_masks(~data.conj_mask,roi_comb_itr) = 0;
        end
        
        
        
    case 'CONTRASTS'
        
        num_of_voxels = size(data.contrasts_inclusion_masks,1);
        num_of_new_rois = length(settings.contrasts_masks.exclusion_names); 

        data.roi_masks = ones(num_of_voxels, num_of_new_rois);
        data.roi_masks_headers = settings.contrasts_masks.exclusion_names;

        for roi_ind = 1:num_of_new_rois

            curr_contrasts = settings.contrasts_masks.new_masks{roi_ind};

            inclusion_mask = ones(num_of_voxels, 1);
            exclusion_mask = zeros(num_of_voxels, 1);

            % first we go through all +1 contrasts and calculate their
            % conjunction
            for cont_itr = 1:length(curr_contrasts)
                if curr_contrasts(cont_itr) == 1
                    inclusion_mask = inclusion_mask & data.contrasts_inclusion_masks(:,cont_itr);
                end
            end

            % then we go through all -1 contrasts and union for exclusion
            for cont_itr = 1:length(curr_contrasts)
                if curr_contrasts(cont_itr) == -1
                    exclusion_mask = exclusion_mask | data.contrasts_exclusion_masks(:,cont_itr);
                end
            end

            data.roi_masks(:,roi_ind) = inclusion_mask & ~exclusion_mask;
        end          
       
end

data.roi_masks = logical(data.roi_masks);  