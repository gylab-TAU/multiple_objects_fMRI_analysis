function [] = unite_ROIs_for_prepData(PARAMS)



                    
                    
for subj_itr=1:length(PARAMS.subjects_list)
    
   curr_subj_dir = [PARAMS.SpmDir PARAMS.subjects_list{subj_itr} filesep];
   curr_new_roi_dir = [curr_subj_dir PARAMS.unite_rois.new_rois_dir];
   curr_orig_roi_dir = [curr_subj_dir PARAMS.unite_rois.orig_rois_dir];
%    rmdir (curr_new_roi_dir, 's');
   mkdir (curr_new_roi_dir);
   
   
   roi_dirs = dir (curr_orig_roi_dir);
   
   for roi_dir_itr = 1:length(roi_dirs)
      
       curr_roi_name = roi_dirs(roi_dir_itr).name;
       
       if roi_dirs(roi_dir_itr).isdir && ~strcmp(curr_roi_name,'.') && ~strcmp(curr_roi_name,'..')
          
           roi_list = dir ([curr_orig_roi_dir filesep curr_roi_name filesep '*.mat']);
           if isempty(roi_list)
               continue
           elseif length(roi_list)==1
               source_file = [curr_orig_roi_dir filesep curr_roi_name filesep roi_list(1).name];
               destination_file = [curr_new_roi_dir filesep curr_roi_name '_roi.mat'];
               copyfile (source_file, destination_file);
           else % more than one file per ROI, need to create a single united file

               for file_itr = 1: length(roi_list)
                   source_file = [curr_orig_roi_dir filesep curr_roi_name filesep roi_list(file_itr).name];
                   destination_file = [curr_new_roi_dir filesep curr_roi_name num2str(file_itr) '_roi.mat'];
                   copyfile (source_file, destination_file);
               end
               
              
               
               
               
           end
           
       end
   end
    
end
                    