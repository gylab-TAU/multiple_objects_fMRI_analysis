function copy_design_files(Settings)
% copy_design_files     copy the design files from the rad data (subjects
% folder) to the analyzed data fokder (spm)

    currPath=pwd;

    numOfSubjects = length(Settings.Sessions);


    for subjInd=1:numOfSubjects
        % original and spm folder name
        currSubjectData = Settings.Sessions{subjInd};
        
        currSubjSpmName = currSubjectData{2};    
        this_subject_root_spm_path = [Settings.SpmDir filesep currSubjSpmName];
        currSubjRawPath = [Settings.SubjDir filesep currSubjectData{1}];
        
        orig_folder = [currSubjRawPath filesep 'design_files'];
        dest_folder = [this_subject_root_spm_path filesep 'design_files'];
        copyfile( orig_folder, dest_folder);
        
    end
end