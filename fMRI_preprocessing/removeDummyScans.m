function removeDummyScans(Settings)
% removeDummyScans      removes dummy scan from further analysis.
% The first volumes from each run (specified with Settings.NumOfDummyScans) are moved
% from the run folder and moved to a dummyscan folder

    currPath=pwd;

    numOfSubjects = length(Settings.Sessions);


    for subjInd=1:numOfSubjects
        % original and spm folder name
        currSubjectData = Settings.Sessions{subjInd};
        currSubjSpmName = currSubjectData{2};    
        this_subject_root_spm_path = [Settings.SpmDir filesep currSubjSpmName];
        
        dummy_scans_dir = [this_subject_root_spm_path filesep 'dummyscans'];

        if ~exist(dummy_scans_dir,'dir'),
            system(['mkdir ' dummy_scans_dir]);            
        end

        d = dir(this_subject_root_spm_path);

        numOfRuns = length(currSubjectData) - 2; % the first 2 cells are for the original folder and spm folder names
        for runs_ind=1:numOfRuns

            session_id_str = currSubjectData{runs_ind+2};


                session_images_number = length(dir([this_subject_root_spm_path,filesep,...
                                               session_id_str,filesep,'*.img']));

                if (session_images_number<10) % MPRAGE
                    continue;
                end    

                if (session_images_number >= 100)
                    file_pattern = 's0*_00';
                else
                    file_pattern = 's0*_0';
                end

                dummy_scans_removed_counter = 0;
                for dummy_scan_id=1:Settings.NumOfDummyScans
                    current_dummy_scan_str = int2str(dummy_scan_id);

                    mrDir = dir([this_subject_root_spm_path,filesep,session_id_str,filesep,file_pattern current_dummy_scan_str '.*']);

                    for j=1:length(mrDir)
                        movefile([this_subject_root_spm_path,'/',session_id_str,'/',mrDir(j).name],...
                                 [this_subject_root_spm_path,'/dummyscans/',mrDir(j).name])

                        dummy_scans_removed_counter = dummy_scans_removed_counter + 1;
                    end
                end

                if (Settings.NumOfDummyScans * 2 ~= dummy_scans_removed_counter) % there are two files for every TR to be removed: .img and .hdr
                    warning(['For sessionid: ' d(runs_ind).name ' only '...
                            num2str(dummy_scans_removed_counter) ' files were removed ('...
                            'instead of ' num2str(number_of_dummy_scans * 2) ')']); 
                end
%             end
        cd(currPath);  

        end
    end
end