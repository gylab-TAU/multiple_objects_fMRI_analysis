function [times] = convertDicoms(Settings)
% convertDicoms     convert dcm files to nii files and arrange folders.
% original data is at ./data/subjects and converted to nii (img/hdr - spm 8 format) and saved to
% ./data/spm .
% Data from each run is saved to a folder titled with a number. the number
% matches the number of data folder from the original raw data (dcm files
% from  the scanner)
% This function uses dcm2nii.exe (https://www.nitrc.org/projects/dcm2nii/) and a dcm2nii.ini file.

    %% dcm2nii

    try
        CheckGzip=importdata(Settings.dcm2nii_ini_path); %Check proper configuartion of mricron file.
        if (~strcmp(CheckGzip{4},'4D=0') || ~strcmp(CheckGzip{7},'Gzip=0'))
            errordlg('Your mricron settings file is misconfigured. You should go to D:\tools\mricron\dcm2nii.ini and edit 4D and Gzip lines to equal zero (4D=0, Gzip=0). Please note that the file may be hidden.')
            return;
        else
            disp('MRIcron settings file OK');
        end
    catch
        errordlg ('There was an error checking MRIcron configuartion file. Proceed at your own risk');
    end
    
    %%
    currPath=pwd;

    numOfSubjects = length(Settings.Sessions);
    times = zeros(1,numOfSubjects);

    for subjInd=1:numOfSubjects
        tic
        % original and spm folder name
        currSubjectData = Settings.Sessions{subjInd};
        currSubjDirName = currSubjectData{1};
        currSubjSpmName = currSubjectData{2};

        disp(['Starting subject: ' currSubjDirName]);

        disp('Starting dcm2nii conversion...');            
    %     dataFolders={sessions{subjInd,3:end}};
        numOfRuns = length(currSubjectData) - 2; % the first 2 cells are for the original folder and spm folder names
        subjectsOriginalFolders = dir ([Settings.SubjDir currSubjDirName]);
        subjectsOriginalFoldersNames = cell(1,length(subjectsOriginalFolders));
        [subjectsOriginalFoldersNames{:}] = deal(subjectsOriginalFolders.name);

        for runs_ind=1:numOfRuns
            session_id_str = currSubjectData{runs_ind+2};
            if strcmp(session_id_str,'')
                continue
            end
            originalFolderInd = find(strncmp(session_id_str, subjectsOriginalFoldersNames,2));
            fullCurrRunDirName = subjectsOriginalFoldersNames{originalFolderInd};

            session_dir = [Settings.SubjDir currSubjDirName filesep fullCurrRunDirName];


            disp(['Subject: ' currSubjSpmName ' Session: ' session_id_str]);      


            if ~exist(session_dir,'dir'),
                uiwait(errordlg(['Directory ' session_dir ' doesn''t exist. Error in ' sessionsMatFile '. Can not proceed']));          
            end

            this_subject_root_spm_path = [Settings.SpmDir filesep currSubjSpmName];

            if ~exist(this_subject_root_spm_path,'dir'),
               mkdir (this_subject_root_spm_path);   
            end



            convert_path = [this_subject_root_spm_path filesep session_id_str];            

            if ~exist(convert_path,'dir'),        
                mkdir (convert_path);
            else
                warning(['Session spm path: ' convert_path ' already exists']);
            end


            %Pick up this subjects' info
            waitfor(cd(session_dir));
            DirsList=dir('*.dcm');
            if length(DirsList)==0
                % try to see if the dir has IMA files (Seimens)
                DirsList=dir('*.IMA');
            end
            if length(DirsList)==0
                error(['No files found in ' currSubjDirName filesep 'MR' session_id_str]);
            end
            DicomInfo=dicominfo(DirsList(1).name);

            try % if it's a GE machine
                TR=DicomInfo.RepetitionTime/1000;
                SLICES=DicomInfo.ImagesInAcquisition;
%                 SeriesDescription=FileList.SeriesDescription;
  
            catch % for seimens
                DicomInfo=dicominfo(DirsList(end).name);
                TR=DicomInfo.RepetitionTime/1000;
               %  SLICES=DicomInfo.EchoTrainLength;
               SLICES = [];
%                 SeriesDescription=FileList.SeriesDescription;
            end
            waitfor(cd([this_subject_root_spm_path filesep session_id_str]));
            save 'Info.mat' 'TR' 'SLICES';


            exec_path = [Settings.dcm2nii_path ' -b ' Settings.dcm2nii_ini_path ' -o ' convert_path ' ' session_dir];
            system(exec_path); 
        end
    %     waitbar(1,Waiting,['Removing ' int2str(number_of_dummy_scans) ' scans at each session...']);
    times(subjInd) = toc;
    end
    cd(currPath);
end
