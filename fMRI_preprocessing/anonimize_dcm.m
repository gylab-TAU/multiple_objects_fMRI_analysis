function anonimize_dcm(root_dir)

folders = dir(root_dir);

for folder_itr=3:length(folders)
    if folders(folder_itr).isdir
        
        files = dir ([root_dir filesep folders(folder_itr).name filesep '*.dcm']);
        
        if ~isempty(files)
           
            for file_itr = 1:length(files)
               
                dicomanon ([root_dir filesep folders(folder_itr).name filesep files(file_itr).name]);
                
            end
            
        end
   
    end
end
