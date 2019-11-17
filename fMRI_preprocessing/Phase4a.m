function Phase4a(Settings)
% Phase4a_new   Create folders to later save ROIs 


dirs = Settings.RoiDirs;


for DesignNum=1:length(Settings.ExpDesign)
    ParentDir=[Settings.SpmDir filesep 'ROI_Analysis-' Settings.ExpDesign{DesignNum}.Name '_m'];
    
    for dirindex=1:length(dirs)
        mkdir(ParentDir,[dirs{dirindex} '_right']);
        mkdir(ParentDir,[dirs{dirindex} '_left']);
%         rmdir([ParentDir filesep 'FFA_right']);
%         rmdir([ParentDir filesep 'FFA_left']);
    end
    mkdir(ParentDir,'Excluded');
end
SpmInfo=dir(Settings.SpmDir);

for index=1:length(SpmInfo)
    if strcmp(SpmInfo(index).name,'.') || strcmp(SpmInfo(index).name,'..') || ~isempty(strfind(SpmInfo(index).name,'ROI_Analysis'))
        continue;
    end
    for DesignNum=1:size(Settings.ExpDesign,2)
        mkdir([Settings.SpmDir filesep SpmInfo(index).name filesep 'Results'  Settings.ExpDesign{DesignNum}.Name '_m' filesep  'ROI_Analysis']);
        ParentDir=[Settings.SpmDir filesep SpmInfo(index).name filesep 'Results'  Settings.ExpDesign{DesignNum}.Name '_m' filesep  'ROI_Analysis'];
        
        for dirindex=1:length(dirs)
            mkdir(ParentDir,[dirs{dirindex} '_right']);
            mkdir(ParentDir,[dirs{dirindex} '_left']);
%             rmdir([ParentDir filesep 'FFA_right']);
%         rmdir([ParentDir filesep 'FFA_left']);
        end
        mkdir(ParentDir,'Excluded');
    end
end

helpdlg('Done!','Creating ROI directories');

end