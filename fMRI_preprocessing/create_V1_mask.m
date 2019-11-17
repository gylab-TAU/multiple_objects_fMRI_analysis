function [] = create_aal_masks(Settings)

brodmann_original_file = 'brodmann.nii';
V = spm_vol(brodmann_original_file);

V_data = spm_read_vols(V);

% create_masks
mask_numbers = 17;
            
mask_headers = {'EVC'};

V_mask = V;
V_mask.dt = [2 0];

for mask_itr = 1:length(mask_numbers)
    
    logical_mask_data = [];
    logical_mask_data = V_data == mask_numbers(mask_itr);
    V_mask.fname = ['Brodmann_' mask_headers{mask_itr} '_mask.nii'];
    
    spm_write_vol (V_mask, logical_mask_data);
    
end