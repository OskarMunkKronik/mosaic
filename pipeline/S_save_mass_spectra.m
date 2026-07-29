%% Save_mass spectra
mass_spectra_cell = {feature_groups_all.massSpec};

% Remove zeros and map mz
valid_idx = cellfun(@(s) s > 0, mass_spectra_cell, 'UniformOutput', false);

mz_in_mass_spectra = cellfun(@(mz,idx) mz(idx), ...
    repmat({mzroi_aug}, size(valid_idx)), valid_idx, 'UniformOutput', false);

mass_spectra_cell = cellfun(@(spec,idx) spec(idx), ...
    mass_spectra_cell, valid_idx, 'UniformOutput', false);

len_spectra = cellfun(@numel, mz_in_mass_spectra);

save(fullfile(Options.Paths.save2mat,'mass_spectra.mat'),"len_spectra","mass_spectra_cell","mz_in_mass_spectra")

%% Save MSP file 
saveMSP(fullfile(Options.Paths.save2mat,"mass_spectra.msp"), feature_groups_all, mzroi_aug, feature_table,3,Options)
