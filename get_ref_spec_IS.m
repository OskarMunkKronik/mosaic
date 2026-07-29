ref_spec_dir = input_parameters{2,2}{:};%'C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\SuspectScreening\ref_mass_Spectra';
if exist(ref_spec_dir,"dir") == 7
cd(ref_spec_dir)

refSpec = cell(length(mz_IS),1);
refSpec_int= cell(length(mz_IS),1);

% refSpecList = dir();
%%
for nIS = 1:length(mz_IS)
    refSpecList = dir([IS_mz.Inchikey{nIS},'\*mat' ]);
    mz_diff = mz_IS(nIS)-mz_nonIS(nIS);
    for n = 1:length(refSpecList)
        mass_spectrum_tmp = load(fullfile(refSpecList(n).folder,refSpecList(n).name));
        refSpec{nIS} = [refSpec{nIS};mass_spectrum_tmp.peaks(:,1)+mz_diff];
        refSpec_int{nIS} = [refSpec_int{nIS};mass_spectrum_tmp.peaks(:,2)];
    end
end
else 
    refSpec = [];
end 
%%
% nIS = 1;
%  mz_diff = mz_IS(nIS)-mz_nonIS(nIS);
% refSpec{nIS} = [ 199.07368,299.16235 57.07019,101.09615 ,55.05457]';
% 
% nIS = 10;
%  mz_diff = mz_IS(nIS)-mz_nonIS(nIS);
% refSpec{nIS} = [51.0231, 77.0382 ,50.0151,39.0232,65.0385]';
