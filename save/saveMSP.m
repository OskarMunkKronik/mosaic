function saveMSP(filename_out, feature_groups_all, mzroi_aug, feature_table, num_fragments, Options)
% saveMSP - Save mass spectra to MSP format file for import into MZmine 3
%
% INPUTS:
%   filename_out       - Output .msp file path (string)
%   feature_groups_all - Struct array with .massSpec field
%   mzroi_aug          - m/z values vector
%   feature_table      - Matrix where col 3 = RT1 index, col 4 = RT2 index,
%                        col 5 = sample count, col 6 = MS1 height
%   num_fragments      - Minimum number of fragment peaks required
%   Options            - Struct with Options.coords.X and Options.coords.Y

%% ================== MASS SPECTRA PREP ==================
[f_sort,idx] = sortrows(feature_table,[2,6],'descend');
[~,ia] = unique(f_sort(:,2));

mass_spectra_cell = {feature_groups_all(idx(ia)).massSpec};

% Remove zero-intensity entries and map to mz axis
valid_idx = cellfun(@(s) s > 0, mass_spectra_cell, 'UniformOutput', false);
mz_in_mass_spectra = cellfun(@(mz, idx) mz(idx), ...
    repmat({mzroi_aug}, size(valid_idx)), valid_idx, 'UniformOutput', false);
mass_spectra_cell = cellfun(@(spec, idx) spec(idx), ...
    mass_spectra_cell, valid_idx, 'UniformOutput', false);


% % Normalize each spectrum to base peak = 999 (MSP convention)
% mass_spectra_norm = cellfun(@(spec) spec / max(spec) * 999, ...
%     mass_spectra_cell, 'UniformOutput', false);

% Filter by minimum fragment count
len_spectra = cellfun(@numel, mz_in_mass_spectra);
vec = find(len_spectra >= num_fragments);

% Coordinate grids
rt2_dim  = size(Options.coords.X, 1);
rt1_time = Options.coords.X(1, :);   % RT1 in minutes
rt2_time = Options.coords.Y(:, 1);   % RT2 in seconds
n_samples = max(feature_table(:, 5));

%% ================== WRITE MSP FILE ==================
% Force .msp extension
% [fpath, fname, ~] = fileparts(filename_out);
% filename_out = fullfile(fpath, [fname '.msp']);

fid = fopen(filename_out, 'w');
if fid == -1
    error('Could not open file: %s', filename_out);
end

n_features = length(vec);

for i_count = 1:n_features
    i = vec(i_count);

    mz_vals  = full(mz_in_mass_spectra{i});
    int_vals = full(mass_spectra_cell{i});   % use normalised intensities
    n_peaks  = numel(mz_vals);

    if n_peaks == 0
        continue
    end

    % RT indices
    rt1 = feature_table(i, 3);
    rt2 = feature_table(i, 4);

    scan_id = (rt1 - 1) * rt2_dim + rt2;

    % Precursor m/z = m/z of most intense peak
    [~, base_idx] = max(int_vals);
    precursor_mz  = mz_vals(base_idx);

    rt_sec = rt1_time(rt1) * 60;   % convert RT1 from minutes to seconds

    % ------------------------------------------------------------------
    % MSP block — MZmine 3 recognises these field names:
    %   Name, PrecursorMZ, RT, ADDUCT, Comment, Num Peaks
    % Fields are "Key: Value" (colon + space).
    % The peak list follows immediately after "Num Peaks:" with one
    % "mz\tintensity" pair per line.  A blank line ends the record.
    % ------------------------------------------------------------------
    fprintf(fid, 'Name: row%d_mz%.4f_rt%.2f\n',   i, precursor_mz, rt_sec);
    fprintf(fid, 'PrecursorMZ: %.4f\n',            precursor_mz);
    fprintf(fid, 'RT: %.4f\n',                     rt_sec);          % seconds
    fprintf(fid, 'RT2: %.4f\n',                    rt2_time(rt2));
    fprintf(fid, 'ADDUCT: [M+]+\n');
    fprintf(fid, 'MSLEVEL: 2\n');
    fprintf(fid, 'IONMODE: Positive\n');
    fprintf(fid, 'CHARGE: 1\n');
    fprintf(fid, 'SCANS: %d\n',                    scan_id);
    fprintf(fid, 'FEATURE_ID: %d\n',               i);
    fprintf(fid, 'FEATURE_MS1_HEIGHT: %d\n',       feature_table(i, 6));
    fprintf(fid, 'MERGED_ACROSS_N_SAMPLES: %d\n',  n_samples);
    fprintf(fid, 'SPECTYPE: CORRELATED_MS\n');
    fprintf(fid, 'SOURCE: Mosaic_pipeline\n');
    fprintf(fid, 'Comment: FEATURE_ID=%d SCANS=%d MS1_HEIGHT=%d\n', ...
        i, scan_id, feature_table(i, 6));

    % Peak list — MUST be last, immediately preceded by "Num Peaks:"
    fprintf(fid, 'Num Peaks: %d\n', n_peaks);
    for j = 1:n_peaks
        fprintf(fid, '%.5f\t%.1f\n', mz_vals(j), int_vals(j));
    end

    fprintf(fid, '\n');   % blank line separates records
end

fclose(fid);
fprintf('Saved %d spectra to: %s\n', n_features, filename_out);
end

% function saveMSP(filename_out, feature_groups_all, mzroi_aug, feature_table,num_fragments,Options)
% % saveMSP - Save mass spectra to MSP format file
% %
% % INPUTS:
% %   filename_out       - Output .msp file path (string)
% %   feature_groups_all - Struct array with .massSpec field
% %   mzroi_aug          - m/z values vector
% %   feature_table      - Table/matrix where column 3 = RT start, column 4 = RT end (seconds)
% 
% %% ================== MASS SPECTRA PREP ==================
% mass_spectra_cell = {feature_groups_all.massSpec};
% 
% % Remove zeros and map mz
% valid_idx = cellfun(@(s) s > 0, mass_spectra_cell, 'UniformOutput', false);
% mz_in_mass_spectra = cellfun(@(mz,idx) mz(idx), ...
%     repmat({mzroi_aug}, size(valid_idx)), valid_idx, 'UniformOutput', false);
% mass_spectra_cell = cellfun(@(spec,idx) spec(idx), ...
%     mass_spectra_cell, valid_idx, 'UniformOutput', false);
% 
% % Normalize each spectrum to base peak = 100
% mass_spectra_norm = cellfun(@(spec) spec / max(spec) * 100, ...
%     mass_spectra_cell, 'UniformOutput', false);
% 
% len_spectra = cellfun(@numel, mz_in_mass_spectra);
% vec = find(len_spectra >= num_fragments) ;
% 
% rt2_dim = size(Options.coords.X,1);
% rt1_time = Options.coords.X(1,:);
% rt2_time = Options.coords.Y(:,1);
% n_samples = max(feature_table(:,5));
% %% ================== WRITE MGF FILE ==================
% 
% fid = fopen(filename_out, 'w');
% if fid == -1
%     error('Could not open file: %s', filename_out);
% end
% 
% n_features = length(vec);%numel(feature_groups_all);
% 
% for i_count = 1:n_features
%     i = vec(i_count);
%     mz_vals  = full(mz_in_mass_spectra{i});
%     int_vals = full(mass_spectra_cell{i});
%     n_peaks  = numel(mz_vals);
% 
%     if n_peaks == 0 
%         continue  % Skip empty spectra
%     end
% 
%     % if feature_table(i,6) < 10000 
%     %     continue  % Skip empty spectra
%     % end
%     % RT start and end from columns 3 and 4
%     rt1 = feature_table(i, 3);
%     rt2 = feature_table(i, 4);
% 
%     % Representative m/z: m/z of the most intense peak
%     [~, base_idx] = max(int_vals);
%     pepmass = mz_vals(base_idx);
% 
%     fprintf(fid, 'BEGIN IONS\n');
%     fprintf(fid, 'FEATURE_ID=%d\n',       i);
%     fprintf(fid, 'FEATURE_FULL_ID=row%d_mz%.4f_rt%d_id\n', i, pepmass, (rt1-1) * rt2_dim + rt2);
%     fprintf(fid, 'FEATURELIST_FEATURE_ID=Mosaic_pipeline\n');
%     fprintf(fid, 'MSLEVEL=2\n');
%     fprintf(fid, 'RTINSECONDS=%.2f\n',    rt1_time(rt1)*60);
%     fprintf(fid, 'RT2_RTINSECONDS=%.3f\n',  rt2_time(rt2));
%     fprintf(fid, 'ADDUCT=[M+]\n');
%     fprintf(fid, 'PEPMASS=%.4f\n',        pepmass);
%     fprintf(fid, 'CHARGE=1\n');
%     fprintf(fid, 'FEATURE_MS1_HEIGHT=%d\n',       feature_table(i,6));
%     fprintf(fid, 'MERGED_ACROSS_N_SAMPLES=%d\n',       n_samples);
%     fprintf(fid, 'SPECTYPE=''CORRELATED MS''\n');
%     % fprintf(fid, 'FILENAME=\n');
%     fprintf(fid, 'SCANS=%d\n',(rt1-1) * rt2_dim + rt2);
%     fprintf(fid, 'Num peaks=%d\n',        n_peaks);
% 
%     for j = 1:n_peaks
%         fprintf(fid, '%.5f %.1f\n', mz_vals(j), int_vals(j));
%     end
% 
%     fprintf(fid, 'END IONS\n\n');
% end
% 
% fclose(fid);
% fprintf('Saved %d spectra to: %s\n', n_features, filename_out);
% end
% 
