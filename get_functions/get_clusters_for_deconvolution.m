function get_clusters_for_deconvolution(matFile_list, x_range, y_range, mz_range, ...
    Filter, num_clusters, n_samples, Options)
% EXTRACT_CLUSTERS_FOR_DECONVOLUTION
% Efficiently extracts and filters sparse cluster regions from multiple samples.
%
% INPUTS:
%   matFile_list : struct array of .mat files containing 'Z'
%   x_range      : cell array of x-index ranges for each cluster
%   y_range      : cell array of y-index ranges for each cluster
%   mz_range     : cell array of m/z-index ranges for each cluster
%   Filter       : 3D convolution filter
%   num_clusters : number of clusters to extract
%   n_samples    : number of samples
%   Options      : struct with Options.dir.results path

fprintf('\n=== Extracting clusters for deconvolution ===\n');

q = progressParfor(n_samples);

save_dir = fullfile(Options.dir.results, 'clusters_for_curve_resolution');
f = dir(save_dir);
for n = 1:length(f)
delete(fullfile(f(n).folder,f(n).name));
dir(save_dir)
end 
% Zero padding
n_digits_samples  = ceil(log10(n_samples ) + 1);
n_digits_clusters = ceil(log10(num_clusters) + 1);

if ~exist(save_dir, 'dir')
    
    mkdir(save_dir);
end

%% === LOOP OVER SAMPLES ===
parfor k = 1:n_samples
    send(q, k);
    process_sample(matFile_list, k, x_range, y_range, mz_range, ...
        Filter, num_clusters, n_samples, n_digits_clusters, ...
        save_dir, n_digits_samples);
end

fprintf('\n=== All samples processed successfully ===\n');

end


function process_sample(matFile_list, k, x_range, y_range, mz_range, ...
    Filter, num_clusters, n_samples, n_digits_clusters, ...
    save_dir, n_digits_samples)
% PROCESS_SAMPLE Extract and save clusters for one sample into a single HDF5 file.

%% === LOAD SPARSE TENSOR ===
S = load(fullfile(matFile_list(k).folder, matFile_list(k).name), 'Z');
[subs, V] = find(S.Z);
I = subs(:,1); J = subs(:,2); K = subs(:,3);
clear subs S

%% === SORT BY K ===
[K, sortIdx] = sort(K);
I = I(sortIdx); J = J(sortIdx); V = V(sortIdx);

%% === PRECOMPUTE K SLICES ===
[uniqueK, ~, idxK] = unique(K, 'stable');
startK = [1; find(diff(idxK)) + 1];
endK   = [startK(2:end)-1; numel(K)];

%% === FILE PATHS ===
sample_str   = sprintf(['%0', num2str(n_digits_samples), 'd'], k);
filename = fullfile(save_dir, sprintf('sample_%s.h5', sample_str));
% temp_name    = sprintf('tmp_sample_%s_%d_%d.h5', sample_str, k, feature('getpid'));
% filename     =final_filename;

%% === OPEN HDF5 FILE ONCE ===

fid = H5F.create(filename, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
%% === FILE PATHS ===
% sample_str = sprintf(['%0', num2str(n_digits_samples), 'd'], k);
% filename   = fullfile(save_dir, sprintf('sample_%s.h5', sample_str));

% %% === OPEN HDF5 FILE ONCE ===
% if isfile(filename)
%     fid = H5F.open(filename, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
% else
%     fid = H5F.create(filename, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
% end

% H5P.close(fcpl);
% H5P.close(fapl);

%% === CLUSTER EXTRACTION ===
for n_cluster = 1:num_clusters
    xr = x_range{n_cluster};
    yr = y_range{n_cluster};
    zr = mz_range{n_cluster};

    cluster_str = sprintf(['%0', num2str(n_digits_clusters), 'd'], n_cluster);
    group_name  = sprintf('cluster_%s', cluster_str);

    %% --- Match K indices ---
    [isK, locK] = ismember(zr, uniqueK);
    locK = locK(isK);

    if isempty(locK)
        gid = H5G.create(fid, group_name, 'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');
        write_scalar_dataset(gid, 'empty', 1);
        H5G.close(gid);
        continue;
    end

    %% --- Gather indices ---
    idxRange = arrayfun(@(s,e) s:e, startK(locK), endK(locK), 'UniformOutput', false);
    idxRange = [idxRange{:}];
    Ii_all = I(idxRange); Ji_all = J(idxRange);
    Ki_all = K(idxRange); Vi_all = V(idxRange);

    %% --- Restrict to X/Y ---
    valid = ismember(Ii_all, xr) & ismember(Ji_all, yr);
    if ~any(valid)
        gid = H5G.create(fid, group_name, 'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');
        write_scalar_dataset(gid, 'empty', 1);
        H5G.close(gid);
        continue;
    end

    Ii = Ii_all(valid); Ji = Ji_all(valid);
    Ki = Ki_all(valid); Vi = Vi_all(valid);
    clear Ii_all Ji_all Ki_all Vi_all

    %% --- Map to local coordinates ---
    [~, Ii] = ismember(Ii, xr);
    [~, Ji] = ismember(Ji, yr);
    [~, Ki] = ismember(Ki, zr);

    %% --- Build dense block and convolve ---
    Zk     = accumarray([Ii, Ji, Ki], Vi, [numel(xr), numel(yr), numel(zr)], [], 0);
    Z_conv = convn(Zk, Filter, 'same');
    clear Zk

    %% --- Store as sparse COO ---
    sp          = sptensor(Z_conv);
    subs_out    = uint32(sp.subs); % integer
    vals_out    = double(sp.vals);
    tensor_size = uint32(sp.size); % integer
    clear Z_conv sp

    %% --- Write group and datasets ---
    gid = H5G.create(fid, group_name, 'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');
    write_double_dataset(gid, 'subs',        subs_out);
    write_double_dataset(gid, 'vals',        vals_out);
    write_double_dataset(gid, 'tensor_size', tensor_size);
    write_scalar_dataset(gid, 'sample_ID',   uint32(k));         % integer
    write_scalar_dataset(gid, 'n_cluster',   uint32(n_cluster)); % integer
    write_scalar_dataset(gid, 'n_samples',   uint32(n_samples)); % integer
    H5G.close(gid);

    %% --- Flush after each cluster ---
    H5F.flush(fid, 'H5F_SCOPE_GLOBAL');
end

%% === CLOSE FILE ONCE ===
H5F.close(fid);

%% === RENAME TO FINAL FILENAME WITH RETRY ===


end


function write_double_dataset(gid, name, data)
% Write a double array as a dataset into an open HDF5 group
    dims  = fliplr(size(data));
    space = H5S.create_simple(numel(dims), dims, dims);
    type  = H5T.copy('H5T_NATIVE_DOUBLE');
    dset  = H5D.create(gid, name, type, space, 'H5P_DEFAULT');
    H5D.write(dset, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', data);
    H5D.close(dset);
    H5S.close(space);
    H5T.close(type);
end


function write_scalar_dataset(gid, name, value)
% Write a scalar double as a dataset into an open HDF5 group
    space = H5S.create('H5S_SCALAR');
    type  = H5T.copy('H5T_NATIVE_DOUBLE');
    dset  = H5D.create(gid, name, type, space, 'H5P_DEFAULT');
    H5D.write(dset, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', value);
    H5D.close(dset);
    H5S.close(space);
    H5T.close(type);
end