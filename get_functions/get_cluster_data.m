function data = get_cluster_data(Options, n_cluster, n_samples, num_clusters)
%LOAD_CLUSTER_DATA Load all sample files for a given cluster.
%
%   data = LOAD_CLUSTER_DATA(Options, n_cluster, n_samples, num_clusters)
%
%   Inputs:
%       Options.dir.results  - Base directory containing cluster folders
%       n_cluster            - Cluster index to load (numeric)
%       n_samples            - Total number of samples
%       num_clusters         - Total number of clusters
%
%   Output:
%       data - cell array of structs (each cell is a loaded .mat file)
%
%   Example:
%       data = load_cluster_data(Options, 3, 99, 20);
    
    
    % === Compute zero-padding digits ===
    n_digits_samples  = ceil(log10(n_samples + 1));   % e.g. 2 → "01", 3 → "001"
    n_digits_clusters = ceil(log10(num_clusters + 1));

    % === Build paths ===
    save_dir = fullfile(Options.dir.results, 'clusters_for_curve_resolution');

    % --- Get list of cluster directories (exclude '.' and '..') ---
    f_list_tmp = dir(save_dir);
    f_list_tmp = f_list_tmp(3:end);
    % file_path = 'clusters.mat';
    % varName   = 'cluster_0007';

    
    % cluster_list = cluster_list([cluster_list.isdir]);
    % cluster_list = cluster_list(~ismember({cluster_list.name}, {'.', '..'}));

    % --- Find target cluster folder ---
    % cluster_str = sprintf(['%0', num2str(n_digits_clusters), 'd'], n_cluster);
    cluster_str = sprintf(['cluster_%0', num2str(n_digits_clusters), 'd'], n_cluster);
    for k = 1:numel(f_list_tmp)
        file_path = fullfile(f_list_tmp(k).folder,f_list_tmp(k).name);
        try
            % tmp = matfile('myfile.mat'); oneElem = m.s(5); 
            % tic 
            tmp = load(file_path,cluster_str);
            % toc
            % % benchmark_mat_loading(file_path, cluster_str, 5);
            data{k,1} = tmp.Z;
            data{k,2} = tmp.sample_ID;
            data{k,3} = tmp.n_samples;
            data{k,4} = tmp.n_cluster;
           
            

        catch ME
            % warning('Could not load file: %s\nReason: %s', file_path, ME.message);
            % data{k} = struct();
        end
    end

    % fprintf('✅ Loaded %d samples for cluster %s from %s\n', ...
    %         numel(f_list_tmp), cluster_str, save_dir_sample);
end

%%
% function data = get_cluster_data(Options, n_cluster, n_samples, num_clusters)
% %LOAD_CLUSTER_DATA Load all sample files for a given cluster.
% %
% %   data = LOAD_CLUSTER_DATA(Options, n_cluster, n_samples, num_clusters)
% %
% %   Inputs:
% %       Options.dir.results  - Base directory containing cluster folders
% %       n_cluster            - Cluster index to load (numeric)
% %       n_samples            - Total number of samples
% %       num_clusters         - Total number of clusters
% %
% %   Output:
% %       data - cell array of structs (each cell is a loaded .mat file)
% %
% %   Example:
% %       data = load_cluster_data(Options, 3, 99, 20);
% 
% 
%     % === Compute zero-padding digits ===
%     n_digits_samples  = ceil(log10(n_samples + 1));   % e.g. 2 → "01", 3 → "001"
%     n_digits_clusters = ceil(log10(num_clusters + 1));
% 
%     % === Build paths ===
%     save_dir = fullfile(Options.dir.results, 'clusters_for_curve_resolution');
% 
%     % --- Get list of cluster directories (exclude '.' and '..') ---
%     cluster_list = dir(save_dir);
%     % cluster_list = cluster_list([cluster_list.isdir]);
%     % cluster_list = cluster_list(~ismember({cluster_list.name}, {'.', '..'}));
% 
%     % --- Find target cluster folder ---
%     cluster_str = sprintf(['%0', num2str(n_digits_clusters), 'd'], n_cluster);
%     ind_cluster = find(contains({cluster_list.name}, cluster_str), 1);
%    if isempty(ind_cluster)
%        data = [];
%        return
%         % error('Cluster folder "%s" not found in %s', cluster_str, save_dir);
%     end
% 
%     % --- List .mat files for this cluster ---
%     save_dir_sample = fullfile(save_dir, cluster_list(ind_cluster).name);
%     f_list_tmp = dir(fullfile(save_dir_sample, '*.mat'));
% 
%     if isempty(f_list_tmp)
%         warning('No .mat files found in cluster folder: %s', save_dir_sample);
%         data = {};
%         return;
%     end
% 
%     % === Load all .mat files safely ===
%     data = cell(numel(f_list_tmp), 1);
% 
%     for k = 1:numel(f_list_tmp)
%         file_path = fullfile(f_list_tmp(k).folder, f_list_tmp(k).name);
%         try
%             % tmp = matfile('myfile.mat'); oneElem = m.s(5); 
%             tmp = load(file_path);
%             data{k,1} = tmp.Z;
%             data{k,2} = tmp.sample_ID;
%             data{k,3} = tmp.n_samples;
%             data{k,4} = tmp.n_cluster;
% 
% 
% 
%         catch ME
%             % warning('Could not load file: %s\nReason: %s', file_path, ME.message);
%             % data{k} = struct();
%         end
%     end
% 
%     % fprintf('✅ Loaded %d samples for cluster %s from %s\n', ...
%     %         numel(f_list_tmp), cluster_str, save_dir_sample);
% end

function benchmark_mat_loading(file_path, varName, n_iter)

if nargin < 3
    n_iter = 50;
end

fprintf('\n=== Benchmarking MAT-file loading ===\n');
fprintf('File: %s\n', file_path);
fprintf('Variable: %s\n\n', varName);

%% --- 1. Partial load ---
t_partial = zeros(n_iter,1);

for i = 1:n_iter
    tic
    tmp = load(file_path, varName);
    data = tmp.(varName);
    clear tmp data
    t_partial(i) = toc;
end

fprintf('Partial load (load file,var): %.4f ± %.4f sec\n', ...
    mean(t_partial), std(t_partial));

%% --- 2. Full load ---
t_full = zeros(n_iter,1);

for i = 1:n_iter
    tic
    tmp = load(file_path);
    clear tmp
    t_full(i) = toc;
end

fprintf('Full load (load file):      %.4f ± %.4f sec\n', ...
    mean(t_full), std(t_full));

%% --- 3. matfile access ---
t_matfile = zeros(n_iter,1);

try
    m = matfile(file_path);
    
    for i = 1:n_iter
        tic
        data = m.(varName);
        clear data
        t_matfile(i) = toc;
    end

    fprintf('matfile access:            %.4f ± %.4f sec\n', ...
        mean(t_matfile), std(t_matfile));
catch
    fprintf('matfile not supported (likely not v7.3)\n');
end

%% --- Summary ---
fprintf('\n=== Speed ratios ===\n');
fprintf('Full / Partial: %.2fx slower\n', mean(t_full)/mean(t_partial));

if exist('t_matfile','var') && any(t_matfile)
    fprintf('matfile / Partial: %.2fx slower\n', ...
        mean(t_matfile)/mean(t_partial));
end

fprintf('\nDone.\n\n');

end