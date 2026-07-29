function [Z_refolded,results,elu1,elu2,mz] = get_SEP_massSpecIdx(n_cluster,load_dir,Options,n_digits_clusters)


% Zero-padded cluster folder name
cluster_str = sprintf(['%0', num2str(n_digits_clusters), 'd'], n_cluster);
% Output filename (zero-padded and structured)
filename = sprintf("results_%s.mat", cluster_str);

if ~exist(fullfile(load_dir, filename),"file") 
    Z_refolded = [];
    results = struct(); % Initialize results structure if file does not exist
    elu1 = [];
    elu2 = [];
    mz = [];
    return
end 
% === Save ===
results = load(fullfile(load_dir, filename));
% --- Setup indices ---
% [X,Y] = meshgrid(Options.coords.X(1,results.x_range),Options.coords.Y(results.y_range,1));

% [X,Y] = meshgrid(results.x_range,results.y_range);




% --- Refold all components once per cluster ---
Z_refolded = refold_Z_blocks(results.W, results.H, ...
    results.x_range, results.y_range, results.mz_range, ...
    results.n_samples, []);

elu1 = accumarray(results.x_range', max(Z_refolded,[],2:4)'   ,[size(Options.coords.X,2) 1]);
elu2 = accumarray(results.y_range', max(Z_refolded,[],[1,3:4]),[size(Options.coords.X,1) 1]);
mz   = results.mz_range;
end

