function Z_cube = build_4D_cluster_cube(Z_data, n_cluster, n_samples)
%BUILD_4D_CLUSTER_CUBE Build a 4D cube (X×Y×MZ×sample) from cell or struct input.
%
%   Z_cube = BUILD_4D_CLUSTER_CUBE(Z_data, n_cluster, n_samples)
%
%   Inputs:
%       Z_data    - Cell array {n_samples,2} with {Z_block, sample_ID}
%                   OR struct array with fields: .Z and .sample_ID
%       n_cluster - Cluster index (for labeling/debugging)
%       n_samples - Total number of samples (4th dim size)
%
%   Output:
%       Z_cube    - Numeric 4D cube [X × Y × MZ × n_samples]
%
%   Notes:
%       • Each Z block is inserted at position sample_ID in the 4th dimension.
%       • Missing samples are left as zeros.
%       • Automatically handles both cell or struct inputs.

    % --- Extract Z blocks and sample IDs ---
    if iscell(Z_data)
        % Case: {Z_block, sample_ID}
        Z_cells = Z_data(:,1);
        sample_ids = cell2mat(Z_data(:,2));
    elseif isstruct(Z_data)
        % Case: struct array with fields .Z and .sample_ID
        Z_cells = {Z_data.Z}.';
        sample_ids = [Z_data.sample_ID].';
    else
        error('Z_data must be a cell array or struct array.');
    end

    % --- Convert all to double arrays ---
    Z_cells = cellfun(@double, Z_cells, 'UniformOutput', false);

    % --- Remove empty entries ---
    valid_mask = ~cellfun(@isempty, Z_cells);
    Z_cells = Z_cells(valid_mask);
    sample_ids = sample_ids(valid_mask);

    % --- Check that something is left ---
    if isempty(Z_cells)
        error('No valid Z blocks found for cluster %d.', n_cluster);
    end

    % --- Get cube size from first nonempty block ---
    [nx, ny, nz] = size(Z_cells{1});

    % --- Preallocate full 4D cube ---
    Z_cube = zeros(nx, ny, nz, n_samples, 'like', Z_cells{1});

    % --- Assign each Z block to its correct sample position ---
    cellfun(@(Zb, sid) assign_to_cube(Zb, sid), Z_cells, num2cell(sample_ids));

    % --- Nested helper for assigning data efficiently ---
    function assign_to_cube(Zb, sid)
        if sid >= 1 && sid <= n_samples
            Z_cube(:,:,:,sid) = Zb;
        else
            warning('Sample ID %d out of range [1, %d]. Skipped.', sid, n_samples);
        end
    end

end
