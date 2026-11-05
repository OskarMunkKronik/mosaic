function [elu, rt_mov, I0] = align_elution_profiles_all(elution_profile)
%ALIGN_ELUTION_PROFILES Aligns elution profiles by common apex.
%
%   [elu, rt_mov, I0] = ALIGN_ELUTION_PROFILES(elution_profile)
%
%   Input
%   -----
%   elution_profile : n × D numeric array
%       Each column is an elution profile.
%
%   Output
%   ------
%   elu    : m × D numeric array
%       Aligned elution profiles with zero rows removed.
%   rt_mov : 1 × D integer array
%       Row shifts applied to each column.
%   I0     : scalar integer
%       Common apex index used for alignment.

    P = elution_profile;
    [n0, D] = size(P);

    % Apex indices
    [~, ind] = max(P, [], 1);

    % Common apex index
    I0 = max(ind);

    % Shifts
    rt_mov = I0 - ind;

    % Output length
    L = n0 + max(abs(rt_mov));

    % Preallocate output
    elu = zeros(L, D, 'like', P);
    q = progressParfor(D);
    % Column-wise assignment (memory efficient)
    parfor c = 1:D
    % if mod(c, max(1, round(0.05 * D))) == 0 || c == 1 || c == D
           send(q, c);
    % end
    
    shift = rt_mov(c);

    col = zeros(L, 1, 'like', P);

    src_start = max(1, 1 - shift);
    src_end   = min(n0, L - shift);

    dst_start = src_start + shift;
    dst_end   = src_end + shift;

    col(dst_start:dst_end) = P(src_start:src_end, c);

    elu(:, c) = col;
    end

    % Remove empty rows
    elu(sum(elu, 2) == 0, :) = [];

    % Optional sanity check (comment out if not needed)
    % [~, check] = max(elu, [], 1);
    % assert(all(check == I0), 'Apex alignment failed');
end
