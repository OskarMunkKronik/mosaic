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

    % Column-wise assignment (memory efficient)
    for c = 1:D
        shift = rt_mov(c);

        src_start = max(1, 1 - shift);
        src_end   = min(n0, L - shift);

        dst_start = src_start + shift;
        dst_end   = src_end + shift;

        elu(dst_start:dst_end, c) = P(src_start:src_end, c);
    end

    % Remove empty rows
    elu(sum(elu, 2) == 0, :) = [];

    % Optional sanity check (comment out if not needed)
    % [~, check] = max(elu, [], 1);
    % assert(all(check == I0), 'Apex alignment failed');
end
