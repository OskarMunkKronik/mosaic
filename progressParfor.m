function q = progressParfor(totalIters)
%PROGRESSPARFOR Creates a progress display for parfor loops
%   Usage:
%       q = progressParfor(totalIterations);
%       parfor i = 1:totalIterations
%           % do stuff
%           send(q, i); % report progress
%       end

    q = parallel.pool.DataQueue;
    afterEach(q, @updateProgress);

    % Store progress in persistent variables
    persistent count nIter step nextThresh
    count = 0;
    nIter = totalIters;
    step = floor(nIter * 0.01); % every 5%
    nextThresh = step;

    function updateProgress(~)
        count = count + 1;
        if count >= nextThresh || count == nIter
            fprintf('Progress: %.0f%% (%d/%d)\n', ...
                round(100 * count / nIter), count, nIter);
            nextThresh = nextThresh + step;
        end
    end

    assignin('base', 'progressQueue', q); % make queue available in caller
end
