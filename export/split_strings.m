function [strings] = split_strings(suspect_results,unique_idx,column,n_digits)

        % This returns a cell array of strings
        vals = suspect_results(unique_idx, column);
        if isnumeric(vals{:,1})
            strings = sprintf(['%.',num2str(n_digits),'f\n'], vals{:,1});
        else 
            % vals = strings
        vals{:,1} = regexp(vals{:,1}, '^[^/]*', 'match', 'once');
        % Convert to a single comma-separated string
        strings = strjoin(vals{:,1}, '\n');
        end 
end 