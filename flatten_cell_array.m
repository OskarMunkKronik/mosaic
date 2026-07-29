%% Auxillery functions
function [output] = flatten_cell_array(input)
input = vertcat(input);
out = cell2mat(cellfun(@size, input, 'UniformOutput', false));
ind = out(:,1)>1;
output = vertcat(input(~ind),vertcat(input{ind}));
output = [output{:}];
output(cellfun(@isempty,output)) = [];
end
