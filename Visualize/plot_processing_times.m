clf 
S = Options.processing_time;
fn = fieldnames(S);

% figure;
tiledlayout('flow')

for k = 1:numel(fn)
    d = struct2table(S.(fn{k}).process);

    nexttile
    p = bar(d{1,:});
    xticks(1:width(d))
    % p.th
    xticklabels(d.Properties.VariableNames)
    ax = gca;
    ax.TickLabelInterpreter = 'none';
    xtickangle(45)
    title(strrep(fn{k}, '_', '\_'))
end
% p = gca;y
saveas(gcf,fullfile(Options.Paths.save2mat,[ 'Options_nSamples_',num2str(n_samples_speed_test),'.tif']))