
%% --- Visualization ---
close all
% Build reference-based adjacency for graph
A = zeros(numPoints);
for i = 1:numPoints
    for j = i+1:numPoints
        if clusters_test(i) == clusters_test(j)
            A(i,j) = 1;
            A(j,i) = 1;
        end
    end
end
G = graph(A);
figure('Name','Clustering Graph','Position',[100 100 600 500]);
p = plot(G, 'Layout','force');
title('Clustering Graph');
colors = lines(numClusters);
% for i = 1:length(clusters_test)
%     highlight(p, i, 'NodeColor', colors(clusters_test(i),:), 'MarkerSize', 5);
% end

% Elution profiles with reference overlay
figure('Name','Elution Profiles','Position',[100 100 1200 500]);
% chosenClusters = 1:min(3,numClusters);
[~,a] = min(abs([points.basePeak]'-mzTarget));
 chosenClusters = find(clusters_test(a) == clusters_test);
for k = 1:length(chosenClusters)
    c = chosenClusters(k);
    subplot(1,length(chosenClusters),k);
    hold on;
    members = find(clusters_test == c);
    ref = points(find(clusters_test==c,1)); % reference sample
    for m = members
        plot(points(m).elution1, '-', 'Color', colors(c,:)*0.7, 'HandleVisibility','off');
        plot(points(m).elution2, '--', 'Color', colors(c,:)*0.7, 'HandleVisibility','off');
    end
    % overlay reference
    plot(ref.elution1, '-k','LineWidth',2,'DisplayName','Reference e1');
    plot(ref.elution2, '--k','LineWidth',2,'DisplayName','Reference e2');
    hold off;
    title(sprintf('Cluster %d - Elution Profiles',c));
    xlabel('Time index'); ylabel('Intensity');
    legend show;
end

% Mass spectra with reference overlay
figure('Name','Mass Spectra','Position',[100 100 1200 500]);
for k = 1:length(chosenClusters)
    c = chosenClusters(k);
    subplot(1,length(chosenClusters),k);
    hold on;
    members = find(clusters_test == c);
    ref = points(find(clusters_test==c,1));
    for m = members
        stem(mzroi_aug, points(m).massSpec, 'Color', colors(c,:)*0.7, 'HandleVisibility','off','Marker','none');
    end
    stem(mzroi_aug, ref.massSpec, 'k','LineWidth',2,'DisplayName','Reference','Marker','none');
    hold off;
    title(sprintf('Cluster %d - Mass Spectra',c));
    xlabel('m/z index'); ylabel('Intensity');
    legend show;
end
