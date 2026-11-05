function MSroi = remove_baseline(MSroi,n_smoothpoints)
% y =MSroi{k,T}(n,:);
for n = 1:size(MSroi,2)
Y = full(MSroi);

baseline = smoothdata(Y, 2, 'movmedian', n_smoothpoints);

Y = Y - baseline;
Y(Y < 0) = 0;

MSroi = sparse(Y);
end