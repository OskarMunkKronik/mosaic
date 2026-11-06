data_dir = [Options.Paths.save2mat,'/ROI/Augmented'];
% load(fullfile(data_dir,'MSroi_aug.mat'),'MSroi_aug')
load(fullfile(data_dir,'mzroi_aug.mat'),'mzroi_aug','Rt_aug','trace_aug')
mz_candidate = find(mzroi_aug >= input_parameters{13,"Options_value"} & mzroi_aug <= input_parameters{14,"Options_value"})';
Options.processing_time.convert2sparse_tensor.process.all = tic;
%% make directories BPC and TIC 
if exist(fullfile(Options.Paths.save2mat,'BPC'),'dir')  ~= 7
mkdir(fullfile(Options.Paths.save2mat,'BPC'))
end
if exist(fullfile(Options.Paths.save2mat,'TIC'),'dir')  ~= 7
mkdir(fullfile(Options.Paths.save2mat,'TIC'))
end
if exist(fullfile(Options.Paths.save2mat,'SparseTensor'),'dir') ~= 7
    mkdir(fullfile(Options.Paths.save2mat,'SparseTensor'))
end 
%% Fold data
sample_vec = 1:length(fileList_mat);
[mod,uRt,Rt2] = deal(cell(max(sample_vec),1));

for k = sample_vec
    fprintf(1,'Folding sample: %i/%i\n',k,length(sample_vec))

    Options.firstRt(k)      =   find(Rt_aug{k} >= Options.PhaseShift + Rt_aug{k}(1) ,1,'first');

    % Adjust the Rt vector
    Rt2{k} = Rt_aug{k}(Options.firstRt(k):end);
    Rt2{k} = Rt2{k}-Rt2{k}(1);

    %Find modulateoions and 2D Rts
    [mod{k},uRt{k}] = FoldChrom(Rt2{k},Options.modTime,trace_aug{k});
end 

[Options.coords.X,Options.coords.Y] = meshgrid([1:max(cellfun(@max,mod(:)))].*Options.modTime/60 +  median(cellfun(@(x)x(1),Rt_aug(sample_vec)))/60,[1:max(cellfun(@max,uRt(:)))].*median(diff(Rt_aug{k})));


%% Plot data
Options.processing_time.convert2sparse_tensor.save = tic;
loadDir = fullfile(Options.Paths.save2mat,'ROI','Augmented','Samples');
fprintf('Plotting and constructing 3D sparse tensor...\n');
% for k = sample_vec
% 
%     send(q, k);
%     load data 
%      d =     load(fullfile(loadDir, fileList_mat(k).name), 'MSroi_aug');
%        Map column indices to mod and RT coordinates
%     [ionSlice, c, val] = find(d.MSroi_aug(:,Options.firstRt(k):end));  % find nonzero entries
% 
%     modCoord = mod{k}(c);   % mod{k} is length numRt
%     rtCoord = uRt{k}(c);   % uRt{k} is length numRt
%     Plot BPC
%     p = surf(Options.coords.X,Options.coords.Y,accumarray([modCoord, rtCoord], val,[max(cellfun(@max,mod(:))), max(cellfun(@max,uRt(:)))],@max)','EdgeColor','interp','FaceColor','interp');
%     view([0 90])
%     xlabel('^1t_R (min)')
%     ylabel('^2t_R (s)')
%     title(['BPC: ',fileList(k).name(1:end-4)] ,'Interpreter', 'none')
%     axis tight
%     colorbar
%     saveas(p,fullfile(Options.Paths.save2mat,'BPC',[fileList(k).name(1:end-4),'.tif']))
% 
%     Plot TIC
%     p = surf(Options.coords.X,Options.coords.Y,accumarray([modCoord, rtCoord], val,[ max(cellfun(@max,mod(:))), max(cellfun(@max,uRt(:)))],@sum)','EdgeColor','interp','FaceColor','interp');
%         view([0 90])
%     xlabel('^1t_R (min)')
%     ylabel('^2t_R (s)')
%     title(['TIC: ', fileList(k).name(1:end-4) ],'Interpreter', 'none')
%     axis tight
%     colorbar
%     saveas(p,fullfile(Options.Paths.save2mat,'TIC',[fileList(k).name(1:end-4),'.tif']))
% 
%     Make sparse 3D tensor
%     Z = sptensor([modCoord, rtCoord, ionSlice], val,[ max(cellfun(@max,mod(:))), max(cellfun(@max,uRt(:))),length(mzroi_aug)]);
% 
%     Save data
%     save(fullfile(Options.Paths.save2mat,'SparseTensor',[fileList(k).name,'.mat']), ...
%      'Z', 'mzroi_aug', 'fileList',"Options", '-v7.3');
% end
% --- once, before the loop ---
dims = [max(cellfun(@max, mod(:))), max(cellfun(@max, uRt(:)))];
for sub = ["BPC","TIC","SparseTensor"]          % avoid mkdir races inside parfor
    if ~exist(fullfile(Options.Paths.save2mat, sub), 'dir'), mkdir(fullfile(Options.Paths.save2mat, sub)); end
end

q = progressParfor(length(sample_vec));
% afterEach(q, @(k) fprintf('Finished sample %d\n', k));

% make sure workers can find sptensor (Tensor Toolbox)
% pool = gcp; addAttachedFiles(pool, {'path/to/tensor_toolbox'});   % or addpath on workers

parfor k = sample_vec
    send(q, k);
    convert2Ztensor_and_plot(k, q, loadDir, fileList_mat, fileList, mod{k}, uRt{k}, mzroi_aug, dims, Options);
end

%% take time

Options.processing_time.convert2sparse_tensor.save = toc(Options.processing_time.convert2sparse_tensor.save);
Options.processing_time.convert2sparse_tensor.process.all = ...
    toc(Options.processing_time.convert2sparse_tensor.process.all);




%% Old


% data_dir = [Options.Paths.save2mat,'/ROI/Samples_aug'];
% % load(fullfile(data_dir,'MSroi_aug.mat'),'MSroi_aug')
% load(fullfile(data_dir,'mzroi_aug.mat'),'mzroi_aug','Rt_aug','trace_auig')
% mz_candidate = find(mzroi_aug >= input_parameters{13,"Options_value"} & mzroi_aug <= input_parameters{14,"Options_value"})';
% Options.processing_time.convert2sparse_tensor.process.all = tic;
% %%
% sample_vec = 1:length(MSroi_aug);
% [Z,mod,uRt,Rt2] = deal(cell(max(sample_vec),1));
% for k = sample_vec
%                 fprintf(1,'Folding sample: %i/%i\n',k,length(sample_vec))
% 
%     Options.firstRt(k)      =   find(Rt_aug{k} >= Options.PhaseShift + Rt_aug{k}(1) ,1,'first');
% 
%     % Adjust the Rt vector
%     Rt2{k} = Rt_aug{k}(Options.firstRt(k):end);
%     Rt2{k} = Rt2{k}-Rt2{k}(1);
% 
%     %Find modulateoions and 2D Rts
%     [mod{k},uRt{k}] = FoldChrom(Rt2{k},Options.modTime,trace_aug{k});
% end 
% 
% [Options.coords.X,Options.coords.Y] = meshgrid([1:max(cellfun(@max,mod(:)))].*Options.modTime/60 +  median(cellfun(@(x)x(1),Rt_aug))/60,[1:max(cellfun(@max,uRt(:)))].*median(diff(Rt_aug{k})));
% 
% 
% %% save BPC and TIC 
% if exist(fullfile(Options.Paths.save2mat,'BPC'),'dir')  ~= 7
% mkdir(fullfile(Options.Paths.save2mat,'BPC'))
% end
% if exist(fullfile(Options.Paths.save2mat,'TIC'),'dir')  ~= 7
% mkdir(fullfile(Options.Paths.save2mat,'TIC'))
% end
% 
% for k = sample_vec
%       fprintf(1,'Plotting TIC and BPC, sample: %i/%i\n',k,length(sample_vec))
%     % BPC
%     [r, c, val] = find(MSroi_aug{k}(:,Options.firstRt(k):end));  % find nonzero entries
%      % Map column indices to mod and RT coordinates
%     modCoord = mod{k}(c);   % mod{k} is length numRt
%     rtCoord = uRt{k}(c);   % uRt{k} is length numRt
%     p = surf(Options.coords.X,Options.coords.Y,accumarray([modCoord, rtCoord], val,[ max(cellfun(@max,mod(:))), max(cellfun(@max,uRt(:)))],@max)','EdgeColor','interp','FaceColor','interp');
%     view([0 90])
%     xlabel('^1t_R (min)')
%     ylabel('^2t_R (s)')
%     title(['BPC: ',fileList(k).name(1:end-4)] ,'Interpreter', 'none')
%     axis tight
%     colorbar
%     saveas(p,fullfile(Options.Paths.save2mat,'BPC',[fileList(k).name(1:end-4),'.tif']))
% 
%     % TIC
%     p = surf(Options.coords.X,Options.coords.Y,accumarray([modCoord, rtCoord], val,[ max(cellfun(@max,mod(:))), max(cellfun(@max,uRt(:)))],@sum)','EdgeColor','interp','FaceColor','interp');
%         view([0 90])
%     xlabel('^1t_R (min)')
%     ylabel('^2t_R (s)')
%     title(['TIC: ', fileList(k).name(1:end-4) ],'Interpreter', 'none')
%     axis tight
%     colorbar
%     saveas(p,fullfile(Options.Paths.save2mat,'TIC',[fileList(k).name(1:end-4),'.tif']))
% end 
% clear val
% clear u_ind
% clear trace_aug
% clear r
% %% 
% 
% 
% for k = sample_vec
%                 fprintf(1,'Organize to sparse tensor sample: %i/%i\n',k,length(sample_vec))
%     % k = 3
%     [r, c, val] = find(MSroi_aug{k}(:,Options.firstRt(k):end));  % find nonzero entries
%     [~,~,u_ind] = unique(c);
%     MSroi_aug{k} = [];
%     % Map column indices to mod and RT coordinates
%     modCoord = mod{k}(c);   % mod{k} is length numRt
%     rtCoord = uRt{k}(c);   % uRt{k} is length numRt
% 
%     % Use row index as 3rd mode (ion slice)
%     ionSlice = r;
% 
%     % Build sparse 3-way tensor
%     Z{k} = sptensor([modCoord, rtCoord, ionSlice], val,[ max(cellfun(@max,mod(:))), max(cellfun(@max,uRt(:))),length(mzroi_aug)]);
% 
% end
% 
% 
% %%
% Options.processing_time.convert2sparse_tensor.save = tic;
% Z_all = Z;
% if exist(fullfile(Options.Paths.save2mat,'SparseTensor'),'dir') ~= 7
%     mkdir(fullfile(Options.Paths.save2mat,'SparseTensor'))
% end 
% for k = sample_vec
%       fprintf(1,'save sample: %i/%i\n',k,length(sample_vec))
%     Z = Z_all{k};
%     save(fullfile(Options.Paths.save2mat,'SparseTensor',[fileList(k).name,'.mat']), ...
%      'Z', 'mzroi_aug', 'fileList','uRt','mod',"Options", '-v7.3');
% end 
% % save(fullfile(Options.Paths.save2mat,'Rts_2D.mat'))
% 
% Options.processing_time.convert2sparse_tensor.save = toc(Options.processing_time.convert2sparse_tensor.save);
% clear MSroi_aug
% clear Z_all
% Options.processing_time.convert2sparse_tensor.process.all = ...
%     toc(Options.processing_time.convert2sparse_tensor.process.all);