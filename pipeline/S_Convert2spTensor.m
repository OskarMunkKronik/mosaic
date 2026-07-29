data_dir = [Options.Paths.save2mat,'/ROI/Samples_aug'];
load(fullfile(data_dir,'MSroi_aug.mat'))
load(fullfile(data_dir,'mzroi_aug.mat'))
mz_candidate = find(mzroi_aug >= input_parameters{13,"Options_value"} & mzroi_aug <= input_parameters{14,"Options_value"})';
tic;
%%
sample_vec = 1:length(MSroi_aug);
[Z,mod,uRt,Rt2] = deal(cell(max(sample_vec),1));
for k = sample_vec
                fprintf(1,'Folding sample: %i/%i\n',k,length(sample_vec))

    Options.firstRt(k)      =   find(Rt_aug{k} >= Options.PhaseShift  ,1,'first');

    % Adjust the Rt vector
    Rt2{k} = Rt_aug{k}(Options.firstRt(k):end);
    Rt2{k} = Rt2{k}-Rt2{k}(1);

    %Find modulateoions and 2D Rts
    [mod{k},uRt{k}] = FoldChrom(round(Rt2{k},2),Options.modTime,trace_aug{k});
end 

for k = sample_vec
                fprintf(1,'Organize to sparse tensor sample: %i/%i\n',k,length(sample_vec))
    % k = 3
    [r, c, val] = find(MSroi_aug{k}(:,Options.firstRt(k):end));  % find nonzero entries
    [~,~,u_ind] = unique(c);
    % Map column indices to mod and RT coordinates
    modCoord = mod{k}(c);   % mod{k} is length numRt
    rtCoord = uRt{k}(c);   % uRt{k} is length numRt

    % Use row index as 3rd mode (ion slice)
    ionSlice = r;

    % Build sparse 3-way tensor
    Z{k} = sptensor([modCoord, rtCoord, ionSlice], val,[ max(cellfun(@max,mod(:))), max(cellfun(@max,uRt(:))),length(mzroi_aug)]);
end
[Options.coords.X,Options.coords.Y] = meshgrid([1:max(cellfun(@max,mod(:)))].*Options.modTime/60,[1:max(cellfun(@max,uRt(:)))].*median(diff(Rt_aug{k})));

process_time.convert2sparse_tensor = toc;
%%
tic
Z_all = Z;
if exist(fullfile(Options.Paths.save2mat,'SparseTensor'),'dir') ~= 7
    mkdir(fullfile(Options.Paths.save2mat,'SparseTensor'))
end 
for k = sample_vec
      fprintf(1,'save sample: %i/%i\n',k,length(sample_vec))
    Z = Z_all{k};
    save(fullfile(Options.Paths.save2mat,'SparseTensor',[fileList(k).name,'.mat']), ...
     'Z', 'mzroi_aug', 'fileList','uRt','mod',"Options", '-v7.3');
end 
% save(fullfile(Options.Paths.save2mat,'Rts_2D.mat'))

process_time.convert2sparse_tensor_save = toc;
clear MSroi_aug