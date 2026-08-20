% 

%% Load Options file 
ROI_all_time = tic;
% path2OptionsFile = 'C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab'; %Should go into an excel file
% cd(path2OptionsFile)
% S_OptionsStruct  

%% Get list of CDF files that needs to be subjected to ROI 
cd([Options.Paths.CDF])
% fileList                                          = dir(['*',Options.ROI.MS1_Suffix]);
Options.ROI.mzOptimization.SamplesForOptimization = randi([1,length(fileList)],Options.ROI.mzOptimization.nSamples,1);
if exist(fullfile(Options.Paths.save2mat,'ROI'),'dir') ~= 7
    mkdir(fullfile(Options.Paths.save2mat,'ROI'))
end 
save([Options.Paths.save2mat,'/ROI/OptionsFile_Samples.mat'],"Options")
save([Options.Paths.save2mat,'/ROI/fileList_Samples.mat'],"fileList")

%% Optimize 
%Pre-allocate space 
Options.processing_time.ROI.process.optimization = tic;
% load(fullfile(Options.Paths.save2mat,'ROI','OptionsFile_Samples.mat'))
% load(fullfile(Options.Paths.save2mat,'ROI','fileList_Samples.mat'))

Options.ROI.mzOptimization.nSamples = min(Options.ROI.mzOptimization.nSamples, length(fileList));
[mzerror_Sample,nmz] = deal(zeros(length(Options.ROI.mzOptimization.MZmultiplyFactors),Options.ROI.mzOptimization.nSamples));

parfor n = 1:Options.ROI.mzOptimization.nSamples
    fprintf(1,'Optimization sample: %i/%i\n',n,Options.ROI.mzOptimization.nSamples)
    k                                   =  Options.ROI.mzOptimization.SamplesForOptimization(n);
    [mzerror_Sample(:,n),~,nmz(:,n)]    =  OptParamRoi(fileList(k).name,Options.ROI.mzOptimization.mzerror,Options.ROI.mzOptimization.MZmultiplyFactors,Options.ROI.minroi,Options.ROI);
end

%Calculate the median mzerrror of all the samples used in the optimizationn
Options.ROI.mzerror = median(mzerror_Sample,'all');

Options.processing_time.ROI.process.optimization = toc(Options.processing_time.ROI.process.optimization);
%% Plot for visualisation
figure
clear p
for n = 1:Options.ROI.mzOptimization.nSamples
p(n)  = plot( Options.ROI.mzOptimization.MZmultiplyFactors*Options.ROI.mzOptimization.mzerror,nmz(:,n));
hold on
scatter(mzerror_Sample(1,n),max(nmz(:,n)),'red','filled')

end 
legend(p,{fileList(Options.ROI.mzOptimization.SamplesForOptimization).name})
p = plot(Options.ROI.mzOptimization.MZmultiplyFactors*Options.ROI.mzOptimization.mzerror,median(nmz,2),'LineWidth',2,'LineStyle','--');
% ps = ;

xlabel('\delta_m_/_z','FontAngle','italic')
ylabel('Number of m/z traces / ROIs')

saveas(p,[Options.Paths.save2mat,'/ROI_mz_optim.tif'])


%% ROI processing 
close all
Options.processing_time.ROI.process.sample_wise = tic;

% cd C:\Users\mht541\Documents\MATLAB
% parpool(8)
% [mzroi,MSroi,Rt,Dt,trace] = deal(cell(length(fileList),Options.ROI.NumTrace));
NbrePts          = length(fileList)*Options.ROI.NumTrace;

% Waitbar's Message Field
Msg = ['Sample ',num2str(1),' - ROI Progress...!'];
% Create ParFor Waitbar
[hWaitbar,hWaitbarMsgQueue]= ParForWaitbarCreateMH(Msg,NbrePts);

Options.ROI

addpath (Options.Paths.path2ROI)
fileList_MS2 = dir([[Options.Paths.CDF_MS2,'\*',Options.ROI.MS2_Suffix]]);
if exist(fullfile(Options.Paths.save2mat,'ROI','Samples'),'dir') ~= 7
    mkdir(fullfile(Options.Paths.save2mat,'ROI','Samples'))
end 
for T = 1:Options.ROI.NumTrace
    fprintf(1,'Trace: %i/%i\n',T,Options.ROI.NumTrace)

    parfor k = 1:length(fileList)
        fprintf(1,'Sample: %i/%i\n',k,length(fileList))
%         if T == 1
%             cd(Options.Paths.CDF)
%             FileName = fileList(k).name;
%         else 
%             cd(Options.Paths.CDF_MS2)
%             FileName = fileList_MS2(contains({fileList_MS2.name},fileList(k).name(1:end-Options.ROI.RemoveSuffix))).name;
%         end 
%         hWaitbarMsgQueue.send(0);              
%         % ROI
% %         [mzroi{k,T},MSroi{k,T},~,     ~,      ~,     Rt{k,T},Dt{k,T},trace{k,T}] = ROIpeaks_ACN(FileName,Options.ROI);
%         [mzroi,MSroi,~,     ~,      ~,     Rt,~,trace] = ROIpeaks_ACN(FileName,Options.ROI);
%         save(fullfile(Options.Paths.save2mat,'ROI','Samples', [fileList(k).name(1:end-4),'_MSroi.mat']),'mzroi','MSroi',"Rt","trace","Options")
    processROI(k, T, fileList, fileList_MS2, Options, hWaitbarMsgQueue)
    end

end
%Save to mat
% save([Options.Paths.save2mat,'mzroi.mat'],'mzroi',"Options")
% save([Options.Paths.save2mat,'MSroi.mat'],'mzroi','MSroi',"Rt","Options")

% delete(gcp('nocreate'))
Options.processing_time.ROI.process.sample_wise = toc(Options.processing_time.ROI.process.sample_wise);
%% Augment
Options.processing_time.ROI.process.augmentation = tic;
%Load RoIs for augmentations
[mzroi_aug,MSroi_aug,Rt_aug,~,trace_aug] = deal(cell(length(fileList),Options.ROI.NumTrace));
fileList_mat = dir([Options.Paths.save2mat,'/ROI/Samples/*.mat']);
for k = 1:length(fileList)
    fprintf(1,'Load sample data for augmentation: %i/%i\n',k,length(fileList))
    ind = contains({fileList_mat.name},fileList(k).name(1:end-4));  
    d = load(fullfile(Options.Paths.save2mat,'ROI','Samples',fileList_mat(ind).name));
    % Options.Paths.save2mat = 'C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles';
    mzroi_aug{k} = d.mzroi;
    MSroi_aug{k} = d.MSroi;
    Rt_aug{k} = d.Rt;
    trace_aug{k} = d.trace;
end

% Augmenting 
[mzroi_aug,MSroi_aug]=MSroiaug_ACN2(mzroi_aug,MSroi_aug,Options.ROI);

Options.processing_time.ROI.process.augmentation = toc(Options.processing_time.ROI.process.augmentation);

%%
Options.processing_time.ROI.save = tic; 
%Save to mat
if exist([Options.Paths.save2mat,'/ROI/Samples_aug'],'dir')  ~= 7
    mkdir([Options.Paths.save2mat,'/ROI/Samples_aug'])
end
save([Options.Paths.save2mat,'/ROI/Samples_aug/mzroi_aug.mat'],'mzroi_aug',"Options")
save([Options.Paths.save2mat, '/ROI/Samples_aug/MSroi_aug.mat'], ...
     'mzroi_aug', 'MSroi_aug', 'Rt_aug', 'trace_aug', ...
     'Options', 'fileList_mat', 'fileList', '-v7.3');

Options.processing_time.ROI.save = toc(Options.processing_time.ROI.save);
Options.processing_time.ROI.process.all = toc(ROI_all_time);
clear ROI_all_time


%% functions 
function processROI(k, T, fileList, fileList_MS2, Options, hWaitbarMsgQueue)
% processROI Handles ROI extraction and saving for MS1 / MS2 data
%
% Inputs:
%   k                   - index of current file
%   T                   - mode (1 = MS1, otherwise MS2)
%   fileList            - struct array of MS1 files
%   fileList_MS2        - struct array of MS2 files
%   Options             - struct with paths and ROI settings
%   hWaitbarMsgQueue    - progress queue object

    %% Select correct file and directory
    if T == 1
        cd(Options.Paths.CDF);
        FileName = fileList(k).name;
    else
        cd(Options.Paths.CDF_MS2);

        matchName = fileList(k).name(1:end - Options.ROI.RemoveSuffix);
        idx = contains({fileList_MS2.name}, matchName);

        if ~any(idx)
            error('No matching MS2 file found for %s', fileList(k).name);
        end

        FileName = fileList_MS2(idx).name;
    end

    %% Update waitbar / progress
    if ~isempty(hWaitbarMsgQueue)
        send(hWaitbarMsgQueue, 0);
    end

    %% Run ROI extraction
    [mzroi, MSroi, ~, ~, ~, Rt, ~, trace] = ROIpeaks_ACN(FileName, Options.ROI);

    %% Prepare save path
    saveDir = fullfile(Options.Paths.save2mat, 'ROI', 'Samples');
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end

    %% Output filename
    outName = [fileList(k).name(1:end-4), '_MSroi.mat'];

    %% Save results
    save(fullfile(saveDir, outName), 'mzroi', 'MSroi', 'Rt', 'trace', 'Options');

end
% %% Load Options file 
% tic
% % path2OptionsFile = 'C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab'; %Should go into an excel file
% % cd(path2OptionsFile)
% % S_OptionsStruct  
% 
% %% Get list of CDF files that needs to be subjected to ROI 
% cd([Options.Paths.CDF])
% fileList                                          = dir(['*',Options.ROI.MS1_Suffix]);
% Options.ROI.mzOptimization.SamplesForOptimization = randi([1,length(fileList)],Options.ROI.mzOptimization.nSamples,1);
% if exist(fullfile(Options.Paths.save2mat,'ROI'),'dir') ~= 7
%     mkdir(fullfile(Options.Paths.save2mat,'ROI'))
% end 
% save([Options.Paths.save2mat,'/ROI/OptionsFile_Samples.mat'],"Options")
% save([Options.Paths.save2mat,'/ROI/fileList_Samples.mat'],"fileList")
% 
% %% Optimize 
% %Pre-allocate space 
% load(fullfile(Options.Paths.save2mat,'ROI','OptionsFile_Samples.mat'))
% load(fullfile(Options.Paths.save2mat,'ROI','fileList_Samples.mat'))
% 
% Options.ROI.mzOptimization.nSamples = min(Options.ROI.mzOptimization.nSamples, length(fileList));
% [mzerror_Sample,nmz] = deal(zeros(length(Options.ROI.mzOptimization.MZmultiplyFactors),length(Options.ROI.mzOptimization.SamplesForOptimization)));
% % parpool(3)
% 
% for n = 1:Options.ROI.mzOptimization.nSamples
%     fprintf(1,'Optimization sample: %i/%i\n',n,Options.ROI.mzOptimization.nSamples)
%     k                                   =  Options.ROI.mzOptimization.SamplesForOptimization(n);
%     [mzerror_Sample(:,n),~,nmz(:,n)]    =  OptParamRoi(fileList(k).name,Options.ROI.mzOptimization.mzerror,Options.ROI.mzOptimization.MZmultiplyFactors,Options.ROI.minroi,Options.ROI);
% end
% 
% %Calculate the median mzerrror of all the samples used in the optimizationn
% Options.ROI.mzerror = median(mzerror_Sample,'all');
% 
% %% Plot for visualisation
% figure
% clear p
% for n = 1:Options.ROI.mzOptimization.nSamples
% p(n)  = plot( Options.ROI.mzOptimization.MZmultiplyFactors*Options.ROI.mzOptimization.mzerror,nmz(:,n));
% hold on
% scatter(mzerror_Sample(1,n),max(nmz(:,n)),'red','filled')
% 
% end 
% legend(p,{fileList(Options.ROI.mzOptimization.SamplesForOptimization).name})
% plot(Options.ROI.mzOptimization.MZmultiplyFactors*Options.ROI.mzOptimization.mzerror,median(nmz,2),'LineWidth',2,'LineStyle','--')
% % ps = ;
% 
% xlabel('\delta_m_/_z','FontAngle','italic')
% ylabel('Number of m/z traces / ROIs')
% 
% 
% 
% %% ROI processing 
% close all
% 
% % [mzroi,MSroi,Rt,Dt,trace] = deal(cell(length(fileList),Options.ROI.NumTrace));
% NbrePts          = length(fileList)*Options.ROI.NumTrace;
% 
% % Waitbar's Message Field
% Msg = ['Sample ',num2str(1),' - ROI Progress...!'];
% % Create ParFor Waitbar
% [hWaitbar,hWaitbarMsgQueue]= ParForWaitbarCreateMH(Msg,NbrePts);
% 
% Options.ROI
% 
% addpath (Options.Paths.path2ROI)
% fileList_MS2 = dir([[Options.Paths.CDF_MS2,'\*',Options.ROI.MS2_Suffix]]);
% if exist(fullfile(Options.Paths.save2mat,'ROI','Samples'),'dir') ~= 7
%     mkdir(fullfile(Options.Paths.save2mat,'ROI','Samples'))
% end 
% for T = 1:Options.ROI.NumTrace
%     fprintf(1,'Trace: %i/%i\n',T,Options.ROI.NumTrace)
% 
%     for k = 1:length(fileList)
%         fprintf(1,'Sample: %i/%i\n',k,length(fileList))
%         if T == 1
%             cd(Options.Paths.CDF)
%             FileName = fileList(k).name;
%         else 
%             cd(Options.Paths.CDF_MS2)
%             FileName = fileList_MS2(contains({fileList_MS2.name},fileList(k).name(1:end-Options.ROI.RemoveSuffix))).name;
%         end 
%         hWaitbarMsgQueue.send(0);              
%         % ROI
% %         [mzroi{k,T},MSroi{k,T},~,     ~,      ~,     Rt{k,T},Dt{k,T},trace{k,T}] = ROIpeaks_ACN(FileName,Options.ROI);
%         [mzroi,MSroi,~,     ~,      ~,     Rt,~,trace] = ROIpeaks_ACN(FileName,Options.ROI);
%         save(fullfile(Options.Paths.save2mat,'ROI','Samples', [fileList(k).name(1:end-4),'_MSroi.mat']),'mzroi','MSroi',"Rt","trace","Options")
%     end
% 
% end
% %Save to mat
% % save([Options.Paths.save2mat,'mzroi.mat'],'mzroi',"Options")
% % save([Options.Paths.save2mat,'MSroi.mat'],'mzroi','MSroi',"Rt","Options")
% 
% 
% 
% %% Augment
% 
% %Load RoIs for augmentations
% [mzroi_aug,MSroi_aug,Rt_aug,~,trace_aug] = deal(cell(length(fileList),Options.ROI.NumTrace));
% fileList_mat = dir([Options.Paths.save2mat,'/ROI/Samples/*.mat']);
% for k = 1:length(fileList)
%     fprintf(1,'Load sample data for augmentation: %i/%i\n',k,length(fileList))
%     ind = contains({fileList_mat.name},fileList(k).name(1:end-4));  
%     d = load(fullfile(Options.Paths.save2mat,'ROI','Samples',fileList_mat(ind).name));
%     % Options.Paths.save2mat = 'C:\Users\mht541\Documents\Python_Scripts\2DLC_workflow\matlab\DataFiles';
%     mzroi_aug{k} = d.mzroi;
%     MSroi_aug{k} = d.MSroi;
%     Rt_aug{k} = d.Rt;
%     trace_aug{k} = d.trace;
% end
% 
% % Augmenting 
% [mzroi_aug,MSroi_aug]=MSroiaug_ACN2(mzroi_aug,MSroi_aug,Options.ROI);
% 
% process_time.ROI_process = toc;
% %%
% tic 
% %Save to mat
% if exist([Options.Paths.save2mat,'/ROI/Samples_aug'],'dir')  ~= 7
%     mkdir([Options.Paths.save2mat,'/ROI/Samples_aug'])
% end
% save([Options.Paths.save2mat,'/ROI/Samples_aug/mzroi_aug.mat'],'mzroi_aug',"Options")
% save([Options.Paths.save2mat, '/ROI/Samples_aug/MSroi_aug.mat'], ...
%      'mzroi_aug', 'MSroi_aug', 'Rt_aug', 'trace_aug', ...
%      'Options', 'fileList_mat', 'fileList', '-v7.3');
% 
% process_time.ROI_save = toc;
% 
% 