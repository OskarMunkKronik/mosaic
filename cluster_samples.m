Options.rt_shift =[2,5];
 n_final_group = 1;
for k = 3
    basePeak_vec = find(allPeaks{k}(:,6));
    cos_combined_elu =[];
    cos_mass_spectrum_all = [];
    for n = 1:length(basePeak_vec)
        fG_ref= basePeak_vec(n);
        fC_ref = allPeaks{k}(fG_ref,6);
        rt_k = allPeaks{k}(fG_ref,[2,3]);

        for kk = 4

            candidate_fG =  sum(abs(rt_k-allPeaks{kk}(:,[2,3]))<[Options.rt_shift],2)>0;
            %           n_test_fG = nnz(allPeaks{kk}(candidate_fG,6));
            basePeak_candidate_vec =find(sum([allPeaks{kk}(:,6),candidate_fG]>0,2)>1);
            for nn = 1:length(basePeak_candidate_vec)
                fG_kk =basePeak_candidate_vec(nn);
                fC_kk =allPeaks{kk}(basePeak_candidate_vec(nn),6);
                clf
                %                 fG_kk = basePeak_candidate_vec(nn);
                %                 elu_profile = { sum_elution_profile{k}{fG_ref},...
                %                     sum_elution_profile{kk}{fG_kk} };
                %                 for dim = 1:2
                %                     [~,a] = max(sum_elution_profile{k}{fG_ref}{dim});
                %                     [~,aa] = max(sum_elution_profile{kk}{fG_kk}{dim});
                %                     rt_mov = a-aa; elu = zeros(length(sum_elution_profile{k}{fG_ref}{dim}),1);
                %                     if rt_mov == 0
                %                         elu = sum_elution_profile{k}{fG_ref}{dim};
                %                     elseif rt_mov<0
                %                         elu(abs(rt_mov)+1:end) = sum_elution_profile{k}{fG_ref}{dim}(1:end+rt_mov);
                %                     else elu(1:end-rt_mov) = sum_elution_profile{k}{fG_ref}{dim}(abs(rt_mov)+1:end);
                %                     end
                %                     x = 1:length(sum_elution_profile{k}{fG_ref}{dim});
                %                     subplot(2,1,dim)
                %                     plot(sum_elution_profile{k}{fG_ref}{dim}./max(sum_elution_profile{k}{fG_ref}{dim}),'b')
                %                     hold on
                %                     plot(elu./max(sum_elution_profile{k}{fG_ref}{dim}),'--b')
                %                     plot(sum_elution_profile{kk}{fG_kk}{dim}./max(sum_elution_profile{kk}{fG_kk}{dim}),'r')
                %                     elu_profile{1}{dim} = elu;
                %                 end
                profile_ref = sum_elution_profile{k}{fG_ref};
                profile_cmp = sum_elution_profile{kk}{fG_kk};

                [elu_profile, shifts] = align_elution_profiles(profile_ref, profile_cmp, true);
                elu_profile = {profile_ref,elu_profile};
                [cosMat_tmp,cosMat_combined_tmp] = calculate_cos_matrix(elu_profile);%{sum_elution_profile{k}{n},sum_elution_profile{kk}{fG_kk});

                cos_combined_elu(n,nn) = cosMat_combined_tmp(1,2);

                % Mass spectrum
    clf
                mass_spectrum_temp = zeros(length(mzroi_aug),2);
                mass_spectrum_temp(:,1)= mass_spectrum{k}(:,fC_ref);
                mass_spectrum_temp(:,2)= mass_spectrum{kk}(:,fC_kk);
%                 mass_spectrum_temp = mass_spectrum_temp./max(mass_spectrum_temp);
                stem(mzroi_aug,     mass_spectrum_temp(:,1))
                hold on
                title(allPeaks{k}(n,[2,3,4]))
                stem(mzroi_aug,    - mass_spectrum_temp(:,2))
                cos_mass_spectrum_all(n,nn)= dot(mass_spectrum_temp(:,1),mass_spectrum_temp(:,2)) / (norm(mass_spectrum_temp(:,1))*norm(mass_spectrum_temp(:,2)) + eps);

            end
        end


    end
    cos_final = (cos_combined_elu.*cos_mass_spectrum_all);
    cos_final_filtered = cos_final.*double(cos_final>(1-Options.Clustering.cutoff))
% [val,ind_final_G] =max((cos_final_filtered))
[r,c,val] = find(cos_final_filtered);
for n_group_count = 1:length(r)
    allPeaks{k}(basePeak_vec(r(n_group_count)),7) = n_final_group;
    allPeaks{kk}(basePeak_candidate_vec(c(n_group_count)),7) = n_final_group;
    n_final_group = n_final_group +1;
end 
% find(val)
% allPeaks{k}(basePeak_vec(ind_final_G(val>0)),7) = ind_final_G(val>0)+n_final_group;
% [val,ind_final_G] =max((cos_final_filtered),[],2)
% allPeaks{kk}(basePeak_candidate_vec(basePeak_candidate_vec(val>0)),7) = ind_final_G(val>0)+n_final_group;

end

%%
function [elu_profile, rt_mov_all] = align_elution_profiles(profile_ref, profile_cmp, doPlot)
%ALIGN_ELUTION_PROFILES Aligns comparison profile to match reference
%
% Input:
%   profile_ref : {1x2} cell, reference elution profile (two dimensions)
%   profile_cmp : {1x2} cell, comparison elution profile (two dimensions)
%   doPlot      : logical, true to plot alignment (default = false)
%
% Output:
%   elu_profile : {1x2} cell, aligned version of profile_cmp
%   rt_mov_all  : [1x2] vector, shift amount applied in each dimension

if nargin < 3
    doPlot = false;
end

elu_profile = cell(1,2);
rt_mov_all = zeros(1,2);

for dim = 1:2
    % Find maxima positions
    [~,a]  = max(profile_ref{dim});
    [~,aa] = max(profile_cmp{dim});
    rt_mov = a - aa;           % how much to move cmp to match ref
    rt_mov_all(dim) = rt_mov;

    % Initialize aligned signal
    n = length(profile_cmp{dim});
    elu = zeros(n,1);

    % Shift comparison
    if rt_mov == 0
        elu = profile_cmp{dim};
    elseif rt_mov < 0
        elu(1:end+rt_mov) = profile_cmp{dim}(-rt_mov+1:end);
    else
        elu(rt_mov+1:end) = profile_cmp{dim}(1:end-rt_mov);
    end

    % Save aligned comparison profile
    elu_profile{dim} = elu;

    % --- Optional plotting ---
    if doPlot
        subplot(2,1,dim)
        hold on
        plot(profile_ref{dim}./max(profile_ref{dim}), 'b')
        plot(profile_cmp{dim}./max(profile_cmp{dim}), 'r')
        plot(elu./max(profile_cmp{dim}), '--r')
        hold off
        title(sprintf('Dim %d alignment (shift=%d)', dim, rt_mov))
        legend('Reference','Comparison','Comparison aligned')
    end
end
end
