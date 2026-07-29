dist_mat = zeros(size(allPeaks{3},1),size(allPeaks{4},1));
k = 3;
ref = 4;
[~,fG]=min(abs(mzroi_aug(allPeaks{ref}(:,4))-mzTarget));
[dist_mat,cosMat_combined_all,cos_mass_spectrum_all] = deal([]);
% cos_mass_spectrum_all =[];
for k = 3:5
    [~,k_feature,k_vec] = unique(clusters{k});
    for ref = 3:5
        [~,ref_feature,ref_vec] = unique(clusters{ref});
        for n_count = 1:length(k_feature)
            n = k_feature(n_count);

            for nn_count = 1:length(ref_feature)
            nn = ref_feature(nn_count);


                y = allPeaks{k}(n,[2,3])';
                x = allPeaks{ref}(nn,[2,3])';
                dist_mat(n,nn,k,ref) = sqrt(sum((x - y).^2));
                if dist_mat(n,nn,k,ref) <Options.Clustering.distance_componentization
                    clf
                    %Calculate cosine similarity
                    %             clf
                    elu_profile = {sum_elution_profile{k}{n},sum_elution_profile{ref}{nn}};
                    for dim = 1:2
                        [~,a] = max(sum_elution_profile{k}{n}{dim});
                        [~,aa] = max(sum_elution_profile{ref}{nn}{dim});
                        rt_mov = a-aa;
                        elu  = zeros(length(sum_elution_profile{k}{n}{dim}),1);
                        if rt_mov == 0
                            elu =  sum_elution_profile{k}{n}{dim};
                        elseif rt_mov<0
                            elu(abs(rt_mov)+1:end) =  sum_elution_profile{k}{n}{dim}(1:end+rt_mov);
                        else
                            elu(1:end-rt_mov) =  sum_elution_profile{k}{n}{dim}(abs(rt_mov)+1:end);

                        end
                        x = 1:length(sum_elution_profile{k}{n}{dim});
                        %                 subplot(2,1,dim)
                        %                 plot(sum_elution_profile{k}{n}{dim}./max(sum_elution_profile{k}{n}{dim}),'b')
                        %                 hold on
                        %                 plot(elu./max(sum_elution_profile{k}{n}{dim}),'--b')
                        %                 plot(sum_elution_profile{ref}{nn}{dim}./max(sum_elution_profile{ref}{nn}{dim}),'r')
                        elu_profile{1}{dim} = elu;
                    end

                    [cosMat_tmp,cosMat_combined_tmp] = calculate_cos_matrix(elu_profile);%{sum_elution_profile{k}{n},sum_elution_profile{ref}{nn}});
                    cosMat_all{k,ref}{n} = cosMat_tmp;
                    cosMat_combined_all(n,nn,k,ref) = cosMat_combined_tmp(1,2);

                    % mass spectrum 
                    mass_spectrum_temp = zeros(length(mzroi_aug),2);
                     mass_spectrum_temp(:,1)= mass_spectrum{k}(:,clusters{k}(n));
                   mass_spectrum_temp(:,2)= mass_spectrum{ref}(:,clusters{ref}(nn));
                       mass_spectrum_temp = mass_spectrum_temp./max(mass_spectrum_temp);
%                     if k_feature(n_count) == 151 && ref_feature(nn_count) == 158 && ref == 4
% clf
%                     stem(mzroi_aug,     mass_spectrum_temp(:,1),'Marker','none')
%                     hold on 
%                     title(allPeaks{k}(n,[2,3,4]))
%                     stem(mzroi_aug,    - mass_spectrum_temp(:,2),'Marker','none')
% %                     [k,ref,n,nn]
%                clc
%                     end 
                    
                    
                    cos_mass_spectrum_all(n,nn,k,ref) = dot(mass_spectrum_temp(:,1),mass_spectrum_temp(:,2)) / (norm(mass_spectrum_temp(:,1))*norm(mass_spectrum_temp(:,2)) + eps);
                end
            end

            % clf
%             ind   =  find(clusters{ref}(fG) == clusters{ref},1,'first');
%             [r,c] =  find((dist_mat<5));


        end
    end
end 

%%
cos_final = cos_mass_spectrum_all.*cosMat_combined_all;


    %%
    clf
    hold on
    scatter(allPeaks{ref}(c,2),allPeaks{ref}(c,3),'r')

    text(allPeaks{ref}(c,2), allPeaks{ref}(c,3), num2str(clusters{ref}(c)),'Color','r')
    hold on
    scatter(allPeaks{ref}(ind,2),allPeaks{ref}(ind,3),'*')
    scatter(allPeaks{k}(r,2),allPeaks{k}(r,3),'b')

    text(allPeaks{k}(r,2), allPeaks{k}(r,3), num2str(clusters{k}(r)),'Color','b')


%%
% collapse across k,ref (P,Q) by max or sum
hits_collapsed = max(hits, [], [3 4]);   % keep strongest hit

