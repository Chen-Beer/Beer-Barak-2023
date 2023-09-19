function [IC_com_distance,IC_com_corr,com_distance,com_corr,spatial_profile,IC_spatial_profile,temporal_profile,spatiotemporal_corr,IC_corr] = create_metrics_matrices(events1,events1_com,events2,events2_com,T,electrodes,chosen_label)
    
    n_events1 = length(events1);
    if chosen_label
        n_events2 = size(events2,1);
    else
        n_events2 = length(events2);
    end
    
    spatiotemporal_corr = zeros(n_events1,n_events2);
    IC_corr = zeros(n_events1,n_events2);
    com_distance = zeros(n_events1,n_events2);
    com_corr = zeros(n_events1,n_events2);
    IC_com_distance = zeros(n_events1,n_events2);
    IC_com_corr = zeros(n_events1,n_events2);    
    spatial_profile = zeros(n_events1,n_events2);
    IC_spatial_profile = zeros(n_events1,n_events2);
    temporal_profile = zeros(n_events1,n_events2);
%     diagonal_ordering = zeros(n_events1,n_events2);

    parfor ii = 1:n_events1
        if size(events1{ii},2) >= T
            e1 = events1{ii}(electrodes,1:T);
            com1 = events1_com{ii}(:,1:T);
            for jj = 1:n_events2
                if chosen_label
                    e2 = squeeze(events2(jj,electrodes,1:T));
                    com2 = squeeze(events2_com(jj,:,1:T));
                else
                    if size(events2{jj},2) >= T
                        e2 = events2{jj}(electrodes,1:T);
                        com2 = events2_com{jj}(:,1:T);
                    end
                end

                spatiotemporal_corr(ii,jj) = corr2(e1,e2);
                spatial_profile(ii,jj) = corr(sum(e1,2),sum(e2,2));
                IC_spatial_profile(ii,jj) = corr(e1(:,1),e2(:,1));
                temporal_profile(ii,jj) = corr(sum(e1)',sum(e2)');
%                     diagonal_ordering(i,j) = 

                com_distance(ii,jj) = norm(com1-com2);
                com_corr(ii,jj) = corr2(com1,com2);
                IC_com_distance(ii,jj) = norm(com1(:,1)-com2(:,1));
                IC_com_corr(ii,jj) = corr(com1(:,1),com2(:,1));                
                IC_corr(ii,jj) = corr(e1(:,1),e2(:,1));
            end

        end
    end
    
end

