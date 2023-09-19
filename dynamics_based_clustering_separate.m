clc; clear; %close all;

exps_list = {'26550_1.11','26549_8.11','38428_17.11','38427_1.11','26532_2.3','26550_3.5',...
            '38426_2.11','39740_2.11','26549_11.11','26550_15.11',...
            'control_26550_24.1','control_39740_24.4','control_38427_24.4','control_26536_7.2'};

% exps_list = {'26532_17.11','38428_20.2','broken_8.11','control_38427_7.2'};

expi = 1;

load(['..\data\',exps_list{expi},'\metrics.mat']);

T = 20;

%% *************************************************** BEFORE ************************************************** %%

%% Heatmaps

figure; 
subplot 131
tmp1 = triu(corrs_before_before,1); tmp1 = tmp1(:);
tmp2 = triu(corrs_ICs_before_before,1); tmp2 = tmp2(:); 
histogram2(tmp2(tmp1 ~= 0),tmp1(tmp1 ~= 0),'Normalization','pdf','FaceColor','flat','LineStyle','none'); colormap; view(2);
ylabel('2D corr'); xlabel('IC_{corr}');
title('corr');

subplot 132
tmp1 = triu(com_distance_before_before,1); tmp1 = tmp1(:);
tmp2 = triu(IC_com_distance_before_before,1); tmp2 = tmp2(:); 
histogram2(tmp2(tmp2 ~= 0),tmp1(tmp1 ~= 0),'Normalization','pdf','FaceColor','flat','LineStyle','none'); colormap; view(2);
ylabel('burst_{COM dist}'); xlabel('IC_{COM dist}');
title('COM distance');

subplot 133
tmp1 = triu(spatial_profile_before_before,1); tmp1 = tmp1(:);
tmp2 = triu(IC_spatial_profile_before_before,1); tmp2 = tmp2(:); 
histogram2(tmp2(tmp2 ~= 0),tmp1(tmp1 ~= 0),'Normalization','pdf','FaceColor','flat','LineStyle','none'); colormap; view(2);
ylabel('burst'); xlabel('IC');
title('Spatial profile');

sgtitle('Before');

%% Create distance binary matrix - based on heatmap

factor = 0.9;
% spatiotemporal corr (2D corr)
event_corr_th = 0.81*factor;
ICs_corr_th = 0.75*factor;
adj_matrix_2Dcorr = zeros(n_events_before);
for i = 1:n_events_before
    for j = i+1 : n_events_before
        if (corrs_before_before(i,j) > event_corr_th) && (corrs_ICs_before_before(i,j) > ICs_corr_th)
            adj_matrix_2Dcorr(i,j) = 1;
        end
        adj_matrix_2Dcorr(j,i) = adj_matrix_2Dcorr(i,j); %symmetric
    end
end

% COM
event_com_th = 4*factor;
ICs_com_th = 1*factor;
adj_matrix_COM = zeros(n_events_before);
for i = 1:n_events_before
    for j = i+1 : n_events_before
        if (com_distance_before_before(i,j) < event_com_th) && (IC_com_distance_before_before(i,j) < ICs_com_th)
            adj_matrix_COM(i,j) = 1;
        end
        adj_matrix_COM(j,i) = adj_matrix_COM(i,j); %symmetric
    end
end

% Spatial profile
event_sp_th = 0.87*factor;
ICs_sp_th = 0.7*factor;
adj_matrix_spatial = zeros(n_events_before);
for i = 1:n_events_before
    for j = i+1 : n_events_before
        if (spatial_profile_before_before(i,j) > event_sp_th) && (IC_spatial_profile_before_before(i,j) > ICs_sp_th)
            adj_matrix_spatial(i,j) = 1;
        end
        adj_matrix_spatial(j,i) = adj_matrix_spatial(i,j); %symmetric
    end
end

% graph based on adj matrix
adj_matrix_merged_before = adj_matrix_2Dcorr + adj_matrix_COM + adj_matrix_spatial;

adj_matrix = adj_matrix_merged_before;
S = adj_matrix;

%% plot graph

G = graph(S);

figure; 
plot(G, 'EdgeCData', G.Edges.Weight, 'LineWidth', 3);
colorbar
colormap(flip(gray))

%% Cluster - BEFORE
n_clusters = 15;
valid_cluster_size = 0.02; % [%] of all events

Clusters_before = spectralcluster(S,n_clusters,'Distance','precomputed','LaplacianNormalization','symmetric','ClusterMethod','kmedoids');

clusters_valid_before = {};
figure; 
q = 3;%ceil(sqrt(n_clusters));
i = 1;
for c = 1:n_clusters
    events_idx = find(Clusters_before == c);
    if length(events_idx) > valid_cluster_size*n_events_before
        subplot(q,q,i); hold on; 
        for e = events_idx'
            if size(events_com_before{e},2) >= T
                event = events_com_before{e}(:,1:T);
                plot(event(1,:),event(2,:),'-','color',[.5 .5 .5]);
                plot(event(1,1),event(2,1),'o','color',[.5 .5 .5]);
            end
        end
        xlim([0 b]); ylim([0 a]); grid on; drawnow;
        title([num2str(length(events_idx)),' events (',num2str(100*length(events_idx)/n_events_before),'%)']);
        drawnow;
        clusters_valid_before{i} = events_idx; i = i + 1;
    end
end

n_events_clustered = sum(cellfun(@length,clusters_valid_before));
n_valid_clusters_before = length(clusters_valid_before);

sgtitle({['Valid cluster size: ',num2str(valid_cluster_size*100),'% of all events'];[num2str(100*n_events_clustered/n_events_before),'% of the events are clustered']});

%% Embedding - all spontaneous events - UMAP

events_before_all = zeros(n_events_before,T,N);
for i = 1:n_events_before
    if size(events_before{i},2) >= T
        events_before_all(i,:,:) = events_before{i}(:,1:T)';
    end
end
data = reshape(events_before_all,[],N);

data_umaped = run_umap(data);
data_umaped = reshape(data_umaped,[n_events_before,T,2]); 

hold on;
for c = 1:n_valid_clusters_before
    idx = clusters_valid_before{c}';
    tmp = squeeze(median(data_umaped(idx,:,:)));
    plot(tmp(:,1),tmp(:,2),'-k','LineWidth',1.5);
    plot(tmp(1,1),tmp(1,2),'ok','LineWidth',1.5);
    pause;
end

%% *************************************************** AFTER ************************************************** %%

%% Heatmaps

figure; 
subplot 131
tmp1 = triu(corrs_after_after,1); tmp1 = tmp1(:);
tmp2 = triu(corrs_ICs_after_after,1); tmp2 = tmp2(:); 
histogram2(tmp2(tmp1 ~= 0),tmp1(tmp1 ~= 0),'Normalization','pdf','FaceColor','flat','LineStyle','none'); colormap; view(2);
ylabel('<corr>_t'); xlabel('IC_{corr}');
title('corr');

subplot 132
tmp1 = triu(com_distance_after_after,1); tmp1 = tmp1(:);
tmp2 = triu(IC_com_distance_after_after,1); tmp2 = tmp2(:); 
histogram2(tmp2(tmp2 ~= 0),tmp1(tmp1 ~= 0),'Normalization','pdf','FaceColor','flat','LineStyle','none'); colormap; view(2);
ylabel('burst_{COM dist}'); xlabel('IC_{COM dist}');
title('COM distance');

subplot 133
tmp1 = triu(spatial_profile_after_after,1); tmp1 = tmp1(:);
tmp2 = triu(IC_spatial_profile_after_after,1); tmp2 = tmp2(:); 
histogram2(tmp2(tmp2 ~= 0),tmp1(tmp1 ~= 0),'Normalization','pdf','FaceColor','flat','LineStyle','none'); colormap; view(2);
ylabel('burst'); xlabel('IC');
title('Spatial profile');

sgtitle('After');

%% Create distance binary matrix - based on heatmap 

% spatiotemporal corr (2D corr)
event_corr_th = 0.87;
ICs_corr_th = 0.78;
adj_matrix_2Dcorr = zeros(n_events_after);
for i = 1:n_events_after
    for j = i+1 : n_events_after
        if (corrs_after_after(i,j) > event_corr_th) && (corrs_ICs_after_after(i,j) > ICs_corr_th)
            adj_matrix_2Dcorr(i,j) = 1;
        end
        adj_matrix_2Dcorr(j,i) = adj_matrix_2Dcorr(i,j); %symmetric
    end
end

% COM
event_com_th = 5;
ICs_com_th = 1.5;
adj_matrix_COM = zeros(n_events_after);
for i = 1:n_events_after
    for j = i+1 : n_events_after
        if (com_distance_after_after(i,j) < event_com_th) && (IC_com_distance_after_after(i,j) < ICs_com_th)
            adj_matrix_COM(i,j) = 1;
        end
        adj_matrix_COM(j,i) = adj_matrix_COM(i,j); %symmetric
    end
end

% Spatial profile
event_sp_th = 0.87;
ICs_sp_th = 0.75;
adj_matrix_spatial = zeros(n_events_after);
for i = 1:n_events_after
    for j = i+1 : n_events_after
        if (spatial_profile_after_after(i,j) > event_sp_th) && (IC_spatial_profile_after_after(i,j) > ICs_sp_th)
            adj_matrix_spatial(i,j) = 1;
        end
        adj_matrix_spatial(j,i) = adj_matrix_spatial(i,j); %symmetric
    end
end

% graph based on adj matrix
adj_matrix_merged_after = adj_matrix_2Dcorr + adj_matrix_COM + adj_matrix_spatial;

adj_matrix = adj_matrix_merged_after;
S = adj_matrix;

%% Cluster
n_clusters = 30;
valid_cluster_size = 0.02; % [%] of all events

[Clusters_after,V,D] = spectralcluster(S,n_clusters,'Distance','precomputed','LaplacianNormalization','symmetric','ClusterMethod','kmedoids');

clusters_valid_after = {};
figure; 
q = ceil(sqrt(n_clusters));
i = 1;
for c = 1:n_clusters
    events_idx = find(Clusters_after == c);
    if length(events_idx) > valid_cluster_size*n_events_after
        subplot(q,q,i); hold on; 
        for e = events_idx'
            if size(events_com_after{e},2) >= T
                event = events_com_after{e}(:,1:T);
                plot(event(1,:),event(2,:),'-','color',[.5 .5 .5]);
                plot(event(1,1),event(2,1),'o','color',[.5 .5 .5]);
            end
        end
        xlim([0 b]); ylim([0 a]); grid on; drawnow;
        title([num2str(length(events_idx)),' events (',num2str(100*length(events_idx)/n_events_after),'%)']);
        drawnow;
        clusters_valid_after{i} = events_idx; i = i + 1;
    end
end

n_events_clustered = sum(cellfun(@length,clusters_valid_after));
n_valid_clusters_after = length(clusters_valid_after);

sgtitle({['Valid cluster size: ',num2str(valid_cluster_size*100),'% of all events'];[num2str(100*n_events_clustered/n_events_after),'% of the events are clustered']});

%% ********************************************************** save ***********************************************************************

save(['..\data\left out\',exps_list{expi},'\spontaneous_clusters_separate.mat'],'clusters_valid_before','clusters_valid_after',...
                                 'n_valid_clusters_before','n_valid_clusters_after',...
                                 'n_events_before','n_events_after',...
                                 'adj_matrix_merged_before','adj_matrix_merged_after');
