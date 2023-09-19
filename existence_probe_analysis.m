clc; clear; close all;

exps_list_ = {'26550_1.11','26549_8.11','38428_17.11','38427_1.11','26532_2.3','26550_3.5',...
              '38426_2.11','26549_11.11','26550_15.11',...
              'control_26550_24.1','control_39740_24.4','control_38427_24.4','control_26536_7.2',...
              'left out\38428_20.2','left out\broken_8.11','left out\control_38427_7.2'};

% EXCLUDE NON-EXISTING STIMULATED EVOKED
% existing_evoked_list = [[1  1  1];[0  1  0];[1  0  0];[1  0  1];[0  1  0];[0  0  1];...
%                         [0  1  0];[1  0  0];[0  1  1];...
%                         [1  1  1];[0  0  1];[1  0  1];[1  1  1];...
%                         [1  1  1];[1  0  1];[1  0  1]];

%% load

for expi = 1:length(exps_list_)
    expi
load(['..\data\',exps_list_{expi},'\metrics.mat'],'dir_','events_com_before','events_before');
load(['..\data\',exps_list_{expi},'\part2.mat'],'similarity_corr_th');
load(['..\data\',exps_list_{expi},'\spontaneous_clusters_separate.mat']);
load(['..\..\..\second protocol\',dir_(4:end),'\2 probe\responsiveness_analysis.mat']);

T = 20;
Np = size(responses_smoothed,1);
nn = size(responses_smoothed,2);

%% non stimulating electrodes
non_stimulating_electrodes_idx = cell(Np,1);
all_electrodes = 1:N;
for p = 1:Np
    non_stimulating_electrodes_idx{p} = all_electrodes(~ismember(all_electrodes,stimulus_electrods_idx{p}));
end

%% within 2D correlations
corr_probe_all = zeros(Np,nn,nn);
parfor p = 1:Np
    res_ = squeeze(responses_smoothed(p,:,non_stimulating_electrodes_idx{p},1:T));
    for i = 1:nn
        res1 = squeeze(res_(i,:,:));
        for j = 1:nn
            res2 = squeeze(res_(j,:,:));
            corr_probe_all(p,i,j) = corr2(res1,res2);
        end
    end
end

% plot
median_corr_probe_all = zeros(Np,1);
std_corr_probe_all = zeros(Np,1);
figure; 
for p = 1:Np
    subplot(5,4,p);
    tmp = triu(squeeze(corr_probe_all(p,:,:)),1);
    tmp = tmp(tmp~=0);
    histogram(tmp,0:.05:1);
    title([num2str(p)]);%,' | median corr = ',num2str(median(tmp))]);
    median_corr_probe_all(p) = median(tmp);
    std_corr_probe_all(p) = std(tmp);
end

%% robust responses - plot
factor_th = .95;
robust_probes = find(median_corr_probe_all > factor_th*similarity_corr_th);
non_robust_probes = find(median_corr_probe_all <= factor_th*similarity_corr_th);

% figure; 
% q = ceil(sqrt(length(robust_probes)));
% for pi = 1:length(robust_probes)
%     p = robust_probes(pi);
%     subplot(q,q,pi); hold on; 
%     for e = 1:nn
%         event = squeeze(responses_smoothed(p,e,:,1:T));
%         event_com = create_COM(event,non_stimulating_electrodes_idx{p},electrode_names_ordered,electrodes_names);
%         plot(event_com(1,:),event_com(2,:),'-','color',[.5 .5 .5]);
%         plot(event_com(1,1),event_com(2,1),'o','color',[.5 .5 .5]);
%     end
%     xlim([0 b]); ylim([0 a]); grid on; 
%     title([num2str(p)]);
%     drawnow;
% end

%% plot spont vocab
% figure; 
% q = ceil(sqrt(n_valid_clusters_before));
% for c = 1:n_valid_clusters_before
%     events_idx = clusters_valid_before{c};
%     subplot(q,q,c); hold on; 
%     for e = events_idx'
%         if size(events_com_before{e},2) >= T
%             event = events_com_before{e}(:,1:T);
%             plot(event(1,:),event(2,:),'-','color',[.5 .5 .5]);
%             plot(event(1,1),event(2,1),'o','color',[.5 .5 .5]);
%         end
%     end
%     xlim([0 b]); ylim([0 a]); grid on; 
%     title([num2str(length(events_idx)),' events (',num2str(100*length(events_idx)/n_events_before),'%)']);
%     drawnow;
% end

%% within\across 2D correlations
% for pi = 1:length(robust_probes)
%     p = robust_probes(pi);
%     figure; sgtitle(['Probe #',num2str(p)]);
%     for c = 1:n_valid_clusters_before
%         events_idx = clusters_valid_before{c};
%         cluster_within = zeros(length(events_idx),length(events_idx));
%         cluster_probe = zeros(length(events_idx),nn);
%         for e1 = 1:length(events_idx)
%             if size(events_before{events_idx(e1)},2) >= T
%                 spont_ = events_before{events_idx(e1)}(non_stimulating_electrodes_idx{p},1:T);
%                 for ep = 1:nn
%                     probe_ = squeeze(responses_smoothed(p,ep,non_stimulating_electrodes_idx{p},1:T));
%                     cluster_probe(e1,ep) = corr2(spont_,probe_);
%                 end
%                 for e2 = 1:length(events_idx)
%                     if size(events_before{events_idx(e2)},2) >= T
%                         spont2_ = events_before{events_idx(e2)}(non_stimulating_electrodes_idx{p},1:T);
%                         cluster_within(e1,e2) = corr2(spont_,spont2_);
%                     end
%                 end
%             end
%         end
%         subplot(q,q,c); hold on; 
%         histogram(cluster_within,0:.01:1,'Normalization','pdf');
%         histogram(cluster_probe,0:.01:1,'Normalization','pdf');
%         title([num2str(median(cluster_within(:))),' within | ',num2str(median(cluster_probe(:))),' probe']);
%         drawnow;
%     end
% end

%% cluster idx all

clusters_idx_all = zeros(n_events_before,1);
for c = 1:n_valid_clusters_before
    clust_idx = clusters_valid_before{c};
    clusters_idx_all(clust_idx) = c;
end

%% relate each evoked to spont cluster(s)

similar_clusters_to_probe = cell(Np,1);

for p = 1:Np
%     disp(p);
    similar_clusters_to_probe{p} = [];
    probe_ = squeeze(responses_smoothed(p,:,non_stimulating_electrodes_idx{p},1:T));
    for ep = 1:nn
        evoked_ = squeeze(probe_(ep,:,:));
        labels_ = zeros(n_events_before,1);
        parfor es = 1:n_events_before
            if size(events_before{es},2) >= T
                spont_ = events_before{es}(non_stimulating_electrodes_idx{p},1:T);
                if corr2(evoked_,spont_) > .88*similarity_corr_th
                    labels_(es) = 1;
                end
            end
        end
        if sum(labels_) < 0.005*n_events_before
            similar_clusters_to_probe{p} = cat(1,similar_clusters_to_probe{p},nan);
        else
            [c,g] = groupcounts(clusters_idx_all(labels_ == 1));
            for i = 1:length(g)
                c(i) = c(i)/sum(clusters_idx_all == g(i));
            end
            [~,i] = max(c);
            similar_clusters_to_probe{p} = cat(1,similar_clusters_to_probe{p},g(i));
%             similar_clusters_to_probe{p} = cat(1,similar_clusters_to_probe{p},clusters_idx_all(labels_ == 1));
        end
    end
end

%% histograms - relation to vocab clusters

% figure; 
% for p = 1:Np
%     subplot(5,4,p);
%     histogram(similar_clusters_to_probe{p},-.5:.5:n_valid_clusters_before+.5)
%     title(p); xlabel('Cluster #'); ylabel('similar events')
%     ylim([0 nn])
% end

%% now many clusters explain X% of the probe responses?
th_ = 0.7;
count_clusters = inf(Np,1);
for p = 1:Np
    [c,g] = groupcounts(similar_clusters_to_probe{p});
    if ~all(isnan(g))
        if isnan(g(end))
            c(end) = [];
            g(end) = [];
        end
        [c,i] = sort(c,'descend'); g = g(i);
        gg = 1;
        cc = c(gg);
        while cc < th_*nn && gg < length(g)
            gg = gg + 1;
            cc = cc + c(gg);
        end
        if cc >= th_*nn
            count_clusters(p) = gg;
        end
    end
end

% figure; 
% count_clusters_ = count_clusters;
% count_clusters_(count_clusters == inf) = -1;
% plot(median_corr_probe_all,count_clusters_,'x');

%% Median overall activity - probe
median_overall_activity = zeros(Np,1);
for p = 1:Np
    res_ = squeeze(responses_smoothed(p,:,non_stimulating_electrodes_idx{p},1:T));
    median_overall_activity(p) = median(sum(squeeze(sum(res_,2)),2));
end

%% main_similar_clusters

main_similar_clusters = zeros(length(chosen_patterns),1);
for p = 1:length(chosen_patterns)
    main_similar_clusters(p) = mode(similar_clusters_to_probe{chosen_patterns(p)});
end

%% main non-similar clusters

tmp = unique(cat(1,similar_clusters_to_probe{chosen_patterns}));
tmp = tmp(~isnan(tmp));

clusters_vec = 1:n_valid_clusters_before;

main_not_similar_clusters = clusters_vec(~ismember(clusters_vec,tmp));

%% save 
save(['..\data\',exps_list_{expi},'\evoked_existence_clusters_tmp.mat'],'main_not_similar_clusters','main_similar_clusters','similar_clusters_to_probe','chosen_patterns','median_overall_activity','clusters_idx_all','count_clusters','median_corr_probe_all','std_corr_probe_all','robust_probes','non_robust_probes');
drawnow;

disp(exps_list_{expi});
disp('Similar clusters:')
disp(num2str(main_similar_clusters'))
disp('Not similar clusters:')
disp(num2str(main_not_similar_clusters))
disp('*****************')

end

%% **************************** plot all ****************************************

exps_list_ = {'26550_1.11','26549_8.11','38428_17.11','38427_1.11','26532_2.3','26550_3.5',...
              '38426_2.11','26549_11.11','26550_15.11',...
              'control_26550_24.1','control_39740_24.4','control_38427_24.4','control_26536_7.2',...
              '38428_20.2','broken_8.11','control_38427_7.2'};

figure; hold on; 
for expi = 1:length(exps_list_)
    load(['..\data\',exps_list_{expi},'\evoked_existence_clusters_tmp.mat']);
    load(['..\data\',exps_list_{expi},'\metrics.mat'],'chosen_patterns');
    count_clusters_ = count_clusters + (rand(length(count_clusters),1)-0.5)*.3;
    count_clusters_(count_clusters == inf) = 10 + (rand(sum(count_clusters == inf),1)-0.5)*.3;
    scatter(median_corr_probe_all,count_clusters_,5,'k','filled');
    if ~isempty((median_corr_probe_all > 0.75) & (count_clusters_ > 6))
        disp(['**', exps_list_{expi},': '])
        disp(find((median_corr_probe_all > 0.75) & (count_clusters_ > 6)))
        disp(count_clusters(find((median_corr_probe_all > 0.75) & (count_clusters_ > 6))))
    end
    for p = 1:3
%         if existing_evoked_list(expi,p)
            plot(median_corr_probe_all(chosen_patterns(p)),count_clusters_(chosen_patterns(p)),'or');
%         end
    end
%     check_idx = find(count_clusters < 4 & median_corr_probe_all > 0.75)';
%     disp([exps_list_{expi},' : ',num2str((ismember(chosen_patterns,check_idx)))])
end
xlabel('Robustness: Probe within correlation');
ylabel('# clusters explaining 50% of responses')
