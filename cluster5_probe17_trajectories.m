clc; clear; 

exps_list = {'26550_1.11','26549_8.11','38428_17.11','38427_1.11','26532_2.3','26550_3.5',...
            '38426_2.11','39740_2.11','26549_11.11','26550_15.11',...
            'control_26550_24.1','control_39740_24.4','control_38427_24.4','control_26536_7.2'};

expi = 1;

load(['..\..\..\data\',exps_list{expi},'\metrics.mat']);
load(['..\..\..\..\second protocol\',dir_(4:end),'\2 probe\responsiveness_analysis.mat']);
load(['..\..\..\data\',exps_list{expi},'\spontaneous_clusters_separate.mat']);

T = 20;
events_before_all = zeros(n_events_before,T,N);
for i = 1:n_events_before
    if size(events_before{i},2) >= T
        events_before_all(i,:,:) = events_before{i}(:,1:T)';
    end
end


cluster5 = squeeze(events_before_all(clusters_valid_before{5},:,:));
figure; hold on; 
for i = 1:size(cluster5,1)
    tmp = squeeze(cluster5(i,:,:))';
    tmp_com = create_COM(tmp,electrodes_chosen{3},electrode_names_ordered,electrodes_names);
    plot(tmp_com(1,:),tmp_com(2,:),'color',[.5 .5 .5]);
    plot(tmp_com(1,1),tmp_com(2,1),'o','color',[.5 .5 .5]);
end

probe17 = squeeze(probe_responses_smoothed(17,:,:,1:T));
for i = 1:size(probe17,1)
    tmp = squeeze(probe17(i,:,:));
    tmp_com = create_COM(tmp,electrodes_chosen{3},electrode_names_ordered,electrodes_names);
    plot(tmp_com(1,:),tmp_com(2,:),'color',[.56 .76 .49]);
    plot(tmp_com(1,1),tmp_com(2,1),'o','color',[.56 .76 .49]);
end
xlim([1 b]); ylim([1 a]); grid on;
    