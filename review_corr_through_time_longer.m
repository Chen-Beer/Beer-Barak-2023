clc; clear; 

exps_list = {'26550_1.11','26549_8.11','38428_17.11','38427_1.11','26532_2.3','26550_3.5',...
            '38426_2.11','39740_2.11','26549_11.11','26550_15.11',...
            'control_26550_24.1','control_39740_24.4','control_38427_24.4','control_26536_7.2'};

expi = 1;

load(['..\data\',exps_list{expi},'\metrics.mat']);
load(['..\data\',exps_list{expi},'\spontaneous_clusters_separate.mat']);
load(['..\data\',exps_list{expi},'\part2.mat'],'similarity_corr_th');

%% tmp_binary_mat
min_T = 30;
tmp_binary_mat = zeros(n_events_before);
parfor i = 1:n_events_before
    disp(i)
    for j = 1:n_events_before
        if corrs_events_before_before(i,j) > similarity_corr_th
            if min(size(events_before{i},2),size(events_before{j},2)) > min_T && j>i
                tmp_binary_mat(i,j) = 1;
            end
        end
    end
end

nn = sum(sum(tmp_binary_mat));

%% corr through time for pair of similar bursts (2D corr > th)
temporal_corrs = zeros(nn,min_T);
ni = 1;
for i = 1:n_events_before
    for j = i+1:n_events_before
        if tmp_binary_mat(i,j)
            event1 = events_before{i};
            event2 = events_before{j};
            parfor t = 1:min_T
                temporal_corrs(ni,t) = corr(event1(:,t),event2(:,t));
            end
            ni = ni + 1;
        end
    end
end

%% plot

figure;
shadedErrorBar(dt*(1:min_T),mean(temporal_corrs),std(temporal_corrs));

