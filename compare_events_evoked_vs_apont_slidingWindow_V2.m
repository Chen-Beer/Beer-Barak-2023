% This script extracts events (network spikes) from a given time series
% Events categories:
% 1) spontaneous before experiment
% 2) true evoked - sliding window responses (hand picked "responses" after a stimulus (labeled according to pattern number)
% 3) spontaneous during experiment (1 min between stimulations)
% 3) spontaneous after experiment

clc;  clear;
% close all;
plot_flag = 0;
clean_spont_events = 0;

%% data paths 1
% path_before = '../data/23days_27.10/2 spon/data.mat';
% path_exp = '../data/23days_27.10/3 exp/analysis_6metrics.mat';
% path_after = '../data/23days_27.10/4 spon/data.mat';
% 
% path_save = '../data/23days_27.10/';

%% data paths 2
% path_before = '../data/18days_4.11/1 spont/data.mat';
% path_exp = '../data/18days_4.11/4 exp/analysis_6metrics.mat';
% 
% path_save = '../data/18days_4.11/';

%% data paths 3
% path_before = '../data/39740_1.11/17days_18.11/1 spont/data.mat';
% path_exp = '../data/39740_1.11/17days_18.11/4 exp/analysis_6metrics.mat';
% path_after = '../data/39740_1.11/17days_18.11/5 spont/data.mat';
% 
% path_save = '../data/39740_1.11/17days_18.11/';

%% data paths 4
% path_before = '../data/26550_1.11/19days_20.11/1 spont/data.mat';
% path_exp = '../data/26550_1.11/19days_20.11/4 exp/analysis_6metrics.mat';
% path_after = '../data/26550_1.11/19days_20.11/5 spont/data.mat';
% 
% path_save = '../data/26550_1.11/19days_20.11/';

%% data paths 5
% path_before = '../data/38427_1.11/20days_21.11/1 spont/data.mat';
% path_exp = '../data/38427_1.11/20days_21.11/4 exp/analysis_6metrics.mat';
% path_after = '../data/38427_1.11/20days_21.11/5 spont/data.mat';
% 
% path_save = '../data/38427_1.11/20days_21.11/';

%% data paths 6
% path_before = '../data/26549_8.11/1 spont/data.mat';
% path_exp = '../data/26549_8.11/4 exp/analysis_6metrics.mat';
% path_after = '../data/26549_8.11/5 spont/data.mat';
% 
% path_save = '../data/26549_8.11/';

%% data paths 7
% path_before = '../data/broken_8.11/1 spont/data.mat';
% path_exp = '../data/broken_8.11/4 exp/analysis_6metrics.mat';
% path_after = '../data/broken_8.11/5 spont/data.mat';
% 
% path_save = '../data/broken_8.11/';

%% data paths 8
% path_before = '../data/38428_17.11/1 spont/data.mat';
% path_exp = '../data/38428_17.11/4 exp/analysis_6metrics.mat';
% path_after = '../data/38428_17.11/5 spont/data.mat';
% 
% path_save = '../data/38428_17.11/';

%% data paths 9
% path_before = '../data/26532_17.11/1 spont/data.mat';
% path_exp = '../data/26532_17.11/4 exp/analysis_6metrics.mat';
% path_after = '../data/26532_17.11/5 spont/data.mat';
% 
% path_save = '../data/26532_17.11/';

%% data paths 10
% path_before = '../data/39740_22.11/1 spont/data.mat';
% path_exp = '../data/39740_22.11/4 exp/analysis_6metrics.mat';
% path_after = '../data/39740_22.11/5 spont/data.mat';
% 
% path_save = '../data/39740_22.11/';

%% data paths 11
% path_before = '../data/26550_22.11/1 spont/data.mat';
% path_exp = '../data/26550_22.11/4 exp/analysis_6metrics.mat';
% path_after = '../data/26550_22.11/5 spont/data.mat';
% 
% path_save = '../data/26550_22.11/';

%% data paths 12
% path_before = '../data/39740_13.12/1 spont/data.mat';
% path_exp = '../data/39740_13.12/4 exp/analysis_6metrics.mat';
% path_after = '../data/39740_13.12/5 spont/data.mat';
% 
% path_save = '../data/39740_13.12/';

%% data paths 13
path_before = '../data/38428_13.12/1 spont/data.mat';
path_exp = '../data/38428_13.12/4 exp/analysis_6metrics.mat';
path_after = '../data/38428_13.12/5 spont/data.mat';

path_save = '../data/38428_13.12/';

%% electrode name to 2D index
a = 10; b = 12;
electrode_names_ordered = cell(a,b);
for i = 1:a
    for j = 1:b
        electrode_names_ordered{i,j} = ['A'+j-1,'0','1'+i-1];
        if i == a
            electrode_names_ordered{i,j} = char(['A'+j-1,'10']);
        end
        if j >= 9
            electrode_names_ordered{i,j} = char(['A'+j,'0','1'+i-1]);
            if i == a
                electrode_names_ordered{i,j} = char(['A'+j,'10']);
            end
        end
    end
end
    
%% spontaneous - before
disp('** Spontaneous activity - before: **');
disp('Loading data...');
load(path_before)
activity_smoothed_before = activity_smoothed;
chosen_electrodes = 1:N; % no 
% stimulating electrodes
[time_frames_before,events_smoothed_before,events_2D_before,events_CoM_before] = event_detection(activity_smoothed_before,dt,electrode_names_ordered,electrodes_names,chosen_electrodes,plot_flag,clean_spont_events);

% within: 
disp('Calculating pairwise distances...');
f = waitbar(0);
plot_flag = 0;
n_events = length(time_frames_before);
max_cc_before = nan(n_events);
com_mean_dist_before = nan(n_events);
temporal_corr_before = nan(n_events);
spatial_corr_before = nan(n_events);
diagonal_ordering_before = nan(n_events);
% EMD_before = nan(n_events);
for i = 1:n_events
    waitbar(i/n_events,f,sprintf([num2str(i),'/',num2str(n_events),' completed']))
    event1 = events_smoothed_before{i};
    com1 = events_CoM_before{i};
    parfor j = i:n_events
        event2 = events_smoothed_before{j};
        com2 = events_CoM_before{j};
%         [max_cc_before(i,j),com_mean_dist_before(i,j),temporal_corr_before(i,j),spatial_corr_before(i,j),diagonal_ordering_before(i,j),EMD_before(i,j)] = event_comparison_6metrics(event1,event2,com1,com2,plot_flag)
        [max_cc_before(i,j),com_mean_dist_before(i,j),temporal_corr_before(i,j),spatial_corr_before(i,j),diagonal_ordering_before(i,j)] = event_comparison_6metrics(event1,event2,com1,com2,plot_flag)
    end
end
% flip matrix to be symmetric
max_cc_before = tril(max_cc_before.',-1) + triu(max_cc_before);
com_mean_dist_before = tril(com_mean_dist_before.',-1) + triu(com_mean_dist_before);
temporal_corr_before = tril(temporal_corr_before.',-1) + triu(temporal_corr_before);
spatial_corr_before = tril(spatial_corr_before.',-1) + triu(spatial_corr_before);
diagonal_ordering_before = tril(diagonal_ordering_before.',-1) + triu(diagonal_ordering_before);
% EMD_before = tril(EMD_before.',-1) + triu(EMD_before);
delete(f);

figure; 
subplot 231; hold on; histogram(nonzeros(triu(max_cc_before,1)),'Normalization','probability'); title('max cross-corr');
subplot 232; hold on; histogram(nonzeros(triu(com_mean_dist_before,1)),'Normalization','probability'); title('CoM - mean euclidean distance');
subplot 233; hold on; histogram(nonzeros(triu(temporal_corr_before,1)),'Normalization','probability'); title('Temporal profile - corr');
subplot 234; hold on; histogram(nonzeros(triu(spatial_corr_before,1)),'Normalization','probability'); title('Spatial profile - corr');
subplot 235; hold on; histogram(nonzeros(triu(diagonal_ordering_before,1)),'Normalization','probability'); title('Diagonal ordering');
% subplot 236; hold on; histogram(nonzeros(triu(EMD_before,1)),'Normalization','probability'); title('EMD');
drawnow;

%% true evoked - experiment (sliding window)
disp('** Evoked (stimulation): **');
disp('Loading data...');
load(path_exp)

% converting the responses (true evoked sliding window events) to the relevant format
% window = 20 repetitions
% step = 10 repetitions (overlap between windows)
disp('Extracting true evoked events...');
events_smoothed_evoked = permute(sliding_responses,[2 1 3 4]);
n_events = size(events_smoothed_evoked,2);

% center of mass
disp('Extracting center of mass trajectories...');
E = length(chosen_electrodes);
T = size(events_smoothed_evoked,4);
events_CoM_evoked = permute(sliding_CoM,[2 1 3 4]);

events_CoM_evoked_ = zeros(n_patterns*n_events,2,T);
events_smoothed_evoked_ = zeros(n_patterns*n_events,E,T);
for p = 1:n_patterns
    idx = (p-1)*n_events + (1:n_events);
    events_CoM_evoked_(idx,:,:) = squeeze(events_CoM_evoked(p,:,:,:));
    events_smoothed_evoked_(idx,:,:) = squeeze(events_smoothed_evoked(p,:,:,:));
end

labels_evoked_patterns = zeros(1,n_events*n_patterns);
for p = 1:n_patterns
    idx = (p-1)*n_events + (1:n_events);
    labels_evoked_patterns(idx) = p;
end
labels_evoked_repetitions = repmat(1:n_events,1,n_patterns);

% within:
disp('Calculating pairwise distances...');
f = waitbar(0);
plot_flag = 0;
max_cc_evoked = nan(n_events*n_patterns);
com_mean_dist_evoked = nan(n_events*n_patterns);
temporal_corr_evoked = nan(n_events*n_patterns);
spatial_corr_evoked = nan(n_events*n_patterns);
diagonal_ordering_evoked = nan(n_events);
% EMD_evoked = nan(n_events);
for i = 1:(n_events*n_patterns)
    waitbar(i/(n_events*n_patterns),f,sprintf([num2str(i),'/',num2str(n_events*n_patterns),' completed']))
    event1 = squeeze(events_smoothed_evoked_(i,:,:));
    com1 = squeeze(events_CoM_evoked_(i,:,:));
    if max(max(abs(com1)))>0 && max(max(abs(event1)))>0
        parfor j = i:(n_events*n_patterns)
            event2 = squeeze(events_smoothed_evoked_(j,:,:));
            com2 = squeeze(events_CoM_evoked_(j,:,:));
            if max(max(abs(com2)))>0 && max(max(abs(event2)))>0
%                 [max_cc_evoked(i,j),com_mean_dist_evoked(i,j),temporal_corr_evoked(i,j),spatial_corr_evoked(i,j),diagonal_ordering_evoked(i,j),EMD_evoked(i,j)] = event_comparison_6metrics(event1,event2,com1,com2,plot_flag);
                [max_cc_evoked(i,j),com_mean_dist_evoked(i,j),temporal_corr_evoked(i,j),spatial_corr_evoked(i,j),diagonal_ordering_evoked(i,j)] = event_comparison_6metrics(event1,event2,com1,com2,plot_flag);                
            end
        end
    end
end
delete(f);

chosen_electrodes_evoked = chosen_electrodes;

% flip matrix to be symmetric
max_cc_evoked = tril(max_cc_evoked.',-1) + triu(max_cc_evoked);
com_mean_dist_evoked = tril(com_mean_dist_evoked.',-1) + triu(com_mean_dist_evoked);
temporal_corr_evoked = tril(temporal_corr_evoked.',-1) + triu(temporal_corr_evoked);
spatial_corr_evoked = tril(spatial_corr_evoked.',-1) + triu(spatial_corr_evoked);
diagonal_ordering_evoked = tril(diagonal_ordering_evoked.',-1) + triu(diagonal_ordering_evoked);
% EMD_evoked = tril(EMD_evoked.',-1) + triu(EMD_evoked);

subplot 231; histogram(nonzeros(triu(max_cc_evoked,1)),'Normalization','probability');
subplot 232; histogram(nonzeros(triu(com_mean_dist_evoked,1)),'Normalization','probability');
subplot 233; histogram(nonzeros(triu(temporal_corr_evoked,1)),'Normalization','probability');
subplot 234; histogram(nonzeros(triu(spatial_corr_evoked,1)),'Normalization','probability');
subplot 235; histogram(nonzeros(triu(diagonal_ordering_evoked,1)),'Normalization','probability');
% subplot 236; histogram(nonzeros(triu(EMD_evoked,1)),'Normalization','probability');
drawnow;

%% spontaneous - during experiment
disp('** Spontaneous activity - during experiment: **');
activity_smoothed_during = spont_activity_smoothed;
chosen_electrodes = 1:N; % no stimulating electrodes
[time_frames_during,events_smoothed_during,events_2D_during,events_CoM_during] = event_detection(activity_smoothed_during,dt,electrode_names_ordered,electrodes_names,chosen_electrodes,plot_flag,clean_spont_events);

% within: 
disp('Calculating pairwise distances...');
f = waitbar(0);
plot_flag = 0;
n_events = length(time_frames_during);
max_cc_during = nan(n_events);
com_mean_dist_during = nan(n_events);
temporal_corr_during = nan(n_events);
spatial_corr_during = nan(n_events);
diagonal_ordering_during = nan(n_events);
% EMD_during = nan(n_events);
for i = 1:n_events
    waitbar(i/n_events,f,sprintf([num2str(i),'/',num2str(n_events),' completed']))
    event1 = events_smoothed_during{i};
    com1 = events_CoM_during{i};
    parfor j = i:n_events
        event2 = events_smoothed_during{j};
        com2 = events_CoM_during{j};
%         [max_cc_during(i,j),com_mean_dist_during(i,j),temporal_corr_during(i,j),spatial_corr_during(i,j),diagonal_ordering_during(i,j),EMD_during(i,j)] = event_comparison_6metrics(event1,event2,com1,com2,plot_flag);
        [max_cc_during(i,j),com_mean_dist_during(i,j),temporal_corr_during(i,j),spatial_corr_during(i,j),diagonal_ordering_during(i,j)] = event_comparison_6metrics(event1,event2,com1,com2,plot_flag);
    end
end
delete(f);

labels_time_during = time_frames_during(:,1);

% flip matrix to be symmetric
max_cc_during = tril(max_cc_during.',-1) + triu(max_cc_during);
com_mean_dist_during = tril(com_mean_dist_during.',-1) + triu(com_mean_dist_during);
temporal_corr_during = tril(temporal_corr_during.',-1) + triu(temporal_corr_during);
spatial_corr_during = tril(spatial_corr_during.',-1) + triu(spatial_corr_during);
diagonal_ordering_during = tril(diagonal_ordering_during.',-1) + triu(diagonal_ordering_during);
% EMD_during = tril(EMD_during.',-1) + triu(EMD_during);

subplot 231; histogram(nonzeros(triu(max_cc_during,1)),'Normalization','probability');
subplot 232; histogram(nonzeros(triu(com_mean_dist_during,1)),'Normalization','probability');
subplot 233; histogram(nonzeros(triu(temporal_corr_during,1)),'Normalization','probability');
subplot 234; histogram(nonzeros(triu(spatial_corr_during,1)),'Normalization','probability');
subplot 235; histogram(nonzeros(triu(diagonal_ordering_during,1)),'Normalization','probability');
% subplot 236; histogram(nonzeros(triu(EMD_during,1)),'Normalization','probability');
drawnow;


%% spontaneous - after experiment
disp('** Spontaneous activity - after: **');
load(path_after)
activity_smoothed_after = activity_smoothed;
chosen_electrodes = 1:N; % no stimulating electrodes
[time_frames_after,events_smoothed_after,events_2D_after,events_CoM_after] = event_detection(activity_smoothed_after,dt,electrode_names_ordered,electrodes_names,chosen_electrodes,plot_flag,clean_spont_events);

% within: 
disp('Calculating pairwise distances...');
f = waitbar(0);
plot_flag = 0;
n_events = length(time_frames_after);
max_cc_after = nan(n_events);
com_mean_dist_after = nan(n_events);
temporal_corr_after = nan(n_events);
spatial_corr_after = nan(n_events);
diagonal_ordering_after = nan(n_events);
% EMD_after = nan(n_events);
for i = 1:n_events
    waitbar(i/n_events,f,sprintf([num2str(i),'/',num2str(n_events),' completed']))
    event1 = events_smoothed_after{i};
    com1 = events_CoM_after{i};
    parfor j = i:n_events
        event2 = events_smoothed_after{j};
        com2 = events_CoM_after{j};
%         [max_cc_after(i,j),com_mean_dist_after(i,j),temporal_corr_after(i,j),spatial_corr_after(i,j),diagonal_ordering_after(i,j),EMD_after(i,j)] = event_comparison_6metrics(event1,event2,com1,com2,plot_flag);
        [max_cc_after(i,j),com_mean_dist_after(i,j),temporal_corr_after(i,j),spatial_corr_after(i,j),diagonal_ordering_after(i,j)] = event_comparison_6metrics(event1,event2,com1,com2,plot_flag);
    end
end
delete(f);
% flip matrix to be symmetric
max_cc_after = tril(max_cc_after.',-1) + triu(max_cc_after);
com_mean_dist_after = tril(com_mean_dist_after.',-1) + triu(com_mean_dist_after);
temporal_corr_after = tril(temporal_corr_after.',-1) + triu(temporal_corr_after);
spatial_corr_after = tril(spatial_corr_after.',-1) + triu(spatial_corr_after);
diagonal_ordering_after = tril(diagonal_ordering_after.',-1) + triu(diagonal_ordering_after);
% EMD_after = tril(EMD_after.',-1) + triu(EMD_after);

subplot 231; histogram(nonzeros(triu(max_cc_after,1)),'Normalization','probability');
subplot 232; histogram(nonzeros(triu(com_mean_dist_after,1)),'Normalization','probability');
subplot 233; histogram(nonzeros(triu(temporal_corr_after,1)),'Normalization','probability');
subplot 234; histogram(nonzeros(triu(spatial_corr_after,1)),'Normalization','probability');
subplot 235; histogram(nonzeros(triu(diagonal_ordering_after,1)),'Normalization','probability');
% subplot 236; histogram(nonzeros(triu(EMD_after,1)),'Normalization','probability');
legend('spont. before','evoked','spont. during','spont. after');
drawnow;

%% across: compare the 3 classes

% before - after
disp('1 Calculating pairwise distances: before & after...');
f = waitbar(0);
max_cc_before_after = nan(length(time_frames_before),length(time_frames_after));
com_mean_dist_before_after = nan(length(time_frames_before),length(time_frames_after));
temporal_corr_before_after = nan(length(time_frames_before),length(time_frames_after));
spatial_corr_before_after = nan(length(time_frames_before),length(time_frames_after));
diagonal_ordering_before_after = nan(length(time_frames_before),length(time_frames_after));
% EMD_before_after = nan(length(time_frames_before),length(time_frames_after));
for i = 1:length(time_frames_before)
    waitbar(i/length(time_frames_before),f,sprintf([num2str(i),'/',num2str(length(time_frames_before)),' completed']))
    event1 = events_smoothed_before{i};
    com1 = events_CoM_before{i};
    parfor j = 1:length(time_frames_after)
        event2 = events_smoothed_after{j};
        com2 = events_CoM_after{j};
%         [max_cc_before_after(i,j),com_mean_dist_before_after(i,j),temporal_corr_before_after(i,j),spatial_corr_before_after(i,j),diagonal_ordering_before_after(i,j),EMD_before_after(i,j)] = event_comparison_6metrics(event1,event2,com1,com2,plot_flag);
        [max_cc_before_after(i,j),com_mean_dist_before_after(i,j),temporal_corr_before_after(i,j),spatial_corr_before_after(i,j),diagonal_ordering_before_after(i,j)] = event_comparison_6metrics(event1,event2,com1,com2,plot_flag);
    end
end

% before - evoked
disp('2 Calculating pairwise distances: before & evoked...');
f = waitbar(0);
max_cc_before_evoked = nan(length(time_frames_before),length(com_mean_dist_evoked));
com_mean_dist_before_evoked = nan(length(time_frames_before),length(com_mean_dist_evoked));
temporal_corr_before_evoked = nan(length(time_frames_before),length(com_mean_dist_evoked));
spatial_corr_before_evoked = nan(length(time_frames_before),length(com_mean_dist_evoked));
diagonal_ordering_before_evoked = nan(length(time_frames_before),length(com_mean_dist_evoked));
% EMD_before_evoked = nan(length(time_frames_before),length(com_mean_dist_evoked));
for i = 1:length(time_frames_before)
    waitbar(i/length(time_frames_before),f,sprintf([num2str(i),'/',num2str(length(time_frames_before)),' completed']))
    event1 = events_smoothed_before{i}(chosen_electrodes_evoked,:);
    com1 = events_CoM_before{i};
    if max(max(abs(com1)))>0 && max(max(abs(event1)))>1e-5
        parfor j = 1:length(com_mean_dist_evoked)
            event2 = squeeze(events_smoothed_evoked_(j,:,:));
            com2 = squeeze(events_CoM_evoked_(j,:,:));
            if max(max(abs(com2)))>0 && max(max(abs(event2)))>1e-5
%                 [max_cc_before_evoked(i,j),com_mean_dist_before_evoked(i,j),temporal_corr_before_evoked(i,j),spatial_corr_before_evoked(i,j),diagonal_ordering_before_evoked(i,j),EMD_before_evoked(i,j)] = event_comparison_6metrics(event1,event2,com1,com2,plot_flag);
                [max_cc_before_evoked(i,j),com_mean_dist_before_evoked(i,j),temporal_corr_before_evoked(i,j),spatial_corr_before_evoked(i,j),diagonal_ordering_before_evoked(i,j)] = event_comparison_6metrics(event1,event2,com1,com2,plot_flag);
            end
        end
    end
end

% after - evoked
disp('3 Calculating pairwise distances: evoked & after...');
f = waitbar(0);
max_cc_after_evoked = nan(length(time_frames_after),length(com_mean_dist_evoked));
com_mean_dist_after_evoked = nan(length(time_frames_after),length(com_mean_dist_evoked));
temporal_corr_after_evoked = nan(length(time_frames_after),length(com_mean_dist_evoked));
spatial_corr_after_evoked = nan(length(time_frames_after),length(com_mean_dist_evoked));
diagonal_ordering_after_evoked = nan(length(time_frames_after),length(com_mean_dist_evoked));
% EMD_after_evoked = nan(length(time_frames_after),length(com_mean_dist_evoked));
for i = 1:length(time_frames_after)
    waitbar(i/length(time_frames_after),f,sprintf([num2str(i),'/',num2str(length(time_frames_after)),' completed']))
    event1 = events_smoothed_after{i}(chosen_electrodes_evoked,:);
    com1 = events_CoM_after{i};
    if max(max(abs(com1)))>0 && max(max(abs(event1)))>1e-5
        parfor j = 1:length(com_mean_dist_evoked)
            event2 = squeeze(events_smoothed_evoked_(j,:,:));
            com2 = squeeze(events_CoM_evoked_(j,:,:));
            if max(max(abs(com2)))>0 && max(max(abs(event2)))>1e-5
%                 [max_cc_after_evoked(i,j),com_mean_dist_after_evoked(i,j),temporal_corr_after_evoked(i,j),spatial_corr_after_evoked(i,j),diagonal_ordering_after_evoked(i,j),EMD_after_evoked(i,j)] = event_comparison_6metrics(event1,event2,com1,com2,plot_flag);
                [max_cc_after_evoked(i,j),com_mean_dist_after_evoked(i,j),temporal_corr_after_evoked(i,j),spatial_corr_after_evoked(i,j),diagonal_ordering_after_evoked(i,j)] = event_comparison_6metrics(event1,event2,com1,com2,plot_flag);
            end
        end
    end
end

% before - during
disp('4 Calculating pairwise distances: before & during...');
f = waitbar(0);
max_cc_before_during = nan(length(time_frames_before),length(com_mean_dist_during));
com_mean_dist_before_during = nan(length(time_frames_before),length(com_mean_dist_during));
temporal_corr_before_during = nan(length(time_frames_before),length(com_mean_dist_during));
spatial_corr_before_during = nan(length(time_frames_before),length(com_mean_dist_during));
diagonal_ordering_before_during = nan(length(time_frames_before),length(com_mean_dist_during));
% EMD_before_during = nan(length(time_frames_before),length(com_mean_dist_during));
for i = 1:length(time_frames_before)
    waitbar(i/length(time_frames_before),f,sprintf([num2str(i),'/',num2str(length(time_frames_before)),' completed']))
    event1 = events_smoothed_before{i};
    com1 = events_CoM_before{i};
    if max(max(abs(com1)))>0 && max(max(abs(event1)))>1e-5
        parfor j = 1:length(com_mean_dist_during)
            event2 = events_smoothed_during{j};
            com2 = events_CoM_during{j};
            if max(max(abs(com2)))>0 && max(max(abs(event2)))>1e-5
%                 [max_cc_before_during(i,j),com_mean_dist_before_during(i,j),temporal_corr_before_during(i,j),spatial_corr_before_during(i,j),diagonal_ordering_before_during(i,j),EMD_before_during(i,j)] = event_comparison_6metrics(event1,event2,com1,com2,plot_flag);
                [max_cc_before_during(i,j),com_mean_dist_before_during(i,j),temporal_corr_before_during(i,j),spatial_corr_before_during(i,j),diagonal_ordering_before_during(i,j)] = event_comparison_6metrics(event1,event2,com1,com2,plot_flag);
            end
        end
    end
end
delete(f);

% after - during
disp('5 Calculating pairwise distances: after & during...');
f = waitbar(0);
max_cc_after_during = nan(length(time_frames_after),length(com_mean_dist_during));
com_mean_dist_after_during = nan(length(time_frames_after),length(com_mean_dist_during));
temporal_corr_after_during = nan(length(time_frames_after),length(com_mean_dist_during));
spatial_corr_after_during = nan(length(time_frames_after),length(com_mean_dist_during));
diagonal_ordering_after_during = nan(length(time_frames_after),length(com_mean_dist_during));
% EMD_after_during = nan(length(time_frames_after),length(com_mean_dist_during));
for i = 1:length(time_frames_after)
    waitbar(i/length(time_frames_after),f,sprintf([num2str(i),'/',num2str(length(time_frames_after)),' completed']))
    event1 = events_smoothed_after{i};
    com1 = events_CoM_after{i};
    if max(max(abs(com1)))>0 && max(max(abs(event1)))>1e-5
        parfor j = 1:length(com_mean_dist_during)
            event2 = events_smoothed_during{j};
            com2 = events_CoM_during{j};
            if max(max(abs(com2)))>0 && max(max(abs(event2)))>1e-5
%                 [max_cc_after_during(i,j),com_mean_dist_after_during(i,j),temporal_corr_after_during(i,j),spatial_corr_after_during(i,j),diagonal_ordering_after_during(i,j),EMD_after_during(i,j)] = event_comparison_6metrics(event1,event2,com1,com2,plot_flag);
                [max_cc_after_during(i,j),com_mean_dist_after_during(i,j),temporal_corr_after_during(i,j),spatial_corr_after_during(i,j),diagonal_ordering_after_during(i,j)] = event_comparison_6metrics(event1,event2,com1,com2,plot_flag);
            end
        end
    end
end
delete(f);

% evoked - during
disp('6 Calculating pairwise distances: evoked & during...');
f = waitbar(0);
max_cc_evoked_during = nan(length(com_mean_dist_evoked),length(com_mean_dist_during));
com_mean_dist_evoked_during = nan(length(com_mean_dist_evoked),length(com_mean_dist_during));
temporal_corr_evoked_during = nan(length(com_mean_dist_evoked),length(com_mean_dist_during));
spatial_corr_evoked_during = nan(length(com_mean_dist_evoked),length(com_mean_dist_during));
diagonal_ordering_evoked_during = nan(length(com_mean_dist_evoked),length(com_mean_dist_during));
% EMD_evoked_during = nan(length(com_mean_dist_evoked),length(com_mean_dist_during));
for i = 1:length(com_mean_dist_evoked)
    waitbar(i/length(com_mean_dist_evoked),f,sprintf([num2str(i),'/',num2str(length(com_mean_dist_evoked)),' completed']))
    event1 = squeeze(events_smoothed_evoked_(i,:,:));
    com1 = squeeze(events_CoM_evoked_(i,:,:));
    if max(max(abs(com1)))>0 && max(max(abs(event1)))>1e-5
        parfor j = 1:length(com_mean_dist_during)
            event2 = events_smoothed_during{j}(chosen_electrodes_evoked,:);
            com2 = events_CoM_during{j};
            if max(max(abs(com2)))>0 && max(max(abs(event2)))>1e-5
%                 [max_cc_evoked_during(i,j),com_mean_dist_evoked_during(i,j),temporal_corr_evoked_during(i,j),spatial_corr_evoked_during(i,j),diagonal_ordering_evoked_during(i,j),EMD_evoked_during(i,j)] = event_comparison_6metrics(event1,event2,com1,com2,plot_flag);
                [max_cc_evoked_during(i,j),com_mean_dist_evoked_during(i,j),temporal_corr_evoked_during(i,j),spatial_corr_evoked_during(i,j),diagonal_ordering_evoked_during(i,j)] = event_comparison_6metrics(event1,event2,com1,com2,plot_flag);
            end
        end
    end
end
delete(f);

%% save

disp('Saving mat files...');

savefig([path_save,'histograms.fig'])
saveas(gcf,[path_save,'histograms.png'])
clear activity* responses* sliding* slop*

% save all 
if clean_spont_events
    save([path_save,'comparing_events_slidingWin_cleaned_V2.mat']);
    file_name = 'distance_metrix_slidingWin_cleaned_V2.mat';
else
    save([path_save,'comparing_events_slidingWin_V2.mat']);
    file_name = 'distance_metrix_slidingWin_V2.mat';
end

disp('* * * ALL done! :) * * *');
