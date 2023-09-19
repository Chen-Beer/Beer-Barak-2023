%% load data

clc; clear; close all;
exp_list = {'26550_1.11','26549_8.11','38428_17.11','38427_1.11','26532_2.3','new batch\38426_2.11','new batch\39740_2.11','new batch\26549_11.11','new batch\26550_15.11',...
            'control_26550_24.1','control_39740_24.4','control_38427_24.4','control_26536_7.2',...
            'broken_8.11','26532_17.11','38428_20.2','control_38427_7.2'};

exp_idx = 16;

dir_ = ['..\data\',exp_list{exp_idx}];
struct_before = load([dir_,'\1 spont\data.mat'],'activity_smoothed','N','units','dt','sigma','kernel_frame');
struct_after = load([dir_,'\5 spont\data.mat'],'activity_smoothed');

N = struct_before.N;
units = struct_before.units;
dt = struct_before.dt;
sigma = struct_before.sigma;
kernel_frame = struct_before.kernel_frame;
activity_smoothed_before = struct_before.activity_smoothed;
activity_smoothed_after = struct_after.activity_smoothed;
n_patters = 3;

struct_events = load([dir_,'\spont_events_individuals.mat'],'time_frames_before','time_frames_after',...
                                                           'events_smoothed_before','events_smoothed_after',...
                                                           'events_2D_before','events_2D_after',...
                                                           'events_CoM_before','events_CoM_after',...
                                                           'events_CoM_evoked','events_smoothed_evoked');
timeframes_before = struct_events.time_frames_before;
events_before = struct_events.events_smoothed_before;
events2D_before = struct_events.events_2D_before;
events_com_before = struct_events.events_CoM_before;

timeframes_after = struct_events.time_frames_after;
events_after = struct_events.events_smoothed_after;
events2D_after = struct_events.events_2D_after;
events_com_after = struct_events.events_CoM_after;

events_evoked = struct_events.events_smoothed_evoked;
events_com_evoked = struct_events.events_CoM_evoked;

struct_probe1 = load([dir_,'\2 probe\responsiveness_analysis.mat'],'chosen_patterns','responses_smoothed','center_of_mass','stimulus_electrods_idx','stimulus_electrods');
chosen_patterns = struct_probe1.chosen_patterns;
probe_responses_smoothed = struct_probe1.responses_smoothed;
probe_com = struct_probe1.center_of_mass;
probe_stimulus_electrods_idx = struct_probe1.stimulus_electrods_idx;

clear struct_before struct_after struct_events struct_during struct_probe1

n_events_before = length(events_com_before);
n_events_after = length(events_com_after);
n_patterns = size(events_evoked,1);

T = 30;

%% for each event - extract the IC

event_ti = 1:1; % temporal index for IC

ICs_all_before = zeros(n_events_before,N*length(event_ti));
for e = 1:n_events_before
    tmp = events_before{e}(:,event_ti)';
    ICs_all_before(e,:) = tmp(:);
end

ICs_all_after = zeros(n_events_after,N*length(event_ti));
for e = 1:n_events_after
    tmp = events_after{e}(:,event_ti)';
    ICs_all_after(e,:) = tmp(:);
end

% calc distance between ICs - before & after
dist_ICs_after_before = zeros(n_events_after,n_events_before);
corr_ICs_after_before = zeros(n_events_after,n_events_before);
parfor i = 1:n_events_after
    for j = 1:n_events_before
%         dist_ICs_after_before(i,j) = norm(ICs_all_after(i,:) - ICs_all_before(j,:));
        dist_ICs_after_before(i,j) = norm(ICs_all_after(i,:) - ICs_all_before(j,:))/(norm(ICs_all_after(i,:))*norm(ICs_all_before(j,:)));
        corr_ICs_after_before(i,j) = corr(ICs_all_after(i,:)',ICs_all_before(j,:)');
    end
end

% calc distance between ICs - before & before
dist_ICs_before_before = zeros(n_events_before,n_events_before);
corr_ICs_before_before = zeros(n_events_before,n_events_before);
parfor i = 1:n_events_before
    for j = 1:n_events_before
        dist_ICs_before_before(i,j) = norm(ICs_all_before(i,:) - ICs_all_before(j,:))/(norm(ICs_all_before(i,:))*norm(ICs_all_before(j,:)));
        corr_ICs_before_before(i,j) = corr(ICs_all_before(i,:)',ICs_all_before(j,:)');
    end
end

%% intial states throughout stimulation

struct_during = load([dir_,'\freshlook_merged_data.mat'],'n_minutes','events_during_minutes');
n_minutes = struct_during.n_minutes;
events_during_minutes = struct_during.events_during_minutes;
win_size_minutes = 30;
n_windows_during = floor(n_minutes/win_size_minutes);

clear struct_during

% divide during events to windows
events_during_windowed = cell(n_windows_during,1);
widx = 1;
for w = 1:n_windows_during
    subset_ = events_during_minutes(widx:widx+win_size_minutes-1);
    events_ = cat(1, subset_{:});  
    events_during_windowed{w} = events_;
    widx = widx + win_size_minutes;
end

% extract ICs
ICs_during_windowed = cell(n_windows_during,1);
parfor w = 1:n_windows_during
    ICs_during_windowed{w} = zeros(length(events_during_windowed{w}),N*length(event_ti));
    for e = 1:length(events_during_windowed{w})
        ICs_during_windowed{w}(e,:) = events_during_windowed{w}{e}(:,event_ti);
    end
end

%% 1D IC predicting 2D burst? - during

% calc corrs - during-during
corrs_events_during_during = cell(n_windows_during,1);
dist_events_during_during = cell(n_windows_during,1);
corrs_ICs_during_during = cell(n_windows_during,1);
dist_ICs_during_during = cell(n_windows_during,1);
parfor w = 1:n_windows_during
    s = length(events_during_windowed{w});
    corrs_events_during_during{w} = zeros(s,s);
    dist_events_during_during{w} = zeros(s,s);
    corrs_ICs_during_during{w} = zeros(s,s);    
    dist_ICs_during_during{w} = zeros(s,s);    
    for i = 1:s
    %     disp(i)
        if size(events_during_windowed{w}{i},2) >= T
            e1 = events_during_windowed{w}{i}(:,1:T);
            for j = 1:s
                if size(events_during_windowed{w}{j},2) >= T
                    e2 = events_during_windowed{w}{j}(:,1:T);
                    corrs_events_during_during{w}(i,j) = corr2(e1,e2);
                    dist_events_during_during{w}(i,j) = norm(e1-e2);
                    corrs_ICs_during_during{w}(i,j) = corr(ICs_during_windowed{w}(i,:)',ICs_during_windowed{w}(j,:)');
                    dist_ICs_during_during{w}(i,j) = norm(ICs_during_windowed{w}(i,:)-ICs_during_windowed{w}(j,:))/(norm(ICs_during_windowed{w}(i,:))*norm(ICs_during_windowed{w}(j,:)));
                end
            end
        end
    end
end

% calc corrs - before-during
corrs_events_before_during = cell(n_windows_during,1);
dist_events_before_during = cell(n_windows_during,1);
corrs_ICs_before_during = cell(n_windows_during,1);
dist_ICs_before_during = cell(n_windows_during,1);
parfor w = 1:n_windows_during
    s = length(events_during_windowed{w});
    corrs_events_before_during{w} = zeros(n_events_before,s);
    dist_events_before_during{w} = zeros(n_events_before,s);
    corrs_ICs_before_during{w} = zeros(n_events_before,s);
    dist_ICs_before_during{w} = zeros(n_events_before,s);
    for i = 1:n_events_before
    %     disp(i)
        if size(events_before{i},2) >= T
            e1 = events_before{i}(:,1:T);
            for j = 1:s
                if size(events_during_windowed{w}{j},2) >= T
                    e2 = events_during_windowed{w}{j}(:,1:T);
                    corrs_events_before_during{w}(i,j) = corr2(e1,e2);
                    dist_events_before_during{w}(i,j) = norm(e1-e2);
                    corrs_ICs_before_during{w}(i,j) = corr(ICs_all_before(i,:)',ICs_during_windowed{w}(j,:)');
                    dist_ICs_before_during{w}(i,j) = norm(ICs_all_before(i,:)-ICs_during_windowed{w}(j,:))/(norm(ICs_all_before(i,:))*norm(ICs_during_windowed{w}(j,:)));
                end
            end
        end
    end
end

%% 1D IC predicting 2D burst? - before

% calc corrs
corrs_events_before_before = zeros(n_events_before,n_events_before);
dist_events_before_before = zeros(n_events_before,n_events_before);
corrs_ICs_before_before = zeros(n_events_before,n_events_before);
dist_ICs_before_before = zeros(n_events_before,n_events_before);
parfor i = 1:n_events_before
%     disp(i)
    if size(events_before{i},2) >= T
        e1 = events_before{i}(:,1:T);
        for j = 1:n_events_before
            if size(events_before{j},2) >= T
                e2 = events_before{j}(:,1:T);
                corrs_events_before_before(i,j) = corr2(e1,e2);
                dist_events_before_before(i,j) = norm(e1-e2);
                corrs_ICs_before_before(i,j) = corr(ICs_all_before(i,:)',ICs_all_before(j,:)');
%                 dist_ICs_before_before(i,j) = norm(ICs_all_before(i,:)-ICs_all_before(j,:));
                dist_ICs_before_before(i,j) = norm(ICs_all_before(i,:)-ICs_all_before(j,:))/(norm(ICs_all_before(i,:))*norm(ICs_all_before(j,:)));
            end
        end
    end
end

nbins = 50*[1 1];

figure; 
subplot 121;
tmp1 = triu(corrs_events_before_before,1); tmp1 = tmp1(:);
tmp2 = triu(dist_ICs_before_before,1); tmp2 = tmp2(:); 
hist3([tmp2(tmp2 ~= 0),tmp1(tmp1 ~= 0)],'LineStyle','none','Nbins',nbins,'CdataMode','auto'); colorbar; view(2);
ylabel('burst_{2D corr}'); xlabel('IC_{distance}');
subplot 122;
tmp1 = triu(dist_events_before_before,1); tmp1 = tmp1(:);
tmp2 = triu(dist_ICs_before_before,1); tmp2 = tmp2(:); 
hist3([tmp1(tmp1 ~= 0),tmp2(tmp2 ~= 0)],'LineStyle','none','Nbins',nbins,'CdataMode','auto'); colorbar; view(2);
xlabel('burst_{distance}'); ylabel('IC_{distance}');
sgtitle('Before')

figure; 
tmp1 = triu(corrs_events_before_before,1); tmp1 = tmp1(:);
tmp2 = triu(corr_ICs_before_before,1); tmp2 = tmp2(:); 
hist3([tmp2(tmp2 ~= 0),tmp1(tmp1 ~= 0)],'LineStyle','none','Nbins',nbins,'CdataMode','auto'); colorbar; view(2);
histogram2(tmp2(tmp2 ~= 0),tmp1(tmp1 ~= 0),'Normalization','pdf','FaceColor','flat','LineStyle','none'); colormap;
ylabel('burst_{2D corr}'); xlabel('IC_{corr}');

%% 1D IC predicting 2D burst? - after

% calc corrs
corrs_events_after_after = zeros(n_events_after,n_events_after);
dist_events_after_after = zeros(n_events_after,n_events_after);
corrs_ICs_after_after = zeros(n_events_after,n_events_after);
dist_ICs_after_after = zeros(n_events_after,n_events_after);
parfor i = 1:n_events_after
%     disp(i)
    if size(events_after{i},2) >= T
        e1 = events_after{i}(:,1:T);
        for j = 1:n_events_after
            if size(events_after{j},2) >= T
                e2 = events_after{j}(:,1:T);
                corrs_events_after_after(i,j) = corr2(e1,e2);
                dist_events_after_after(i,j) = norm(e1-e2);
                corrs_ICs_after_after(i,j) = corr(ICs_all_after(i,:)',ICs_all_after(j,:)');
                dist_ICs_after_after(i,j) = norm(ICs_all_after(i,:)-ICs_all_after(j,:))/(norm(ICs_all_after(i,:))*norm(ICs_all_after(j,:)));
            end
        end
    end
end

figure; 
subplot 121;
tmp1 = triu(corrs_events_after_after,1); tmp1 = tmp1(:);
tmp2 = triu(dist_ICs_after_after,1); tmp2 = tmp2(:); 
hist3([tmp1(tmp1 ~= 0),tmp2(tmp2 ~= 0)],'LineStyle','none','Nbins',nbins,'CdataMode','auto'); colorbar; view(2);
xlabel('burst_{2D corr}'); ylabel('IC_{distance}');
subplot 122;
tmp1 = triu(dist_events_after_after,1); tmp1 = tmp1(:);
tmp2 = triu(dist_ICs_after_after,1); tmp2 = tmp2(:); 
hist3([tmp1(tmp1 ~= 0),tmp2(tmp2 ~= 0)],'LineStyle','none','Nbins',nbins,'CdataMode','auto'); colorbar; view(2);
xlabel('burst_{distance}'); ylabel('IC_{distance}');
sgtitle('After');

%% scatter plot - IC predicting burst? cross before-after

% calc corrs
corrs_events_before_after = zeros(n_events_before,n_events_after);
dist_events_before_after = zeros(n_events_before,n_events_after);
corrs_ICs_before_after = zeros(n_events_before,n_events_after);
dist_ICs_before_after = zeros(n_events_before,n_events_after);
parfor i = 1:n_events_before
    disp(i)
    if size(events_before{i},2) >= T
        e1 = events_before{i}(:,1:T);
        for j = 1:n_events_after
            if size(events_after{j},2) >= T
                e2 = events_after{j}(:,1:T);
                corrs_events_before_after(i,j) = corr2(e1,e2);
                dist_events_before_after(i,j) = norm(e1-e2);
                corrs_ICs_before_after(i,j) = corr(ICs_all_before(i,:)',ICs_all_after(j,:)');
                dist_ICs_before_after(i,j) = norm(ICs_all_before(i,:)-ICs_all_after(j,:))/(norm(ICs_all_before(i,:))*norm(ICs_all_after(j,:)));
            end
        end
    end
end

figure; 
subplot 121;
hist3([corrs_events_before_after(:),dist_ICs_before_after(:)],'LineStyle','none','Nbins',nbins,'CdataMode','auto'); colorbar; view(2);
xlabel('burst_{2D corr}'); ylabel('IC_{distance}');
subplot 122;
hist3([dist_events_before_after(:),dist_ICs_before_after(:)],'LineStyle','none','Nbins',nbins,'CdataMode','auto'); colorbar; view(2);
xlabel('burst_{distance}'); ylabel('IC_{distance}');
sgtitle('Before - After');

%% for the chosen probes - find similar spontaneous events in before 

% calc corrs
all_electrodes = 0;
Np = length(probe_stimulus_electrods_idx);
R = size(probe_responses_smoothed,2);
all_e = 1:N;
corrs_probe_before = zeros(n_patterns,R,n_events_before);
for p = 1:n_patterns
%     disp(p)
    probe_ = squeeze(probe_responses_smoothed(chosen_patterns(p),:,:,:));
    ei = 1:N;
    if ~all_electrodes
        stim_e = probe_stimulus_electrods_idx{chosen_patterns(p)};
        ei = all_e(~ismember(all_e,stim_e));
        probe_ = probe_(:,ei,:);
    end
    for i = 1:R
        e1 = squeeze(probe_(i,:,1:T));
        for e = 1:n_events_before
            if size(events_before{e},2) >= T
                e2 = events_before{e}(ei,1:T);
                corrs_probe_before(p,i,e) = corr2(e1,e2);
            end
        end
    end
end

% find similar event for each probe 
similar_events_probe_before_idx = cell(n_patterns,1);
corr_th = 0.73;
median_corr = zeros(n_patterns,n_events_before);
for p = 1:n_patterns
    idx = [];
    for e = 1:n_events_before
        tmp = squeeze(corrs_probe_before(p,:,e)); tmp = tmp(:);
        median_corr(p,e) = median(tmp);
        if median(tmp) > corr_th
            idx = [idx e];
        end
    end
    similar_events_probe_before_idx{p} = idx;
end

figure; 
imagesc(median_corr); colorbar;
title('median corr(probe,before)');
xlabel('Events'); ylabel('Chosen');

%% find similar ICs & 2D events in before & after (based on similar_events_probe_before_idx)

event_corr_th = 0.8;
IC_dist_th = 1.2;

figure; 
subplot 121; hold on; 
histogram(dist_ICs_before_before,'Normalization','Probability')
histogram(dist_ICs_after_before,'Normalization','Probability')
histogram(dist_ICs_after_after,'Normalization','Probability')
plot([IC_dist_th IC_dist_th],[0 0.02],'--k','linewidth',2);
title('Distance between ICs');
legend('before-before','after-before','after-after');
subplot 122; hold on; 
histogram(corrs_events_before_before,'Normalization','Probability')
histogram(corrs_events_before_after,'Normalization','Probability')
histogram(corrs_events_after_after,'Normalization','Probability')
plot([event_corr_th event_corr_th],[0 0.02],'--k','linewidth',2);
title('2D corr between events');
legend('before-before','after-before','after-after');

figure; 
subplot 121; hold on; 
for w = 1:n_windows_during
    histogram(dist_ICs_before_during{w},'Normalization','Probability')
end
plot([IC_dist_th IC_dist_th],[0 0.02],'--k','linewidth',2);
title('Distance between ICs - before\during');
subplot 122; hold on; 
for w = 1:n_windows_during
    histogram(corrs_events_before_during{w},'Normalization','Probability')
end
plot([event_corr_th event_corr_th],[0 0.02],'--k','linewidth',2);
title('2D corr between events - before\during');

% after - events
similar_events_after_idx = cell(n_patterns,1);
for i = 1:n_events_after
    for p = 1:n_patterns
        if ~isempty(similar_events_probe_before_idx{p})
            idx_ = similar_events_probe_before_idx{p};
            if median(corrs_events_before_after(idx_,i)) >= event_corr_th
                similar_events_after_idx{p} = [similar_events_after_idx{p} i];
            end
        end
    end
end
% after - ICs
similar_ICs_after_idx = cell(n_patterns,1);
for i = 1:n_events_after
    for p = 1:n_patterns
        if ~isempty(similar_events_probe_before_idx{p})
            idx_ = similar_events_probe_before_idx{p};
            if median(dist_ICs_after_before(i,idx_)) <= IC_dist_th
                similar_ICs_after_idx{p} = [similar_ICs_after_idx{p} i];
            end
        end
    end
end

% before - events
similar_events_before_idx = cell(n_patterns,1);
for i = 1:n_events_before
    for p = 1:n_patterns
        if ~isempty(similar_events_probe_before_idx{p})
            idx_ = similar_events_probe_before_idx{p};
            if median(corrs_events_before_before(i,idx_)) >= event_corr_th
                similar_events_before_idx{p} = [similar_events_before_idx{p} i];
            end
        end
    end
end
% before - ICs
similar_ICs_before_idx = cell(n_patterns,1);
for i = 1:n_events_before
    for p = 1:n_patterns
        if ~isempty(similar_events_probe_before_idx{p})
            idx_ = similar_events_probe_before_idx{p};
            if median(dist_ICs_before_before(i,idx_)) <= IC_dist_th
                similar_ICs_before_idx{p} = [similar_ICs_before_idx{p} i];
            end
        end
    end
end

% during - events
similar_events_during_idx = cell(n_patterns,n_windows_during);
for w = 1:n_windows_during
    for i = 1:length(events_during_windowed{w})
        for p = 1:n_patterns
            if ~isempty(similar_events_probe_before_idx{p})
                idx_ = similar_events_probe_before_idx{p};
                if median(corrs_events_before_during{w}(idx_,i)) >= event_corr_th
                    similar_events_during_idx{p,w} = [similar_events_during_idx{p,w} i];
                end
            end
        end
    end
end
% during - ICs
similar_ICs_during_idx = cell(n_patterns,n_windows_during);
for w = 1:n_windows_during
    for i = 1:length(events_during_windowed{w})
        for p = 1:n_patterns
            if ~isempty(similar_events_probe_before_idx{p})
                idx_ = similar_events_probe_before_idx{p};
                if median(dist_ICs_before_during{w}(idx_,i)) <= IC_dist_th
                    similar_ICs_during_idx{p,w} = [similar_ICs_during_idx{p,w} i];
                end
            end
        end
    end
end

%% convergence plots - similarity through time - events similar to chosen - before & after
% for each chosen pattern plot convergence plot for all its similar events

corr_th2 = event_corr_th;
figure; 
sgtitle(['Correlation & distance through time - similar events | corr > ',num2str(corr_th2)]);
for p = 1:n_patterns
    % before
    idx = similar_events_before_idx{p};
    if ~isempty(idx)
        events_ = events_before(idx);
        corrs = [];
        dists = [];
        for i = 1:length(idx)
            e1 = events_{i}(:,1:T);        
            for j = i+1 : length(idx)
                e2 = events_{j}(:,1:T);        
                if corr2(e1,e2) > corr_th2
                    corr_ = zeros(1,T);
                    dist_ = zeros(1,T);
                    for t = 1:T
                        corr_(t) = corr(e1(:,t),e2(:,t));
                        dist_(t) = norm(e1(:,t)-e2(:,t))/(norm(e1(:,t))*norm(e2(:,t)));
                    end
                    corrs = [corrs; corr_];
                    dists = [dists; dist_];
                end
            end
        end
        subplot (2,n_patterns,p); hold on; title(p);    
        shadedErrorBar(dt*(1:T),mean(corrs),std(corrs));
        xlabel('Time [sec]'); ylim([0 1]); ylabel('corr');
        drawnow;
        subplot (2,n_patterns,n_patterns+p); hold on; 
        shadedErrorBar(dt*(1:T),mean(dists),std(dists));
        xlabel('Time [sec]'); ylabel('dist');
        drawnow; 
    end
    % after
    idx = similar_events_after_idx{p};
    if ~isempty(idx)
        events_ = events_after(idx);
        corrs = [];
        dists = [];
        for i = 1:length(idx)
            e1 = events_{i}(:,1:T);        
            for j = i+1 : length(idx)
                e2 = events_{j}(:,1:T);        
                if corr2(e1,e2) > corr_th2
                    corr_ = zeros(1,T);
                    dist_ = zeros(1,T);
                    for t = 1:T
                        corr_(t) = corr(e1(:,t),e2(:,t));
                        dist_(t) = norm(e1(:,t)-e2(:,t))/(norm(e1(:,t))*norm(e2(:,t)));
                    end
                    corrs = [corrs; corr_];
                    dists = [dists; dist_];
                end
            end
        end
        subplot (2,n_patterns,p);
        shadedErrorBar(dt*(1:T),mean(corrs),std(corrs),'lineProps','-b');
        xlabel('Time [sec]'); ylim([0 1]); ylabel('corr');
        drawnow;
        subplot (2,n_patterns,n_patterns+p);
        shadedErrorBar(dt*(1:T),mean(dists),std(dists),'lineProps','-b');
        xlabel('Time [sec]'); ylabel('dist');
        drawnow; 
    end
end

%% convergence plots - similarity through time - ICs similar to chosen - before & after
% for each chosen pattern plot convergence plot for all its similar ICs

dist_th2 = IC_dist_th;
color_vec = summer(n_windows_during);
figure; 
sgtitle(['Correlation & distance through time - similar ICs | dist < ',num2str(dist_th2)]);
for p = 1:n_patterns
    % before
    idx = similar_ICs_before_idx{p};
    if ~isempty(idx)
        events_ = events_before(idx);
        ICs_ = ICs_all_before(idx,:);
        corrs = [];
        dists = [];
        for i = 1:length(idx)
            if size(events_{i},2) >= T
                e1 = events_{i}(:,1:T);        
                for j = i+1 : length(idx)
                    if size(events_{j},2) >= T
                        e2 = events_{j}(:,1:T);        
                        if norm(ICs_(i,:)-ICs_(j,:))/(norm(ICs_(i,:))*norm(ICs_(j,:))) < dist_th2
                            corr_ = zeros(1,T);
                            dist_ = zeros(1,T);
                            for t = 1:T
                                corr_(t) = corr(e1(:,t),e2(:,t));
                                dist_(t) = norm(e1(:,t)-e2(:,t))/(norm(e1(:,t))*norm(e2(:,t)));
                            end
                            corrs = [corrs; corr_];
                            dists = [dists; dist_];
                        end
                    end
                end
            end
        end
        if ~isempty(corrs)
            subplot (2,n_patterns,p); hold on; title(p);    
            shadedErrorBar(dt*(1:T),mean(corrs),std(corrs));
            xlabel('Time [sec]'); ylim([0 1]); ylabel('corr');
            drawnow;
            subplot (2,n_patterns,n_patterns+p); hold on; 
            shadedErrorBar(dt*(1:T),mean(dists),std(dists));
            xlabel('Time [sec]'); ylabel('dist');
            drawnow; 
        end
    end
    % after
    idx = similar_ICs_after_idx{p};
    if ~isempty(idx)
        events_ = events_after(idx);
        ICs_ = ICs_all_after(idx,:);        
        corrs = [];
        dists = [];
        for i = 1:length(idx)
            if size(events_{i},2) >= T
                e1 = events_{i}(:,1:T);        
                for j = i+1 : length(idx)
                    if size(events_{j},2) >= T
                        e2 = events_{j}(:,1:T);        
                        if norm(ICs_(i,:)-ICs_(j,:))/(norm(ICs_(i,:))*norm(ICs_(j,:))) < dist_th2
                            corr_ = zeros(1,T);
                            dist_ = zeros(1,T);
                            for t = 1:T
                                corr_(t) = corr(e1(:,t),e2(:,t));
                                dist_(t) = norm(e1(:,t)-e2(:,t))/(norm(e1(:,t))*norm(e2(:,t)));
                            end
                            corrs = [corrs; corr_];
                            dists = [dists; dist_];
                        end
                    end
                end
            end
        end
        if ~isempty(corrs)
            subplot (2,n_patterns,p);  
            shadedErrorBar(dt*(1:T),mean(corrs),std(corrs),'lineProps','-b');
            xlabel('Time [sec]'); ylim([0 1]); ylabel('corr');
            drawnow;
            subplot (2,n_patterns,n_patterns+p); 
            shadedErrorBar(dt*(1:T),mean(dists),std(dists),'lineProps','-b');
            xlabel('Time [sec]'); ylabel('dist');
            drawnow; 
        end
    end
    % during
    for w = 1:n_windows_during
        idx = similar_ICs_during_idx{p,w};
        if ~isempty(idx)
            events_ = events_during_windowed{w}(idx);
            ICs_ = ICs_during_windowed{w}(idx,:);        
            corrs = [];
            dists = [];
            for i = 1:length(idx)
                if size(events_{i},2) >= T
                    e1 = events_{i}(:,1:T);        
                    for j = i+1 : length(idx)
                        if size(events_{j},2) >= T
                            e2 = events_{j}(:,1:T);        
                            if norm(ICs_(i,:)-ICs_(j,:))/(norm(ICs_(i,:))*norm(ICs_(j,:))) < dist_th2
                                corr_ = zeros(1,T);
                                dist_ = zeros(1,T);
                                for t = 1:T
                                    corr_(t) = corr(e1(:,t),e2(:,t));
                                    dist_(t) = norm(e1(:,t)-e2(:,t))/(norm(e1(:,t))*norm(e2(:,t)));
                                end
                                corrs = [corrs; corr_];
                                dists = [dists; dist_];
                            end
                        end
                    end
                end
            end
            if ~isempty(corrs)
                subplot (2,n_patterns,p);  hold on; 
                plot(dt*(1:T),mean(corrs),'color',color_vec(w,:));
%                 shadedErrorBar(dt*(1:T),mean(corrs),std(corrs),'lineprops',{'color',color_vec(w,:)});
                xlabel('Time [sec]'); ylim([0 1]); ylabel('corr');
                drawnow;
                subplot (2,n_patterns,n_patterns+p); hold on; 
                plot(dt*(1:T),mean(dists),'color',color_vec(w,:));
%                 shadedErrorBar(dt*(1:T),mean(dists),std(dists),'lineprops',{'color',color_vec(w,:)});
                xlabel('Time [sec]'); ylabel('dist');
                drawnow; 
            end
        end 
    end
end

%% convergence plots - similarity through time - evoked (chosen)
% for each chosen pattern plot convergence plot for all its repetitions

figure; 
sgtitle('Correlation & distance through time');
for p = 1:n_patterns
    % evoked responses
    events_ = squeeze(events_evoked(p,end-20:end,:,:));
    corrs = [];
    dists = [];
    for i = 1:size(events_,1)
        e1 = squeeze(events_(i,:,1:T));        
        for j = i+1 : size(events_,1)
            e2 = squeeze(events_(j,:,1:T));       
            corr_ = zeros(1,T);
            dist_ = zeros(1,T);
            for t = 1:T
                corr_(t) = corr(e1(:,t),e2(:,t));
                dist_(t) = norm(e1(:,t)-e2(:,t))/(norm(e1(:,t))*norm(e2(:,t)));
            end
            corrs = [corrs; corr_];
            dists = [dists; dist_];
        end
    end
    subplot (2,n_patterns,p); hold on; title(p);    
    shadedErrorBar(dt*(1:T),mean(corrs),std(corrs));
    xlabel('Time [sec]'); ylim([0 1]); ylabel('corr');
    drawnow;
    subplot (2,n_patterns,n_patterns+p); hold on; 
    shadedErrorBar(dt*(1:T),mean(dists),std(dists));
    xlabel('Time [sec]'); ylabel('dist');
    drawnow; 
end

%% Plot histograms: pairs of close ICs

corrs_before_similar_ = [];
corrs_after_similar_ = [];
corrs_during_similar_ = cell(n_windows_during,1);
corrs_before_NOTsimilar_ = [];
corrs_after_NOTsimilar_ = [];
corrs_during_NOTsimilar_ = cell(n_windows_during,1);

for i = 1:n_events_before
    if ismember(i,cat(2,similar_ICs_before_idx{:})) % similar IC to chosen
        for j = i+1:n_events_before
            if dist_ICs_before_before(i,j) <= IC_dist_th && ismember(j,cat(2,similar_ICs_before_idx{:}))
                corrs_before_similar_ = [corrs_before_similar_ corrs_events_before_before(i,j)];
            end
        end
    else % IC not dimilar to chosen
        if ~ismember(i,cat(2,similar_events_before_idx{:})) % burst also different from chosen
            for j = i+1:n_events_before
                if dist_ICs_before_before(i,j) <= IC_dist_th && ~ismember(j,cat(2,similar_events_before_idx{:}))
                    corrs_before_NOTsimilar_ = [corrs_before_NOTsimilar_ corrs_events_before_before(i,j)];
                end
            end        
        end
    end
end
for i = 1:n_events_after
    if ismember(i,cat(2,similar_ICs_after_idx{:}))    
        for j = i+1:n_events_after
            if dist_ICs_after_after(i,j) <= IC_dist_th && ismember(j,cat(2,similar_ICs_after_idx{:}))    
                corrs_after_similar_ = [corrs_after_similar_ corrs_events_after_after(i,j)];
            end
        end
    else
        if ~ismember(i,cat(2,similar_events_after_idx{:}))
            for j = i+1:n_events_after
                if dist_ICs_after_after(i,j) <= IC_dist_th && ~ismember(j,cat(2,similar_events_after_idx{:}))   
                    corrs_after_NOTsimilar_ = [corrs_after_NOTsimilar_ corrs_events_after_after(i,j)];
                end
            end        
        end
    end
end
parfor w = 1:n_windows_during
    corrs_during_similar_{w} = [];
    corrs_during_NOTsimilar_{w} = [];
    for i = 1:length(events_during_windowed{w})
        if ismember(i,cat(2,similar_ICs_during_idx{:,w}))    
            for j = i+1:length(events_during_windowed{w})
                if dist_ICs_during_during{w}(i,j) <= IC_dist_th && ismember(j,cat(2,similar_ICs_during_idx{:,w}))    
                    corrs_during_similar_{w} = [corrs_during_similar_{w} corrs_events_during_during{w}(i,j)];
                end
            end
        else
            if ~ismember(i,cat(2,similar_ICs_during_idx{:,w}))
                for j = i+1:length(events_during_windowed{w})
                    if dist_ICs_during_during{w}(i,j) <= IC_dist_th && ~ismember(j,cat(2,similar_events_during_idx{:,w}))   
                        corrs_during_NOTsimilar_{w} = [corrs_during_NOTsimilar_{w} corrs_events_during_during{w}(i,j)];
                    end
                end        
            end
        end
    end
end

figure; 
subplot 121; hold on;
histogram(corrs_before_similar_,50,'Normalization','pdf')
histogram(corrs_after_similar_,50,'Normalization','pdf')
% for w = 1:n_windows_during
%     histogram(corrs_during_similar_{w},'Normalization','Probability','facecolor',color_vec(w,:))
% end
xlabel('2D corr(event_1,event_2)');
legend('before','after');
title(['close ICs | dist(IC_1,IC_2) < ',num2str(IC_dist_th),' | ICS similar to chosen']);
subplot 122; hold on;
histogram(corrs_before_NOTsimilar_,50,'Normalization','pdf')
histogram(corrs_after_NOTsimilar_,50,'Normalization','pdf')
% for w = 1:n_windows_during
%     histogram(corrs_during_NOTsimilar_{w},'Normalization','Probability','facecolor',color_vec(w,:))
% end
xlabel('2D corr(event_1,event_2)');
legend('before','after');
title(['close ICs | dist(IC_1,IC_2) < ',num2str(IC_dist_th),' | ICS & events - not similar to chosen']);

figure;
subplot 121; hold on;
plot(1,median(corrs_before_similar_),'ok','markersize',5,'linewidth',1.5)
plot(1+n_windows_during,median(corrs_after_similar_),'xk','markersize',5,'linewidth',1.5)
for w = 1:n_windows_during
    plot(1+w,median(corrs_during_similar_{w}),'.','color',color_vec(w,:),'markersize',15)
end
ylabel('2D corr(event_1,event_2)'); ylim([0.4 1]);
title('ICs similar to chosen');
subplot 122; hold on;
plot(1,median(corrs_before_NOTsimilar_),'ok','markersize',5,'linewidth',1.5)
plot(1+n_windows_during,median(corrs_after_NOTsimilar_),'xk','markersize',5,'linewidth',1.5)
for w = 1:n_windows_during
    plot(1+w,median(corrs_during_NOTsimilar_{w}),'.','color',color_vec(w,:),'markersize',15)
end
ylabel('2D corr(event_1,event_2)'); ylim([0.4 1]);
title('ICs & events not similar to chosen');
sgtitle(['close ICs | dist(IC_1,IC_2) < ',num2str(IC_dist_th),' | median correlation - before, after & during stimulation']);

%% heatmaps - chosen vs. others

% similar to chosen
figure; 
subplot 121; 
tmp1_ = [];
tmp2_ = [];
for p = 1:n_patterns
    tmp1 = corrs_events_before_before(similar_ICs_before_idx{p},similar_ICs_before_idx{p});
    tmp2 = dist_ICs_before_before(similar_ICs_before_idx{p},similar_ICs_before_idx{p});
    tmp1_ = [tmp1_; tmp1(:)];
    tmp2_ = [tmp2_; tmp2(:)];
end
hist3([tmp1_';tmp2_']','LineStyle','none','Nbins',nbins,'CdataMode','auto'); colorbar; view(2);
title('Before');
xlabel('burst_{2D corr}'); ylabel('IC_{distance}');
subplot 122; 
tmp1_ = [];
tmp2_ = [];
for p = 1:n_patterns
    tmp1 = corrs_events_after_after(similar_ICs_after_idx{p},similar_ICs_after_idx{p});
    tmp2 = dist_ICs_after_after(similar_ICs_after_idx{p},similar_ICs_after_idx{p});
    tmp1_ = [tmp1_; tmp1(:)];
    tmp2_ = [tmp2_; tmp2(:)];
end
hist3([tmp1_';tmp2_']','LineStyle','none','Nbins',nbins,'CdataMode','auto'); colorbar; view(2);
xlabel('burst_{2D corr}'); ylabel('IC_{distance}');
title('After');
sgtitle('Similar ICs');

% all others
tmp_vec = 1:n_events_before;
before_others = tmp_vec(~ismember(1:n_events_before,cat(2,similar_ICs_before_idx{:})));
before_others = before_others(~ismember(before_others,cat(2,similar_events_before_idx{:})));
tmp_vec = 1:n_events_after;
after_others = tmp_vec(~ismember(1:n_events_after,cat(2,similar_ICs_after_idx{:})));
after_others = after_others(~ismember(after_others,cat(2,similar_events_after_idx{:})));
figure; 
subplot 121; 
tmp1 = corrs_events_before_before(before_others,before_others);
tmp2 = dist_ICs_before_before(before_others,before_others);
hist3([tmp1(:)';tmp2(:)']','LineStyle','none','Nbins',nbins,'CdataMode','auto'); colorbar; view(2);
title('Before');
xlabel('burst_{2D corr}'); ylabel('IC_{distance}');
subplot 122; 
tmp1 = corrs_events_after_after(after_others,after_others);
tmp2 = dist_ICs_after_after(after_others,after_others);
hist3([tmp1(:)';tmp2(:)']','LineStyle','none','Nbins',nbins,'CdataMode','auto'); colorbar; view(2);
xlabel('burst_{2D corr}'); ylabel('IC_{distance}');
title('After');
sgtitle('Not similar ICs & bursts');

%% same IC before & after
a = 10; b = 12;

figure; 
for p = 1:n_patterns
    after_idx = similar_ICs_after_idx{p};
    before_idx = similar_ICs_before_idx{p};
    subplot(1,n_patterns,p); hold on;
    avg_after = zeros(2,T);
    median_after = zeros(2,T);
    for i = after_idx
        tmp = events_com_after{i};
        if size(tmp,2) >= T
            avg_after = avg_after + tmp(:,1:T);
            median_after = cat(3,median_after,tmp(:,1:T));
%             plot(tmp(1,1:T),tmp(2,1:T),'-b');
%             plot(tmp(1,1),tmp(2,1),'ob');
        end
    end    
    avg_after = avg_after./length(after_idx); 
    median_after = squeeze(median(median_after,3));
    plot(avg_after(1,1:T),avg_after(2,1:T),'-b');
    plot(avg_after(1,1),avg_after(2,1),'ob');    
    avg_before = zeros(2,T);
    median_before = zeros(2,T);
    for i = before_idx
        tmp = events_com_before{i};
        if size(tmp,2) >= T
            avg_before = avg_before + tmp(:,1:T);
            median_before = cat(3,median_before,tmp(:,1:T));
%             plot(tmp(1,1:T),tmp(2,1:T),'-k');
%             plot(tmp(1,1),tmp(2,1),'ok');
        end
    end
    avg_before = avg_before./length(before_idx);
    median_before = squeeze(median(median_before,3));
    plot(avg_before(1,1:T),avg_before(2,1:T),'-k');
    plot(avg_before(1,1),avg_before(2,1),'ok');  
    xlim([1 b]); ylim([1 a]); grid on;
    title(p);
    drawnow;
end

%% ICs - 2D coverage - before vs. after

% center of mass
ICs_com_before = zeros(n_events_before,2);
for i = 1:n_events_before
    ICs_com_before(i,:) =  events_com_before{i}(:,1);
end
ICs_com_after = zeros(n_events_after,2);
for i = 1:n_events_after
    ICs_com_after(i,:) =  events_com_after{i}(:,1);
end

% plot
figure; hold on; 
scatter(ICs_com_before(:,1),ICs_com_before(:,2),2,'filled')
scatter(ICs_com_after(:,1),ICs_com_after(:,2),2,'filled')
grid on; 
% figure;
% subplot 121;
% hist3([ICs_com_before(:,1)';ICs_com_before(:,2)']','LineStyle','none','Nbins',nbins,'CdataMode','auto'); colorbar; view(2);
% subplot 122;
% hist3([ICs_com_after(:,1)';ICs_com_after(:,2)']','LineStyle','none','Nbins',nbins,'CdataMode','auto'); colorbar; view(2);

% PR
tmp = cov(ICs_all_before');
PR_before = (sum(eig(tmp))^2) / sum(eig(tmp).^2)
tmp = cov(ICs_all_after');
PR_after = (sum(eig(tmp))^2) / sum(eig(tmp).^2)

%% 
n_similar_ICs_before = length(unique(cat(2,similar_ICs_before_idx{:})));
n_similar_ICs_after = length(unique(cat(2,similar_ICs_after_idx{:})));
n_similar_events_before = length(unique(cat(2,similar_events_before_idx{:})));
n_similar_events_after = length(unique(cat(2,similar_events_after_idx{:})));
n_similar_both_before = length(intersect(unique(cat(2,similar_ICs_before_idx{:})),unique(cat(2,similar_events_before_idx{:}))));
n_similar_both_after = length(intersect(unique(cat(2,similar_ICs_after_idx{:})),unique(cat(2,similar_events_after_idx{:}))));

disp('  *  *  *   Summary   *  *  *  ');
disp(['spontaneous events before: ',num2str(n_events_before)]);
disp(['Similar ICs to chosen - before: ', num2str(n_similar_ICs_before)]);
disp(['Similar ICs to chosen - before [%]: ',num2str(100*n_similar_ICs_before/n_events_before)]);
disp(['Similar ICs to chosen - before [corr]: ',num2str(median(corrs_before_similar_))]);
disp(['Similar events to chosen - before: ', num2str(n_similar_events_before)]);
disp(['Similar events to chosen - before [%]: ',num2str(100*n_similar_events_before/n_events_before)]);
disp(['Similar events & ICs - before: ', num2str(n_similar_both_before)]);
disp(['Similar events & ICs - before [%]: ',num2str(100*n_similar_both_before/n_similar_ICs_before)]);
% disp(['Similar events & ICs - before [corr]: ',num2str()]);
disp(['Similar ICs, burst different - before: ', num2str(n_similar_ICs_before-n_similar_both_before)]);
disp(['Similar ICs, burst different - before [%]: ',num2str(100*(n_similar_ICs_before-n_similar_both_before)/n_similar_ICs_before)]);
disp(['Similar ICs, burst different - before [corr]: ',num2str(median(corrs_before_NOTsimilar_))]);
disp('* * *');
disp(['spontaneous events after: ',num2str(n_events_after)]);
disp(['Similar ICs to chosen - after: ',num2str(n_similar_ICs_after)]);
disp(['Similar ICs to chosen - after [%]: ',num2str(100*n_similar_ICs_after/n_events_after)]);
disp(['Similar ICs to chosen - after [corr]: ',num2str(median(corrs_after_similar_))]);
disp(['Similar events to chosen - after: ', num2str(n_similar_events_after)]);
disp(['Similar events to chosen - after [%]: ',num2str(100*n_similar_events_after/n_events_after)]);
disp(['Similar events & ICs - after: ', num2str(n_similar_both_after)]);
disp(['Similar events & ICs - after [%]: ',num2str(100*n_similar_both_after/n_similar_ICs_after)]);
% disp(['Similar events & ICs - after [corr]: ',num2str()]);
disp(['Similar ICs, burst different - after: ', num2str(n_similar_ICs_after-n_similar_both_after)]);
disp(['Similar ICs, burst different - after [%]: ',num2str(100*(n_similar_ICs_after-n_similar_both_after)/n_similar_ICs_after)]);
disp(['Similar ICs, burst different - after [corr]: ',num2str(median(corrs_after_NOTsimilar_))]);
disp('  *  *  *  *  *  *  *  *  *  *  *  *  *');


%% save
save([dir_,'\ICs.mat'],'-V7.3')