%% First: load "comparing_events_slidingWin_V2.mat"

avg_duration_spont = round(mean(cell2mat(cellfun(@width,events_smoothed_before,'UniformOutput',false))));
disp(['Average event duration: ',num2str(avg_duration_spont*dt),' sec']);

%% plot evoked 
if avg_duration_spont > t2_
    avg_duration_spont = t2_;
end
plot_cluster_com_trajectories(find(labels_evoked_patterns == 1),events_CoM_evoked_(:,:,1:avg_duration_spont),a,b,dt,1); title('pattern #1');
plot_cluster_com_trajectories(find(labels_evoked_patterns == 2),events_CoM_evoked_(:,:,1:avg_duration_spont),a,b,dt,1); title('pattern #2');
plot_cluster_com_trajectories(find(labels_evoked_patterns == 3),events_CoM_evoked_(:,:,1:avg_duration_spont),a,b,dt,1); title('pattern #3');
drawnow;

avg_com_evoked = zeros(n_patterns,2,length(1:avg_duration_spont));
avg_res_evoked = zeros(n_patterns,length(chosen_electrodes_evoked),length(1:avg_duration_spont));
for p = 1:n_patterns
    avg_com_evoked(p,:,:) = squeeze(mean(events_CoM_evoked_(labels_evoked_patterns == p,:,1:avg_duration_spont)));
    avg_res_evoked(p,:,:) = squeeze(mean(events_smoothed_evoked_(labels_evoked_patterns == p,:,1:avg_duration_spont)));
end

%% find typical spont. trajectories - before stimulation

max_leaf = 70;
method = 'weighted';
min_cluster_size = 50; % [events]
plot_label = 0;

metric = 3; % (1) max cross-corr (2) com (3) diagonal ordering

% spont. before
if metric == 1
    D = max_cc_before; D = 1-D;
    title_ = 'spontaneous before - max cc';
end
if metric == 2
    D = com_mean_dist_before;
    title_ = 'spontaneous before - com';
end
if metric == 3
    D = diagonal_ordering_before;
    D = D + 1; D = 2 - D;
    title_ = 'spontaneous before - diagonal ordering';
end
D(D < 0) = 0;
D(isnan(D)) = max(max(D));   
[avg_com_before,avg_res_before,cluster_sizes_before] = plot_clusters_V2(title_,D,events_CoM_before,events_smoothed_before,method,max_leaf,a,b,dt,min_cluster_size,plot_label,1,avg_duration_spont);
drawnow;
   
%% rate through time

time_before = units*time_frames_before(:,1);
time_during = units*(time_frames_before(end,1)+time_frames_during(:,1));
time_after = units*(time_frames_during(end,1)+time_frames_after(:,1));

%% quantify similarity - timeline (6 metrics)

% [D_6metrics,D_evoked_6metrics] = quantify_similarity_individual_events(avg_com_evoked,avg_res_evoked,avg_com_before,avg_res_before,events_smoothed_evoked,events_CoM_evoked,events_CoM_before,events_CoM_during,events_CoM_after,events_smoothed_before,events_smoothed_during,events_smoothed_after,chosen_electrodes_evoked,avg_duration_spont);
D_6metrics = quantify_similarity_individual_events(avg_com_evoked,avg_res_evoked,avg_com_before,avg_res_before,events_smoothed_evoked,events_CoM_evoked,events_CoM_before,events_CoM_during,events_CoM_after,events_smoothed_before,events_smoothed_during,events_smoothed_after,chosen_electrodes_evoked,avg_duration_spont);
save([path_save,'D_6metrics.mat'],'D_6metrics')

%% plot

plot_avg_spont = 1;
smooth_plot = 1;

n_patterns = size(events_smoothed_evoked,1);
n_before = length(events_smoothed_before);
n_during = length(events_smoothed_during);
n_after = length(events_smoothed_after);

n_clusters_before = length(avg_com_before);
    

titles = {'max cross-correlation','center of mass - mean euclidean distance','temporal profile - corr','spatial profile - corr','diagonal ordering','EMD'};
legend_str = {};
for i = 1:n_patterns
    legend_str{i} = ['E',num2str(i)];
end
if plot_avg_spont
    legend_str{i+1} = '<spont. before>';
else
    for j = 1:n_clusters_before
        i = i+1;
        legend_str{i} = ['S',num2str(j)];
    end
end

% for m = 1:size(D_6metrics,1)
%     D = squeeze(D_6metrics(m,:,:));
%     figure; hold on;
%     imagesc(D');
%     title(titles{m});
%     y_labels = {'evoked 1','evoked 2','evoked 3'};
%     for i = 1:sum(n_before)
%         y_labels = cat(2,y_labels,['before ',num2str(i)]);
%     end
%     set(gca,'YTick',1:1:size(D,2), 'YTickLabel', y_labels)
%     plot((n_before+.5)*[1 1],[.5 size(D,2)],'-k','LineWidth',3);
%     plot((n_before+n_during+.5)*[1 1],[.5 size(D,2)],'-k','LineWidth',3);
%     plot([.5 size(D_6metrics,2)+.5],(n_patterns+.5)*[1 1],'-w','LineWidth',3);
%     ylim([.5 size(D,2)]); xlim([.5 size(D,1)+.5]);
%     colorbar;
% end

for m = 1:size(D_6metrics,1)
    D = squeeze(D_6metrics(m,:,:))';
    figure; hold on;
    if smooth_plot
        plot(smoothdata(D(1:n_patterns,:),2,'movmedian',10)','-','LineWidth',2);
    else       
        plot(D(1:n_patterns,:)','-','LineWidth',2);
    end
    if plot_avg_spont
        if smooth_plot
            plot(smoothdata(mean(D((n_patterns+1):end,:)),2,'movmedian',10)','color',[.5 .5 .5],'LineWidth',2);
        else
            plot(mean(D((n_patterns+1):end,:)),'color',[.5 .5 .5],'LineWidth',2);        
        end
    else
        if smooth_plot
            plot(smoothdata(D((n_patterns+1):end,:),2,'movmedian',10)','-','LineWidth',2);
        else
            plot(D((n_patterns+1):end,:)','-','LineWidth',2);        
        end
    end
    plot((.5+n_before)*[1 1],[0 max(max(max(D)))],'-k','LineWidth',3);
    plot((.5+n_before+n_during)*[1 1],[0 max(max(max(D)))],'-k','LineWidth',3);        
    ylabel('dist'); xlabel('Time'); xlim([1 size(D,2)]);
    legend(legend_str,'Location','best');
    sgtitle(titles{m});
end  
drawnow; 

%% quantify similarity - timeline (6 metrics) - CONTROL - 20 patterns (probing)

probing_analysis = [path_save,'2 probe/responsiveness_analysis.mat'];
probing_struct = load(probing_analysis,'center_of_mass','responses_smoothed','stimulus_electrods');

n_probs = length(probing_struct.center_of_mass);

probing_stim_electrodes_idx = cell(n_probs,1);
for i = 1:n_probs
    tmp = [];
    for j = 1:length(probing_struct.stimulus_electrods{i})
        tmp = [tmp find(strcmp(electrodes_names, probing_struct.stimulus_electrods{i}{j}))];
    end
    probing_stim_electrodes_idx{i} = tmp;
end

avg_CoM_probs = zeros(n_probs,2,length(t1_:t2_));
avg_smoothed_probs = zeros(n_probs,N-length(probing_struct.stimulus_electrods{i}),length(t1_:t2_));
for i = 1:n_probs
    events_CoM_probs(i,:,:) = squeeze(mean(probing_struct.center_of_mass{i}));
    events_smoothed_probs(i,:,:) = squeeze(mean(probing_struct.responses_smoothed(i,:,:)));

end