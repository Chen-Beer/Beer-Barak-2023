clc; clear; 

exps = {'broken_8.11','26532_17.11','38428_20.2','control_38427_7.2'};

exp_id = 4;

load(['..\..\..\_The complete story\data\left out\',exps{exp_id},'\ICs.mat'])
load(['..\..\..\_The complete story\data\left out\',exps{exp_id},'\2 probe\responsiveness_analysis.mat'],'electrodes_names','stimulus_electrods_idx');

T = 20;

%%
all_electrodes = 1:N;

electrodes_chosen = cell(n_patterns,1);
for p = 1:n_patterns
    electrodes_chosen{p} = all_electrodes(~ismember(all_electrodes,stimulus_electrods_idx{chosen_patterns(p)}));
end

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

probe_com = cell(n_patterns,1);
for p = 1:n_patterns
    responses = squeeze(probe_responses_smoothed(chosen_patterns(p),:,:,1:70));  
    probe_com{p} = zeros(size(responses,1),2,70);
    for i = 1:size(responses,1)
        res = squeeze(responses(i,:,:));
        probe_com{p}(i,:,:) = create_COM(res,electrodes_chosen{p},electrode_names_ordered,electrodes_names);
    end
end

%% 5 metrics

% before-before-all
disp('Metrics: before-before...');
[IC_com_distance_before_before,IC_com_corr_before_before,com_distance_before_before,com_corr_before_before,spatial_profile_before_before,IC_spatial_profile_before_before,temporal_profile_before_before,corrs_before_before,IC_corr_before_before] = create_metrics_matrices(events_before,events_com_before,events_before,events_com_before,T,all_electrodes,0);
% before-before1
[IC_com_distance_before_before1,IC_com_corr_before_before1,com_distance_before_before1,com_corr_before_before1,spatial_profile_before_before1,IC_spatial_profile_before_before1,temporal_profile_before_before1,corrs_before_before1,IC_corr_before_before1] = create_metrics_matrices(events_before,events_com_before,events_before,events_com_before,T,electrodes_chosen{1},0);
% before-before2
[IC_com_distance_before_before2,IC_com_corr_before_before2,com_distance_before_before2,com_corr_before_before2,spatial_profile_before_before2,IC_spatial_profile_before_before2,temporal_profile_before_before2,corrs_before_before2,IC_corr_before_before2] = create_metrics_matrices(events_before,events_com_before,events_before,events_com_before,T,electrodes_chosen{2},0);
% before-before3
[IC_com_distance_before_before3,IC_com_corr_before_before3,com_distance_before_before3,com_corr_before_before3,spatial_profile_before_before3,IC_spatial_profile_before_before3,temporal_profile_before_before3,corrs_before_before3,IC_corr_before_before3] = create_metrics_matrices(events_before,events_com_before,events_before,events_com_before,T,electrodes_chosen{3},0);

% probeChosen-before
disp('Metrics: stimulated probe -before...');
IC_com_distance_chosen_before = cell(n_patterns,1);
IC_com_corr_chosen_before = cell(n_patterns,1);
com_distance_chosen_before = cell(n_patterns,1);
com_corr_chosen_before = cell(n_patterns,1);
spatial_profile_chosen_before = cell(n_patterns,1);
IC_spatial_profile_chosen_before = cell(n_patterns,1);
temporal_profile_chosen_before = cell(n_patterns,1);
corrs_chosen_before = cell(n_patterns,1);
IC_corr_chosen_before = cell(n_patterns,1);
for p = 1:n_patterns
    tmp = squeeze(probe_responses_smoothed(p,:,:,:));
%     tmp_com = probe_com_allElectrodes{chosen_patterns(p)};
    tmp_com = probe_com{p};
    [IC_com_distance_chosen_before{p},IC_com_corr_chosen_before{p},com_distance_chosen_before{p},com_corr_chosen_before{p},spatial_profile_chosen_before{p},IC_spatial_profile_chosen_before{p},temporal_profile_chosen_before{p},corrs_chosen_before{p},IC_corr_chosen_before{p}] = create_metrics_matrices(events_before,events_com_before,tmp,tmp_com,T,electrodes_chosen{p},1);
end

% before-after-all
disp('Metrics: before-after...');
[IC_com_distance_before_after,IC_com_corr_before_after,com_distance_before_after,com_corr_before_after,spatial_profile_before_after,IC_spatial_profile_before_after,temporal_profile_before_after,corrs_before_after,IC_corr_before_after] = create_metrics_matrices(events_before,events_com_before,events_after,events_com_after,T,all_electrodes,0);

% after-after-all
disp('Metrics: after-after...');
[IC_com_distance_after_after,IC_com_corr_after_after,com_distance_after_after,com_corr_after_after,spatial_profile_after_after,IC_spatial_profile_after_after,temporal_profile_after_after,corrs_after_after,IC_corr_after_after] = create_metrics_matrices(events_after,events_com_after,events_after,events_com_after,T,all_electrodes,0);
% after-after1
[IC_com_distance_after_after1,IC_com_corr_after_after1,com_distance_after_after1,com_corr_after_after1,spatial_profile_after_after1,IC_spatial_profile_after_after1,temporal_profile_after_after1,corrs_after_after1,IC_corr_after_after1] = create_metrics_matrices(events_after,events_com_after,events_after,events_com_after,T,electrodes_chosen{1},0);
% after-after2
[IC_com_distance_after_after2,IC_com_corr_after_after2,com_distance_after_after2,com_corr_after_after2,spatial_profile_after_after2,IC_spatial_profile_after_after2,temporal_profile_after_after2,corrs_after_after2,IC_corr_after_after2] = create_metrics_matrices(events_after,events_com_after,events_after,events_com_after,T,electrodes_chosen{2},0);
% after-after3
[IC_com_distance_after_after3,IC_com_corr_after_after3,com_distance_after_after3,com_corr_after_after3,spatial_profile_after_after3,IC_spatial_profile_after_after3,temporal_profile_after_after3,corrs_after_after3,IC_corr_after_after3] = create_metrics_matrices(events_after,events_com_after,events_after,events_com_after,T,electrodes_chosen{3},0);

% probeChosen-after
disp('Metrics: stimulated probe -after...');
IC_com_distance_chosen_after = cell(n_patterns,1);
IC_com_corr_chosen_after = cell(n_patterns,1);
com_distance_chosen_after = cell(n_patterns,1);
com_corr_chosen_after = cell(n_patterns,1);
spatial_profile_chosen_after = cell(n_patterns,1);
IC_spatial_profile_chosen_after = cell(n_patterns,1);
temporal_profile_chosen_after = cell(n_patterns,1);
corrs_chosen_after = cell(n_patterns,1);
IC_corr_chosen_after = cell(n_patterns,1);
for p = 1:n_patterns
    tmp = squeeze(probe_responses_smoothed(p,:,:,:));
%     tmp_com = probe_com_allElectrodes{chosen_patterns(p)};
    tmp_com = probe_com{p};
    [IC_com_distance_chosen_after{p},IC_com_corr_chosen_after{p},com_distance_chosen_after{p},com_corr_chosen_after{p},spatial_profile_chosen_after{p},IC_spatial_profile_chosen_after{p},temporal_profile_chosen_after{p},corrs_chosen_after{p},IC_corr_chosen_after{p}] = create_metrics_matrices(events_after,events_com_after,tmp,tmp_com,T,electrodes_chosen{p},1);
end

%% save
save(['..\..\..\_The complete story\data\left out\',exps{exp_id},'\metrics.mat'],'-V7.3');
disp('Done.')