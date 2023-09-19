%% load

clc; clear; close all;

% load('../data/23days_27.10/3 exp/data_excluding5msec.mat');
% load('../data/39740_1.11/17days_18.11/4 exp/data_excluding5msec.mat');
% load('../data/38427_1.11/20days_21.11/4 exp/data_excluding5msec.mat');

% define how many repetitions per pattern there are in a single window 
window = 20; 
step = 10;

% response subset
t1_ = 1;
t2_ = 70; %[dt]

remove_baseline = 0;
baseline_max = 0;
baseline_spatiotemp = 0;

%% valid repetitions & electrodes

min_overall_activity = 1;
min_electrode_activity = 0;

valid_repetitions_and_electrodes;
drawnow;

%% remove baseline activity from each electrode
% based on average activity in the first repetitions
avg_reps = 10;
if remove_baseline 
    baseline_removal;
end

%% generate sliding responses (valid avg)
sliding_responses = zeros(n_windows,n_patterns,E,T);
for w = 1:n_windows
    rep_idx = sliding_win_rep_idx(w,:);
    valid_rep_idx = rep_idx(non_zero_responses(p,rep_idx) == 1);
    valid_rep_idx_end = rep_idx(non_zero_responses_end(p,rep_idx) == 1);
    for p = 1:n_patterns
        sliding_responses(w,p,:,:) = squeeze(mean(responses_smoothed(p,valid_rep_idx,chosen_electrodes,t1_:t2_),2));
%         sliding_responses_end(w,p,:,:) = squeeze(mean(responses_smoothed(p,valid_rep_idx_end,chosen_electrodes,end-t2_+1:end),2));
    end
end

%% data summary
disp('* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *');
disp('Data summary');
disp(['Number of patterns: ',num2str(n_patterns)]);
disp(['Repetitions per pattern: ',num2str(n_repetitions)]);
disp(['Average percentage of valid responses: ', num2str(mean(sum(non_zero_responses,2))*100/n_repetitions)]);
disp(['Number of analyzed electrodes: ', num2str(length(chosen_electrodes))]);
disp(['Duration of analyzed response [sec]: ',num2str(T*dt)]);
disp(['Sliding window size: ',num2str(window)]);
disp(['sliding step size: ',num2str(step)]);
disp(['Number of windows: ',num2str(n_windows)]);
disp('* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *');

%% create ordering vectors for each response (raw data)

responses_ordering_index = nan(n_patterns,n_repetitions,E);
responses_ordering_time = nan(n_patterns,n_repetitions,E);
% responses_ordering_index_end = nan(n_patterns,n_repetitions,E);
% responses_ordering_time_end = nan(n_patterns,n_repetitions,E);
extract_ordering_vectors;

%% metrics - sliding window
disp('calculating distance using 6 metrics:');

% (1) max cross correlation
fprintf('(1) max cross correlation...');
sliding_metric_maxcorr = zeros(n_windows,n_patterns,n_patterns);
% sliding_metric_maxcorr_endstart = zeros(n_windows,n_patterns,n_patterns);
sliding_max_crosscorr;
disp(' done.');

% (2) center of mass - spatio-temporal euclidean distance
fprintf('(2) euclidean distance - center of mass...');
sliding_metric_CoMdist = zeros(n_windows,n_patterns,n_patterns);
% sliding_metric_CoMdist_endstart = zeros(n_windows,n_patterns,n_patterns);
sliding_CoM_euclidean_distance;
disp(' done.');

% (3) sum over time
fprintf('(3) sum over time...');
sliding_metric_temporalSum = zeros(n_windows,n_patterns,n_patterns);
% sliding_metric_temporalSum_endstart = zeros(n_windows,n_patterns,n_patterns);
sliding_temporal_summation;
disp(' done.');

% (4) sum over electrodes
fprintf('(4) sum over electrodes...');
sliding_metric_spatialSum = zeros(n_windows,n_patterns,n_patterns);
% sliding_metric_spatialSum_endstart = zeros(n_windows,n_patterns,n_patterns);
sliding_spatial_summation;
disp(' done.');

% (5) ordering - spearman (avg time)
fprintf('(5) ordering - spearman (averaging time)...');
sliding_metric_orderingSpearman = zeros(n_windows,n_patterns,n_patterns);
% sliding_metric_orderingSpearman_endstart = zeros(n_windows,n_patterns,n_patterns);
sliding_ordering_time_spearman;
disp(' done.');

% (6) ordering - pearson (avg index)
fprintf('(6) ordering - pearson (averaging index)...');
sliding_metric_orderingPearson = zeros(n_windows,n_patterns,n_patterns);
% sliding_metric_orderingPearson_endstart = zeros(n_windows,n_patterns,n_patterns);
sliding_ordering_index_pearson;
disp(' done.');

% (7) ordering - edit distance (avg index)
fprintf('(7) ordering - edit distance (avg index)...');
sliding_metric_ordering_editDist = zeros(n_windows,n_patterns,n_patterns);
% sliding_metric_ordering_editDist_endstart = zeros(n_windows,n_patterns,n_patterns);
sliding_ordering_index_editDist;
disp(' done.');

drawnow;

%% slopes - change sin metric values over time

extract_slopes_2nd_protocol(sliding_metric_maxcorr);
sgtitle('max cross corr');
extract_slopes_2nd_protocol(sliding_metric_CoMdist);
sgtitle('center of mass - euclidean distance');
extract_slopes_2nd_protocol(sliding_metric_temporalSum);
sgtitle('temporal summation - corr');
extract_slopes_2nd_protocol(sliding_metric_spatialSum);
sgtitle('spatial summation - corr');
extract_slopes_2nd_protocol(sliding_metric_orderingSpearman);
sgtitle('first spike ordering - spearman');
extract_slopes_2nd_protocol(sliding_metric_orderingPearson);
sgtitle('first spike ordering - pearson');

drawnow;

%% latency
latency = zeros(n_patterns,n_repetitions);
parfor p = 1:n_patterns
    for r = 1:n_repetitions
        latency(p,r) = min(squeeze(responses_ordering_time(p,r,:)));
    end
end

figure; imagesc(latency);
set(gca,'ColorScale','log'); colorbar;
xlabel('repetition'); ylabel('pattern');
title('latency - time to first spike');

drawnow;

%% save .mat
fprintf('Saving mat...');
if remove_baseline
    if baseline_spatiotemp
        save([dir_,'analysis_6metrics_baselineRemoved.mat'], '-v7.3')
    else
        save([dir_,'analysis_6metrics_baselineRemoved_max.mat'], '-v7.3')
    end
else
    save([dir_,'analysis_6metrics.mat'], '-v7.3')
end
disp('Done.');
