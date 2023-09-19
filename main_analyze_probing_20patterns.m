%% check responsiveness for all patterns 

clc; clear; close all;

% load('..\stimulation\new batch\38427_4.10\exp_24.10_20days\2 probe\data_excluding5msec.mat');
% load('..//data/18days_4.11/2 probe/data_excluding5msec.mat');

% response subset
t1_ = 1;
t2_ = 70; %[dt]

%% relevant electrodes for each pattern: (active & not stimulating)

relevant_electrodes_per_pattern;

%% detect blank responses

min_overall_activity = 1.5;

valid_responses1;

%% (1) metric: cross-correlation

metric_max_cross_correlation;


%% (2) metric: center of mass trajectory in 2D

metric_center_of_mass2;

%% responsiveness 
% average across patterns & repetitions
figure; 
plot(squeeze(mean(mean(squeeze(sum(responses_smoothed(:,:,non_stimulating_electrodes,t1_:t2_),3))))))
title('<\Sigma_e>'); xlabel('time [dt]');

%% hierarchical clustering 

% method = 'average';
% method = 'median';
method = 'weighted';

% clustering_tmp;

clustering_tmp2_;

%% chosen patterns

chosen_patterns = [5 8 16];

% chosen electrodes 
chosen_stimulating_electrodes = [];
for i = 1:length(chosen_patterns)
    chosen_stimulating_electrodes = [chosen_stimulating_electrodes stimulus_electrods_idx{chosen_patterns(i)}];
end

chosen_electrodes = 1:1:N;
for e = chosen_stimulating_electrodes
    chosen_electrodes(chosen_electrodes==e) = [];
end

chosen_patterns_view;

%% save .mat
disp('Saving mat...');
save([dir_,'responsiveness_analysis.mat'], '-v7.3')

disp('Done.');