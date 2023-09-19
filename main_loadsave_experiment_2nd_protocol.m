clc; 
clear; 
close all;

control_label = 0;

% dir_ = '../data/control_38427_24.4/4 exp/'; round_cycles = 3; control_label = 1;
dir_ = '../data/39740_2.11/4 exp/'; round_cycles = 3;

% Protocol:
%       3 patterns - distinct areas (not random)
%       6 electrodes each
%       isi 1sec (constant)
%       after each cycle (1-2-3): 1min of spontaneous activity

%% load data
disp('Loading data...');

load_h5_2nd_protocol;

%% Stimulating electrodes 

% 3 patterns
stimulus_electrods = {{'A01','B01','C01','A02','B02','C02'},...
                      {'K03','L03','M03','K04','L04','M04'},...
                      {'D09','E09','F09','D10','E10','F10'}};

if control_label
    stimulus_electrods = {{},...
                          {},...
                          {}};
end


disp('Stimulating electrodes...');
define_stimulating_electrodes;

%% extract responses & remove stimulation artifact "NS"

pre = 10e-3;
pulse = 400; %[micro-sec]
cycles = 1;

% delete spikes shortly after stimulation offset
exclude = 5000; %[micro-sec]

disp('Extracting raw responses & removing stimulatin artifacts...');
extract_raw_responses_spont;  
% sequence = [1,2,3,1,2,3,1,2,3,1,2,3,1,2];
% extract_raw_responses_spont_cycles;

%% create smoothed & binned time series

% temporal resolusion
dt = 5e-3; 
% kernel
sigma = 2e-2;
kernel_frame = 6; %[*sigma]
units = 1e-6; %[sec]

disp('Creating smoothed & binned time series...');
continuous_smoothed_responses_spont;
% continuous_smoothed_responses_spont_cycles;

%% save .mat
disp('Saving mat...');
save([dir_,'data_excluding5msec.mat'], '-v7.3')

disp('Done.');
        
