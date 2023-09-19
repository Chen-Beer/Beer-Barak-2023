clc; 
clear; 
% close all;

% dir_ = '../data/38426_21.11/2 probe/';
% dir_ = '../data/39740_2.11/2 probe/';
% dir_ = '../data/38428_11.11/2 probe/';
% dir_ = '../data/26549_11.11/2 probe/';
dir_ = '../data/26550_15.11/2 probe/';

%% load data
disp('Loading data...');
load_h5;

%% Stimulating electrodes 

% 20 patterns - 6 electrodes - one channel
stimulus_electrods = {{'A01','B01','C01','A02','B02','C02'},...
                      {'D01','E01','F01','D02','E02','F02'},...
                      {'G01','H01','J01','G02','H02','J02'},...
                      {'K01','L01','M01','K02','L02','M02'},...
                      {'A03','B03','C03','A04','B04','C04'},...
                      {'D03','E03','F03','D04','E04','F04'},...
                      {'G03','H03','J03','G04','H04','J04'},...
                      {'K03','L03','M03','K04','L04','M04'},...
                      {'A05','B05','C05','A06','B06','C06'},...
                      {'D05','E05','F05','D06','E06','F06'},...
                      {'G05','H05','J05','G06','H06','J06'},...
                      {'K05','L05','M05','K06','L06','M06'},...
                      {'A07','B07','C07','A08','B08','C08'},...
                      {'D07','E07','F07','D08','E08','F08'},...
                      {'G07','H07','J07','G08','H08','J08'},...
                      {'K07','L07','M07','K08','L08','M08'},...
                      {'A09','B09','C09','A10','B10','C10'},...
                      {'D09','E09','F09','D10','E10','F10'},...
                      {'G09','H09','J09','G10','H10','J10'},...
                      {'K09','L09','M09','K10','L10','M10'}};

disp('Stimulating electrodes...');
define_stimulating_electrodes;

%% extract responses & remove stimulation artifact "NS"

pre = 100e-3;
pulse = 400; %[micro-sec]
cycles = 1;

% delete spikes shortly after stimulation offset
exclude = 5000; %[micro-sec]

disp('Extracting raw responses & removing stimulatin artifacts...');
extract_raw_responses;

%% create smoothed & binned time series

% temporal resolusion
dt = 5e-3; 
% kernel
sigma = 2e-2;
kernel_frame = 6; %[*sigma]
units = 1e-6; %[sec]

disp('Creating smoothed & binned time series...');
continuous_smoothed_responses;

%% save .mat
disp('Saving mat...');
save([dir_,'data_excluding5msec.mat'], '-v7.3')

disp('Done.');

