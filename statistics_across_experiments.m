clc; clear; 

exps_list__ = {'26550_1.11','26549_8.11','38428_17.11','38427_1.11','26532_2.3','26550_3.5',...
              '38426_2.11','26549_11.11','26550_15.11',...
              'control_26550_24.1','control_39740_24.4','control_38427_24.4','control_26536_7.2',...
              'left out\38428_20.2','left out\broken_8.11','left out\control_38427_7.2'};

for expi = 1:length(exps_list__)

load(['..\data\',exps_list__{expi},'\spontaneous_clusters_separate.mat']);
load(['..\data\',exps_list__{expi},'\metrics.mat'],'dir_','dt');
load(['..\..\second protocol\',dir_(4:end),'\spont_events_individuals.mat'],'time_frames_before','time_frames_after');

T_before = time_frames_before(end,end)*dt/60/60;
T_after = time_frames_after(end,end)*dt/60/60;

disp('************')
disp(exps_list__{expi})
disp(['Before: # of events per hour:  ',num2str(n_events_before/T_before)]);
disp(['Before: # of clusters:  ',num2str(n_valid_clusters_before)]);
disp(['Before: % explained:  ',num2str(sum(cellfun(@length,clusters_valid_before))/n_events_before)]);
disp(['After: # of events per hour:  ',num2str(n_events_after/T_after)]);
disp(['After: # of clusters:  ',num2str(n_valid_clusters_after)]);
disp(['After: % explained:  ',num2str(sum(cellfun(@length,clusters_valid_after))/n_events_after)]);
disp('************')

end