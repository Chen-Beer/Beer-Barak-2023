%% load data

clc; clear; %close all;

exp_list_ = {'26550_1.11','26549_8.11','38428_17.11','38427_1.11','26532_2.3','26550_3.5',...
            '38426_2.11','39740_2.11','26549_11.11','26550_15.11'};

% EXCLUDE NON-EXISTING STIMULATED EVOKED
existing_evoked_list = [[1 1 1];[1 1 0];[1 0 0];[1 0 1];[1 0 1];[1 0 1];...
                        [0 1 0];[1 1 0];[1 0 0];[0 1 1]];


%% load all

n_win = 100;
step = 20;
n_patterns = 3;
T = 40;

evoked_variance_all = nan(length(exp_list_),n_patterns,1500/step-5);
for expi = [1:7 9:length(exp_list_)]
    disp([num2str(expi),':', exp_list_{expi}]);
    load(['..\data\',exp_list_{expi},'\metrics.mat'],'dir_','electrodes_chosen');
    load(['..\..\second protocol\data\',dir_(9:end),'\4 exp\data_excluding5msec.mat'],'responses_smoothed');

    evoked_variance = nan(n_patterns,round(size(responses_smoothed,2)/step));
    for p = 1:n_patterns
        if existing_evoked_list(expi,p)
            i = 1; idx = 1;
            while idx+n_win <= size(responses_smoothed,2)
                tmp = squeeze(responses_smoothed(p,idx:idx+n_win,electrodes_chosen{p},1:T));
                mean_tmp = squeeze(mean(tmp));
                norm_ = 0;
                for j = 1:size(tmp,1)
                    norm_ = norm_ + (norm(mean_tmp-squeeze(tmp(j,:,:))));
                end
                evoked_variance(p,i) = norm_/size(tmp,1);
                i = i + 1;
                idx = idx + step;
            end
        end
    end
    evoked_variance = evoked_variance(:,1:i-1);

    for p = 1:n_patterns
%         evoked_variance(p,:) = (evoked_variance(p,:)-min(evoked_variance(p,:)))./(max(evoked_variance(p,:))-min(evoked_variance(p,:)));
        evoked_variance(p,:) = evoked_variance(p,:)./evoked_variance(p,1);
%         evoked_variance(p,:) = evoked_variance(p,:)-evoked_variance(p,1);
    end

    % figure; plot(evoked_variance');
    figure; hold on; 
    plot(linspace(1,10,size(evoked_variance,2)),evoked_variance');
    shadedErrorBar(linspace(1,10,size(evoked_variance,2)),nanmean(evoked_variance),nanvar(evoked_variance));
    ylabel('Variance'); xlabel('Time [h]'); title([exp_list_{expi}]);
    legend('pattern 1','pattern 2','pattern 3');
    drawnow;
    evoked_variance_all(expi,:,:) = evoked_variance(:,1:size(evoked_variance_all,3));
end

%% plot all
tmp = reshape(evoked_variance_all,[],size(evoked_variance_all,3));
figure; 
shadedErrorBar(linspace(0,10,size(tmp,2)),nanmean(tmp),nanvar(tmp));
ylabel('Variance'); xlabel('Time [h]');