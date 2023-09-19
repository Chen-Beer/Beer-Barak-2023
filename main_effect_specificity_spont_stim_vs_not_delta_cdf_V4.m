clc; clear; close all;
addpath(genpath('distributionPlot'));

exps_list__ = {'26550_1.11','26549_8.11','38428_17.11','38427_1.11','26532_2.3','26550_3.5',...
              '38426_2.11','26549_11.11','26550_15.11',...
              'control_26550_24.1','control_39740_24.4','control_38427_24.4','control_26536_7.2',...
              'left out\38428_20.2','left out\broken_8.11','left out\control_38427_7.2'};

existing_evoked_list_all = [[1  1  1];[0  1  0];[1  0  0];[1  0  1];[0  1  0];[0  0  1];...
                            [0  1  0];[1  0  0];[0  1  1];...
                            [1  1  1];[0  0  1];[1  1  1];[1  1  1];...
                            [1  1  1];[1  0  1];[1  0  1]];


theta_vec = 0.85:0.05:1.1;
corr_vec = 0:.001:1;
option = 1;

delta_cdf_similar = cell(length(exps_list__),length(theta_vec));
delta_cdf_not = cell(length(exps_list__),length(theta_vec));
z_similar = cell(length(exps_list__),length(theta_vec));
z_not = cell(length(exps_list__),length(theta_vec));

for expi = 1:length(exps_list__)

disp(expi)
disp(exps_list__{expi})

    % load 
load(['..\data\',exps_list__{expi},'\metrics.mat'],'dir_')
load(['..\..\second protocol\',dir_(4:end),'\ICs.mat'],'corrs_events_before_before','corrs_events_before_after');
load(['..\data\',exps_list__{expi},'\part2.mat'],'similarity_corr_th')
load(['..\data\',exps_list__{expi},'\evoked_existence_clusters.mat']);   
load(['..\data\',exps_list__{expi},'\spontaneous_clusters_separate.mat'])

main_similar_clusters(existing_evoked_list_all(expi,:) ==0) = [];
main_similar_clusters(main_similar_clusters==0) = [];
main_similar_clusters = unique(main_similar_clusters);
main_not_similar_clusters(main_not_similar_clusters==0) = [];
main_not_similar_clusters = unique(main_not_similar_clusters);

Nc = min(length(main_similar_clusters),length(main_not_similar_clusters));
Nc1 = Nc; %length(main_similar_clusters);
Nc2 = Nc; %length(main_not_similar_clusters);
nn = 25;

% indices
similar_idx = cell(Nc1,1);
for p = 1:Nc1
    similar_idx{p} = clusters_valid_before{main_similar_clusters(p)}(1:nn);
end
non_similar_idx = cell(Nc2,1);
for p = 1:Nc2
    non_similar_idx{p} = clusters_valid_before{main_not_similar_clusters(p)}(1:nn);
end

% distance - similar
for p = 1:Nc1
    tmp_idx = similar_idx{p};
    b_ = corrs_events_before_before(tmp_idx,:); b_ = b_(:);
    a_ = corrs_events_before_after(tmp_idx,:); a_ = a_(:);
    cdf_a = zeros(length(corr_vec),1);
    cdf_b = zeros(length(corr_vec),1);
    for c = 1:length(corr_vec)
        cdf_a(c) = length(find(a_ <= corr_vec(c)))/length(a_);
        cdf_b(c) = length(find(b_ <= corr_vec(c)))/length(b_);
    end
    for t = 1:length(theta_vec)
        if option == 1
            [~,m] = min(abs(corr_vec-theta_vec(t)*similarity_corr_th));
            tmp = cdf_b(m)-cdf_a(m);
        else
            [~,i_before] = min(abs(cdf_b-theta_vec(t)));
            [~,i_after] = min(abs(cdf_a-theta_vec(t)));
            tmp = corr_vec(i_before) - corr_vec(i_after);
        end
        delta_cdf_similar{expi,t} = cat(2,delta_cdf_similar{expi,t},tmp);
    end
    
end

% distance - not
delta_cdf_not{expi} = [];
for p = 1:Nc2
    tmp_idx = non_similar_idx{p};
    b_ = corrs_events_before_before(tmp_idx,:); b_ = b_(:);
    a_ = corrs_events_before_after(tmp_idx,:); a_ = a_(:);
    cdf_a = zeros(length(corr_vec),1);
    cdf_b = zeros(length(corr_vec),1);
    for c = 1:length(corr_vec)
        cdf_a(c) = length(find(a_ <= corr_vec(c)))/length(a_);
        cdf_b(c) = length(find(b_ <= corr_vec(c)))/length(b_);
    end
    for t = 1:length(theta_vec)
        if option == 1
            [~,m] = min(abs(corr_vec-theta_vec(t)*similarity_corr_th));
            tmp = cdf_b(m)-cdf_a(m);
        else
            [~,i_before] = min(abs(cdf_b-theta_vec(t)));
            [~,i_after] = min(abs(cdf_a-theta_vec(t)));
            tmp = corr_vec(i_before) - corr_vec(i_after);
        end        
        delta_cdf_not{expi,t} = cat(2,delta_cdf_not{expi,t},tmp);
    end
end

% z-core
for t = 1:length(theta_vec)
    all_vals = [delta_cdf_not{expi,t},delta_cdf_similar{expi,t}];
    mu = mean(all_vals);
    sig = std(all_vals);
    z_similar{expi,t} = (delta_cdf_similar{expi,t}-mu)./sig;
    z_not{expi,t} = (delta_cdf_not{expi,t}-mu)./sig;
end


end

%% ctrl vs stimulation

exps_subset = [1:9,14:15];
ctrl_subset = [10:13,16];

z_not_stim = z_not(exps_subset,:);
z_not_ctrl = z_not(ctrl_subset,:);
z_similar_stim = z_similar(exps_subset,:);
z_similar_ctrl = z_similar(ctrl_subset,:);
delta_cdf_not_stim = delta_cdf_not(exps_subset,:);
delta_cdf_not_ctrl = delta_cdf_not(ctrl_subset,:);
delta_cdf_similar_stim = delta_cdf_similar(exps_subset,:);
delta_cdf_similar_ctrl = delta_cdf_similar(ctrl_subset,:);

%% *************************************** PLOT ALL  - specificity **********************************************
n_ = length(cat(2,z_not_stim{:}))/length(theta_vec);
not_ = zeros(n_,length(theta_vec));
n_ = length(cat(2,z_similar_stim{:}))/length(theta_vec);
similar_ = zeros(n_,length(theta_vec));
for t = 1:length(theta_vec)
    not_(:,t) = cat(2,z_not_stim{:,t});
    similar_(:,t) = cat(2,z_similar_stim{:,t});
end

% violin
figure; hold on; 
plot([0 length(theta_vec)+1],[0 0],'--','color',[.5 .5 .5])
data1 = not_;
data2 = similar_;
distributionPlot(gca,data1,'widthDiv',[2 1],'histOri','left','addSpread',true,'color','g','showMM',false,'histOpt',1.1)
distributionPlot(data2,'widthDiv',[2 2],'histOri','right','addSpread',true,'color','r','showMM',false,'histOpt',1.1)
set(gca,'XTickLabel',theta_vec);
ylabel('Z(\Delta CDF)'); xlabel('\alpha');

figure; hold on; 
vs = violinplot(data2, theta_vec,...
    'HalfViolin','right','ViolinColor',[.64 .08 .18],...%red
    'QuartileStyle','shadow',... % boxplot, none
    'DataStyle', 'scatter',... % scatter, none
    'ShowNotches', false,...
    'ShowMean', false,...
    'ShowMedian', true);
vs = violinplot(data1, theta_vec,...
    'HalfViolin','left','ViolinColor',[.47 .67 .19],...%green
    'QuartileStyle','shadow',... % boxplot, none
    'DataStyle', 'scatter',... % scatter, none
    'ShowNotches', false,...
    'ShowMean', false,...
    'ShowMedian', true);
title('Specificity: similar vs. not - spontaneous clusters')
ylabel('Z(\Delta CDF)'); xlabel('\alpha');
plot([0 length(theta_vec)+1],[0 0],'--','color',[.5 .5 .5])

% significance - experiment vs. control
n_shuf = 1e5;
pval_vec = zeros(size(theta_vec));
for t = 1:length(theta_vec)
%     disp(t)
    similar__ = similar_(:,t);
    not__ = not_(:,t);
    labels = [ones(length(similar__),1); zeros(length(not__),1)];
    WD_ = [similar__; not__]; 
    delta_true = median(WD_(labels == 1)) - median(WD_(labels == 0));
    delta_shuf_vec = zeros(n_shuf,1);
    parfor i = 1:n_shuf
        labels_ = labels(randperm(length(labels)));
        m = median(WD_(labels_ == 1)) - median(WD_(labels_ == 0));
        delta_shuf_vec(i) = m;
    end
    p_val = length(find(delta_shuf_vec < delta_true))/n_shuf;
    pval_vec(t) = p_val;
end
pval_vec

%% *************************************** PLOT ALL  - control vs. stimulation **********************************************
ctrl_tmp = cat(1,delta_cdf_similar_ctrl,delta_cdf_not_ctrl);
stim_tmp = cat(1,delta_cdf_similar_stim,delta_cdf_not_stim);

n_ = length(cat(2,ctrl_tmp{:}))/length(theta_vec);
not_ = zeros(n_,length(theta_vec));
n_ = length(cat(2,stim_tmp{:}))/length(theta_vec);
similar_ = zeros(n_,length(theta_vec));
for t = 1:length(theta_vec)
    not_(:,t) = cat(2,ctrl_tmp{:,t});
    similar_(:,t) = cat(2,stim_tmp{:,t});
end

% violin
figure; hold on; 
plot([0 length(theta_vec)+1],[0 0],'--','color',[.5 .5 .5])
data1 = not_;
data2 = similar_;
distributionPlot(gca,data1,'widthDiv',[2 1],'histOri','left','addSpread',true,'color','g','showMM',false,'histOpt',1.1)
distributionPlot(data2,'widthDiv',[2 2],'histOri','right','addSpread',true,'color','r','showMM',false,'histOpt',1.1)
set(gca,'XTickLabel',theta_vec);
ylabel('\Delta CDF'); xlabel('\alpha');
plot([0 length(theta_vec)+1],[0 0],'--','color',[.5 .5 .5])

figure; hold on; 
vs = violinplot(data2, theta_vec,...
    'HalfViolin','right','ViolinColor',[.64 .08 .18],...%red
    'QuartileStyle','shadow',... % boxplot, none
    'DataStyle', 'scatter',... % scatter, none
    'ShowNotches', false,...
    'ShowMean', false,...
    'ShowMedian', true);
vs = violinplot(data1, theta_vec,...
    'HalfViolin','left','ViolinColor',[.47 .67 .19],...%green
    'QuartileStyle','shadow',... % boxplot, none
    'DataStyle', 'scatter',... % scatter, none
    'ShowNotches', false,...
    'ShowMean', false,...
    'ShowMedian', true);
title('Main effect: control vs. stimulation - spontaneous clusters')
ylabel('\Delta CDF'); xlabel('\alpha');
plot([0 length(theta_vec)+1],[0 0],'--','color',[.5 .5 .5])

% significance - experiment vs. control
n_shuf = 1e5;
pval_vec = zeros(size(theta_vec));
for t = 1:length(theta_vec)
    similar__ = similar_(:,t);
    not__ = not_(:,t);
    labels = [ones(length(similar__),1); zeros(length(not__),1)];
    WD_ = [similar__; not__]; 
    delta_true = median(WD_(labels == 1)) - median(WD_(labels == 0));
    delta_shuf_vec = zeros(n_shuf,1);
    parfor i = 1:n_shuf
        labels_ = labels(randperm(length(labels)));
        m = median(WD_(labels_ == 1)) - median(WD_(labels_ == 0));
        delta_shuf_vec(i) = m;
    end
    p_val = length(find(delta_shuf_vec < delta_true))/n_shuf;
    pval_vec(t) = p_val;
end
pval_vec
