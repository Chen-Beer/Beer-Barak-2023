clc; clear; %close all;

exps_list__ = {'26550_1.11','26549_8.11','38428_17.11','38427_1.11','26532_2.3','26550_3.5',...
             '38426_2.11','26549_11.11','26550_15.11',...
             'left out\38428_20.2','left out\broken_8.11'};

existing_evoked_list_all = [[1  1  1];[0  1  0];[1  0  0];[1  0  1];[0  1  0];[0  0  1];...
                            [0  1  0];[1  0  0];[0  1  1];...
                            [1  1  1];[1  0  1]];
T = 20;
a = 10; b = 12;
n_patterns = 3;

for expi = 1:length(exps_list__)

    expi 

    % load 
load(['..\data\',exps_list__{expi},'\metrics.mat'],'dir_')
load(['..\..\second protocol\',dir_(4:end),'\ICs.mat'],'corrs_events_before_before','corrs_events_before_after','corrs_events_after_after','corrs_ICs_before_before','corrs_ICs_before_after');
load(['..\data\',exps_list__{expi},'\part2.mat'],'similarity_corr_th')
load(['..\data\',exps_list__{expi},'\evoked_existence_clusters.mat']);   
load(['..\data\',exps_list__{expi},'\spontaneous_clusters_separate.mat'])

main_similar_clusters(existing_evoked_list_all(expi,:) ==0) = [];
main_similar_clusters(main_similar_clusters==0) = [];
main_similar_clusters = unique(main_similar_clusters);
main_not_similar_clusters(main_not_similar_clusters==0) = [];
main_not_similar_clusters = unique(main_not_similar_clusters);

Nc = min(length(main_similar_clusters),length(main_not_similar_clusters));
Nc1 = Nc;%length(main_similar_clusters);
Nc2 = Nc;%length(main_not_similar_clusters);
nn = 25;

% indices
similar_idx = cell(Nc1,1);
for p = 1:Nc1
    similar_idx{p} = clusters_valid_before{main_similar_clusters(p)};%(1:nn);
end
non_similar_idx = cell(Nc2,1);
for p = 1:Nc2
    non_similar_idx{p} = clusters_valid_before{main_not_similar_clusters(p)};%(1:nn);
end

%% ICs

similarity_precentage = 0;
similar_ICs_idx_before = cell(Nc1,1);
similar_ICs_idx_after = cell(Nc1,1);
not_similar_ICs_idx_before = cell(Nc2,1);
not_similar_ICs_idx_after = cell(Nc2,1);

for p = 1:Nc1
    tmp = corrs_ICs_before_before(similar_idx{p},:);
    tmp(tmp >= similarity_corr_th) = 1; 
    tmp(tmp < similarity_corr_th) = 0;
    similar_ICs_idx_before{p} = find(sum(tmp) > similarity_precentage*nn);
    tmp = corrs_ICs_before_after(similar_idx{p},:);
    tmp(tmp >= similarity_corr_th) = 1; 
    tmp(tmp < similarity_corr_th) = 0;
    similar_ICs_idx_after{p} = find(sum(tmp) > similarity_precentage*nn);  
end

for p = 1:Nc2
    tmp = corrs_ICs_before_before(non_similar_idx{p},:);
    tmp(tmp >= similarity_corr_th) = 1; 
    tmp(tmp < similarity_corr_th) = 0;
    not_similar_ICs_idx_before{p} = find(sum(tmp) > similarity_precentage*nn);
    tmp = corrs_ICs_before_after(non_similar_idx{p},:);
    tmp(tmp >= similarity_corr_th) = 1; 
    tmp(tmp < similarity_corr_th) = 0;
    not_similar_ICs_idx_after{p} = find(sum(tmp) > similarity_precentage*nn);    
end

%% CDF

theta_vec = 0.85:0.05:1.1;
corr_vec = 0:.001:1;

delta_similar = zeros(Nc1,length(theta_vec));
delta_not = zeros(Nc2,length(theta_vec));
% similar
for p = 1:Nc1
    h_b = histogram(corrs_events_before_before(similar_ICs_idx_before{p},:),corr_vec,'Normalization','cdf','visible','off');
    cdf_b = h_b.Values;
    h_a = histogram(corrs_events_after_after(similar_ICs_idx_after{p},:),corr_vec,'Normalization','cdf','visible','off');
    cdf_a = h_a.Values;
    for t = 1:length(theta_vec)
        [~,m] = min(abs(corr_vec-theta_vec(t)*similarity_corr_th));
        delta_similar(p,t) = cdf_b(m)-cdf_a(m);
    end
end

% not
for p = 1:Nc2
    h_b = histogram(corrs_events_before_before(not_similar_ICs_idx_before{p},:),corr_vec,'Normalization','cdf','visible','off');
    cdf_b = h_b.Values;
    h_a = histogram(corrs_events_after_after(not_similar_ICs_idx_after{p},:),corr_vec,'Normalization','cdf','visible','off');
    cdf_a = h_a.Values;
    for t = 1:length(theta_vec)
        [~,m] = min(abs(corr_vec-theta_vec(t)*similarity_corr_th));
        delta_not(p,t) = cdf_b(m)-cdf_a(m);
    end    
end

%% 
save(['..\data\',exps_list__{expi},'\IC_flattening_effect_delta_V4.mat'],'delta_similar','delta_not','not_similar_ICs_idx_before','not_similar_ICs_idx_after','similar_ICs_idx_after','similar_ICs_idx_before','n_events_before','n_events_after');


end

%% *************************************** PLOT ALL **********************************************

existing_evoked_list_all = [[1  1  1];[0  1  0];[1  0  0];[1  0  1];[0  1  0];[0  0  1];...
                            [0  1  0];[1  0  0];[0  1  1];...
                            [1  1  1];[1  0  1]];

z_similar = cell(size(existing_evoked_list_all,1),length(theta_vec));
z_not = cell(size(existing_evoked_list_all,1),length(theta_vec));

for ex = 1:length(exps_list__)
    ex
    load(['..\data\',exps_list__{ex},'\IC_flattening_effect_delta_V4.mat'])
    all_ = cat(1,delta_not,delta_similar); all_(all_==0) = [];    
    % z-core
    for t = 1:length(theta_vec)
        all_vals = all_(:,t);
        mu = nanmean(all_vals);
        sig = nanstd(all_vals);
        z_similar{ex,t} = (delta_similar(:,t)-mu)./sig;
        z_not{ex,t} = (delta_not(:,t)-mu)./sig;
%         z_similar{ex,t} = delta_similar(:,t);
%         z_not{ex,t} = delta_not(:,t);        
    end
end

% significance
n_shuf = 1e5;
p_val = zeros(length(theta_vec),1);
for t = 1:length(theta_vec)
    similar_ = cat(1,z_similar{:,t});
    not_ = cat(1,z_not{:,t});
    labels = [ones(length(similar_),1); zeros(length(not_),1)];
    tmp_ = [similar_; not_]; 
    delta_true = nanmedian(tmp_(labels == 1)) - nanmedian(tmp_(labels == 0));
    delta_shuf_vec = zeros(n_shuf,1);
    parfor i = 1:n_shuf
        labels_ = labels(randperm(length(labels)));
        m = nanmedian(tmp_(labels_ == 1)) - nanmedian(tmp_(labels_ == 0));
        delta_shuf_vec(i) = m;
    end
    p_val(t) = length(find(delta_shuf_vec < delta_true))/n_shuf;
end

not_ = cell2mat(cat(2,z_not));
similar_ = cell2mat(cat(2,z_similar));

% violin
figure; hold on; 
plot([0 length(theta_vec)+1],[0 0],'--','color',[.5 .5 .5])
data1 = not_;
data2 = similar_;
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

xlabel('\alpha')
ylabel('Z(\Delta CDF)')
ylim([-1.5 2.5])