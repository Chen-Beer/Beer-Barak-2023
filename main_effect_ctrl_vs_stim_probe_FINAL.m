clc; clear; close all;

addpath(genpath('distributionPlot'));

exp_list1 = {'26550_1.11','26549_8.11','38428_17.11','38427_1.11','26532_2.3','long_26550_3.5','new batch\38426_2.11','new batch\26549_11.11','new batch\26550_15.11','38428_20.2','broken_8.11'};
exp_list2 = {'26550_1.11','26549_8.11','38428_17.11','38427_1.11','26532_2.3','26550_3.5','38426_2.11','26549_11.11','26550_15.11','left out\38428_20.2','left out\broken_8.11'};
ctrl_list1 = {'control_26550_24.1','control_39740_24.4','control_38427_24.4','control_26536_7.2','left out\control_38427_7.2'};
ctrl_list2 = {'control_26550_24.1','control_39740_24.4','control_38427_24.4','control_26536_7.2','control_38427_7.2'};

existing_evoked_stim = [[1  1  1];[0  1  0];[1  0  0];[1  0  1];[0  1  0];[0  0  1];[0  1  0];[1  0  0];[0  1  1];[1  1  1];[1  0  1]];
existing_evoked_ctrl = [[1  1  1];[0  0  1];[0  1  1];[1  1  1];[1  0  1]];

n_probe = 20;
n_patterns = 3;

%% 

theta_vec = .7:.05:1.0;
corr_vec = 0:.005:1;

probe_vec = 1:n_probe;

delta_cdf_exp_chosen = zeros(length(exp_list1),n_patterns,length(theta_vec));
delta_cdf_exp_all = zeros(length(exp_list1),n_probe,length(theta_vec));

for i = 1:length(exp_list1)
    disp(i)
    load(['..\data\',exp_list1{i},'\freshlook_merged_data.mat'],'chosen_patterns','corr_evoked_to_before','corr_evoked1_to_after');
    load(['..\..\_The complete story\data\',exp_list2{i},'\part2.mat'],'similarity_corr_th')
    delta_cdf_theta_vec = zeros(n_probe,length(theta_vec));
    for p = 1:n_probe
        b_ = squeeze(corr_evoked_to_before(p,:,:)); b_ = b_(:); 
        a_ = squeeze(corr_evoked1_to_after(p,:,:)); a_ = a_(:); 
        cdf_a = zeros(length(corr_vec),1);
        cdf_b = zeros(length(corr_vec),1);
        for c = 1:length(corr_vec)
            cdf_a(c) = length(find(a_ <= corr_vec(c)))/length(a_);
            cdf_b(c) = length(find(b_ <= corr_vec(c)))/length(b_);
        end
        for t = 1:length(theta_vec)
            [~,m] = min(abs(corr_vec-similarity_corr_th*theta_vec(t)));
            delta_cdf_theta_vec(p,t) = cdf_b(m) - cdf_a(m);
        end
    end
    delta_cdf_exp_chosen(i,:,:) = delta_cdf_theta_vec(chosen_patterns,:);
    delta_cdf_exp_all(i,:,:) = delta_cdf_theta_vec;
end

delta_cdf_ctrl_chosen = zeros(length(ctrl_list2),n_patterns,length(theta_vec));
delta_cdf_ctrl_all = zeros(length(ctrl_list2),n_probe,length(theta_vec));
for i = 1:length(ctrl_list2)
    load(['..\data\',ctrl_list2{i},'\freshlook_merged_data.mat'],'corr_evoked_to_before','corr_evoked1_to_after');
    load(['..\data\',ctrl_list2{i},'\2 probe\responsiveness_analysis.mat'],'chosen_patterns');
    load(['..\..\_The complete story\data\',ctrl_list1{i},'\part2.mat'],'similarity_corr_th')
    delta_cdf_theta_vec = zeros(n_probe,length(theta_vec));
    for p = 1:n_probe
        b_ = squeeze(corr_evoked_to_before(p,:,:)); b_ = b_(:); 
        a_ = squeeze(corr_evoked1_to_after(p,:,:)); a_ = a_(:); 
        cdf_b = zeros(length(corr_vec),1);
        cdf_a = zeros(length(corr_vec),1);
        for c = 1:length(corr_vec)
            cdf_b(c) = length(find(b_ <= corr_vec(c)))/length(b_);
            cdf_a(c) = length(find(a_ <= corr_vec(c)))/length(a_);
        end
        for t = 1:length(theta_vec)
            [~,m] = min(abs(corr_vec-similarity_corr_th*theta_vec(t)));
            delta_cdf_theta_vec(p,t) = cdf_b(m) - cdf_a(m);
        end
    end
    delta_cdf_ctrl_chosen(i,:,:) = delta_cdf_theta_vec(chosen_patterns,:);
    delta_cdf_ctrl_all(i,:,:) = delta_cdf_theta_vec;
end

%% plot


% violin
figure; hold on; 
plot([0 length(theta_vec)+1],[0 0],'--','color',[.5 .5 .5])
data1 = reshape(delta_cdf_exp_chosen,[],length(theta_vec));
data2 = reshape(delta_cdf_ctrl_chosen,[],length(theta_vec));
% data1 = reshape(delta_cdf_exp_all,[],length(theta_vec));
% data2 = reshape(delta_cdf_ctrl_all,[],length(theta_vec));
distributionPlot(gca,data2,'widthDiv',[2 1],'histOri','left','addSpread',true,'color','g','showMM',false,'histOpt',1.1)
distributionPlot(data1,'widthDiv',[2 2],'histOri','right','addSpread',true,'color','r','showMM',false,'histOpt',1.1)
set(gca,'XTickLabel',theta_vec);
ylabel('\Delta cdf'); xlabel('\theta');

figure; hold on; 
vs = violinplot(data1, theta_vec,...
    'HalfViolin','right','ViolinColor',[.64 .08 .18],...%red
    'QuartileStyle','shadow',... % boxplot, none
    'DataStyle', 'scatter',... % scatter, none
    'ShowNotches', false,...
    'ShowMean', false,...
    'ShowMedian', true);
vs = violinplot(data2, theta_vec,...
    'HalfViolin','left','ViolinColor',[.47 .67 .19],...%green
    'QuartileStyle','shadow',... % boxplot, none
    'DataStyle', 'scatter',... % scatter, none
    'ShowNotches', false,...
    'ShowMean', false,...
    'ShowMedian', true);

% significance - experiment vs. control
n_shuf = 1e5;
pval_vec = zeros(size(theta_vec));
for t = 1:length(theta_vec)
    disp(t)
    exp_ = delta_cdf_exp_chosen(:,:,t); exp_ = exp_(:);
    control_ = delta_cdf_ctrl_chosen(:,:,t); control_ = control_(:);
%     exp_ = delta_cdf_exp_all(:,:,t); exp_ = exp_(:);
%     control_ = delta_cdf_ctrl_all(:,:,t); control_ = control_(:);
    labels = [ones(size(exp_)); zeros(size(control_))];
    WD_ = [exp_; control_]; 
    delta_true = median(WD_(labels == 1)) - median(WD_(labels == 0));
    delta_shuf_vec = zeros(n_shuf,1);
    parfor i = 1:n_shuf
        labels_ = labels(randperm(length(labels)));
        m = median(WD_(labels_ == 1)) - median(WD_(labels_ == 0));
        delta_shuf_vec(i) = m;
    end
    p_val = length(find(delta_shuf_vec < delta_true))/n_shuf;
    pval_vec(t) = p_val;
%     figure; hold on; 
%     histogram(delta_shuf_vec,'normalization','probability');
    
%     title(['p_{val} = ',num2str(p_val), ' | n = ',num2str(n_shuf),' | ',num2str(theta_vec(t)),'\theta']);
%     plot([delta_true delta_true],[0 0.03],'-r','linewidth',2);
%     ylabel('P'); xlabel('\mu_{full xperiment} - \mu_{control}')
%     legend('Shuffled','True');
%     drawnow;
end
pval_vec

%% Plot - excluding non-existing patterns
for i = 1:size(delta_cdf_exp_chosen,1)
    if i == 1
        delta_cdf_exp_existing = squeeze(delta_cdf_exp_chosen(i,existing_evoked_stim(i,:)==1,:));
    else
        if size(squeeze(delta_cdf_exp_chosen(i,existing_evoked_stim(i,:)==1,:)),1) == length(theta_vec)
            delta_cdf_exp_existing = cat(1,delta_cdf_exp_existing,squeeze(delta_cdf_exp_chosen(i,existing_evoked_stim(i,:)==1,:))');
        else
            delta_cdf_exp_existing = cat(1,delta_cdf_exp_existing,squeeze(delta_cdf_exp_chosen(i,existing_evoked_stim(i,:)==1,:)));
        end
    end
end
for i = 1:size(delta_cdf_ctrl_chosen,1)
    if i == 1
        delta_cdf_ctrl_existing = squeeze(delta_cdf_ctrl_chosen(i,existing_evoked_ctrl(i,:)==1,:));
    else
        if size(squeeze(delta_cdf_ctrl_chosen(i,existing_evoked_ctrl(i,:)==1,:)),1) == length(theta_vec)
            delta_cdf_ctrl_existing = cat(1,delta_cdf_ctrl_existing,squeeze(delta_cdf_ctrl_chosen(i,existing_evoked_ctrl(i,:)==1,:))');
        else
            delta_cdf_ctrl_existing = cat(1,delta_cdf_ctrl_existing,squeeze(delta_cdf_ctrl_chosen(i,existing_evoked_ctrl(i,:)==1,:)));
        end
    end
end


% violin
figure; hold on; 
plot([0 length(theta_vec)+1],[0 0],'--','color',[.5 .5 .5])
data1 = delta_cdf_exp_existing;
data2 = delta_cdf_ctrl_existing;
distributionPlot(gca,data2,'widthDiv',[2 1],'histOri','left','addSpread',true,'color','g','showMM',false,'histOpt',1.1)
distributionPlot(data1,'widthDiv',[2 2],'histOri','right','addSpread',true,'color','r','showMM',false,'histOpt',1.1)
set(gca,'XTickLabel',theta_vec);
ylabel('\Delta cdf'); xlabel('\theta');

figure; hold on; 
vs = violinplot(data1, theta_vec,...
    'HalfViolin','right','ViolinColor',[.64 .08 .18],...%red
    'QuartileStyle','shadow',... % boxplot, none
    'DataStyle', 'scatter',... % scatter, none
    'ShowNotches', false,...
    'ShowMean', false,...
    'ShowMedian', true);
vs = violinplot(data2, theta_vec,...
    'HalfViolin','left','ViolinColor',[.47 .67 .19],...%green
    'QuartileStyle','shadow',... % boxplot, none
    'DataStyle', 'scatter',... % scatter, none
    'ShowNotches', false,...
    'ShowMean', false,...
    'ShowMedian', true);

% significance - experiment vs. control
n_shuf = 1e5;
pval_vec = zeros(size(theta_vec));
for t = 1:length(theta_vec)
    disp(t)
    exp_ = delta_cdf_exp_chosen(:,:,t); exp_ = exp_(:);
    control_ = delta_cdf_ctrl_chosen(:,:,t); control_ = control_(:);
    labels = [ones(size(exp_)); zeros(size(control_))];
    WD_ = [exp_; control_]; 
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