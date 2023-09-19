clc; clear; close all;

% dir_ = '..\stimulation\age exploration\26544_11.7\exp_4.8_24days\4 spon\';
% dir_ = '..\stimulation\age exploration\broken_28.7\oldexp_5.8_8days\1 spon\';
% dir_ = '..\stimulation\age exploration\broken_28.7\oldexp_5.8_8days\3 spon\';
% dir_ = '..\stimulation\age exploration\broken_28.7\prob_8.8_11days\1 spon\';
% dir_ = '..\stimulation\age exploration\broken_28.7\prob_9.8_12days\1 spon\';
% dir_ = '..\stimulation\age exploration\broken_28.7\prob_10.8_13days\1 spon\';
% dir_ = '..\stimulation\age exploration\broken_28.7\probe_11.8_14days\3 spon\';
% dir_ = '..\stimulation\age exploration\broken_28.7\probe_12.8_15days\1 spon\';
% dir_ = '..\stimulation\age exploration\broken_28.7\probe_13.8_16days\1 spon\';
% dir_ = '..\stimulation\age exploration\broken_28.7\probe_14.8_17days\1 spon\';
% dir_ = '..\stimulation\responsiveness\26536_15.6\2 spontaneous\';
% dir_ = '..\stimulation\responsiveness\broken_20.6\2 spontaneouse (during analysis)\';
% dir_ = '..\stimulation\responsiveness\broken_20.6\2 spontaneouse (during analysis)\';
% dir_ = '..\stimulation\age exploration\26549_9.8\probe_15.8_6days\1 spon\';
% dir_ = '..\stimulation\age exploration\26549_9.8\probe_16.8_7days\1 spon\';
% dir_ = '..\stimulation\age exploration\26549_9.8\probe_17.8_8days\1 spon\';
% dir_ = '..\stimulation\age exploration\26549_9.8\probe_18.8_9days\1 spon\';
% dir_ = '..\stimulation\responsiveness\26536_15.6\2 spontaneous\';
% dir_ = '..\stimulation\age exploration\26537_9.8\probe_19.8_10days\1 spon\';
% dir_= '..\stimulation\age exploration\26537_9.8\probe_20-21.8_11-12days\1 spon\';
% dir_= '..\stimulation\age exploration\26537_9.8\probe_22.8_13days\1 spon\';
% dir_= '..\stimulation\age exploration\26537_9.8\probe_23.8_14days\1 spon\';
% dir_ = '..\stimulation\age exploration\26537_9.8\probe_24.8_15days\1 spon\';
% dir_ = '..\stimulation\age exploration\26537_9.8\probe_25.8_16days\1 spon\';
% dir_ = '..\stimulation\age exploration\26537_9.8\probe_26.8_17days\1 spon\';
% dir_ = '..\stimulation\age exploration\26537_9.8\probe_27.8_18days\1 spon\';
% dir_ = '..\stimulation\age exploration\26537_9.8\probe_28.8_19days\1 spon\';
% dir_ = '..\stimulation\age exploration\26537_9.8\probe_30.8_21days\1 spon\';
% dir_ = '..\stimulation\age exploration\26537_9.8\probe_31.8_22days\1 spon\';
% dir_ = '..\stimulation\age exploration\26532_13.9\probe_30.9_17days\1 spon\';
% dir_ = '..\stimulation\new batch\38427_4.10\exp_24.10_20days\1 spon\';
% dir_ = '..\stimulation\new batch\38427_4.10\exp_24.10_20days\5 spon\';
% dir_ = '..\second protocol\data\38427_4.10\23days_27.10\2 spon\';
% dir_ = '..\second protocol\data\38427_4.10\23days_27.10\4 spon\';
dir_ = '..\second protocol\data\39740_1.11\17days_18.11\1 spont\';

%% load data
fprintf('Loading data...');
load_h5_spontaneous;
disp('Done.');

%% plot raw data
figure; hold on;
for e = 1:N
    tmp = activity_raw{e};
    if ~isempty(tmp)
        plot(tmp,e*ones(size(tmp)),'.k');
    end
end
ylabel('electrode'); xlabel('time [micro-sec]');
xlim([0 T]);
title('raw activity');

%% create smoothed & binned time series

% temporal resolusion
dt = 5e-3; %[sec] 
% kernel
sigma = 2e-2;
kernel_frame = 6; %[*sigma]
units = 1e-6; %[sec]

T_minutes = T*units/60

fprintf('Creating smoothed & binned time series...');
continuous_smoothed_activity;
disp('Done.');

% plot
figure; 
imagesc(activity_smoothed); colorbar;
xlabel(['dt [',num2str(dt),' sec]']);
ylabel('electrode');
% 
% % plot
% figure; hold on
% for e = 1:N
%     plot(edges(1:end-1),e+activity_smoothed(e,:));
% end
% xlabel(['dt [',num2str(dt),' sec]']);
% ylabel('electrode')

%% PR value
fprintf('PR value: ');
e = eig(activity_smoothed*activity_smoothed');
PR_val = ((sum(e))^2)/(sum(e.^2));
disp(PR_val)

%% Detect events - threshold crossing 
overall_activity = squeeze(sum(activity_smoothed));
threshold = 4*std(overall_activity);

% figure; hold on;
% histogram(overall_activity);
% set(gca,'YScale','log')
% plot([threshold threshold],[1 1e7],'-r','linewidth',2);
% title('overall activity');


%% spatial \ temporal profile
figure; 
subplot 221; hold on;
plot(sum(activity_smoothed));
plot([0 size(activity_smoothed,2)],[threshold threshold],'-r','LineWidth',2);
xlabel('dt'); ylabel('overall activity'); title('\Sigma_e'); xlim([0 size(activity_smoothed,2)]);
legend('',[num2str(threshold/std(overall_activity)),'\sigma']);
subplot 222;
plot(sum(activity_smoothed,2)); title('\Sigma_t');
xlabel('electrode'); ylabel('overall activity');
subplot 223; hold on;
histogram(overall_activity);
set(gca,'YScale','log')
plot([threshold threshold],[1 1e7],'-r','linewidth',2);
subplot 224;
histogram(sum(activity_smoothed,2),N/2);


%% PCA

[coeff,~,~,~,explained_var] = pca(activity_smoothed');
activity_pcs = activity_smoothed'*coeff;

figure; hold on;
plot([PR_val PR_val],[0 100],'-k','LineWidth',2);
plot(cumsum(explained_var),'-o');
ylim([0 100]); xlim([1 N])
title('Explained variance');
xlabel('dimensions');
legend(['PR = ',num2str(PR_val)]);
grid on; 

figure; 
hold on;
plot3(activity_pcs(:,1),activity_pcs(:,2),activity_pcs(:,3),'LineWidth',.1)
% above threshold
idx = find(overall_activity > threshold);
plot3(activity_pcs(idx,1),activity_pcs(idx,2),activity_pcs(idx,3),'.r')
grid on; xlabel('PC1');ylabel('PC2');zlabel('PC3');
legend('all','above threshold');
title('activity - PCs');

%% active electrodes - based on raw data
number_of_spikes = zeros(1,N);
for e = 1:N
    number_of_spikes(e) = length(activity_raw{e});
end
figure; 
histogram(number_of_spikes,N);
xlabel('# of spikes - per electrode');
ylabel('electrodes');

%% NS
dt_hist = dt/units;
h =  histogram([activity_raw{:}],0:dt_hist:T);
NS_threshold = threshold*std(h.BinCounts);
figure; hold on; 
time_vec = h.BinEdges(2:end)-h.BinWidth/2;
counts = h.BinCounts;
plot(time_vec,counts,'k');
plot(time_vec,NS_threshold*ones(size(time_vec)),'r','LineWidth',2);
[NS_pks_all,NS_locs_all] = findpeaks(counts,time_vec,'MinPeakDistance',1/units);
[NS_pks,NS_locs] = findpeaks(counts,time_vec,'MinPeakHeight',NS_threshold,'MinPeakDistance',1/units);
NS_dt = diff(NS_locs);

% plot statistics 
figure; 
subplot 121; hold on;
h = histogram(NS_pks_all); xlabel('peak_{NS} [sec]')
plot([NS_threshold NS_threshold],[0 1.1*max(h.BinCounts)],'r','LineWidth',2);
subplot 122;
histogram(NS_dt*units); xlabel('dt_{NS} [sec]');


disp('* * * * * * * * * * * * * * * * * * * * * * *')
disp(['Firing electrodes: ',num2str(sum(number_of_spikes ~= 0)),'\',num2str(N)]);
disp('Network spikes summary:');
disp(['Threshold (# of spikes): ',num2str(NS_threshold)]);
disp(['Median height (# of spikes): ',num2str(median(NS_pks))]);
disp(['Recording duration [min]: ',num2str(T_minutes)]);
disp(['Count - overall: ',num2str(length(NS_pks))]);
disp(['NS/hour: ',num2str(length(NS_pks)*60/T_minutes)]);
disp(['median interval between NSs [sec]: ',num2str(median(NS_dt)*units)]);
disp(['std interval between NSs [sec]: ',num2str(std(NS_dt)*units)]);
disp('* * * * * * * * * * * * * * * * * * * * * * *')

%% animation
% figure; hold on; 
% for i = 1:size(activity_pcs,1)
%    if ismember(i,idx) 
%         plot3(activity_pcs(i,1),activity_pcs(i,2),activity_pcs(i,3),'.r')
%    else
%        plot3(activity_pcs(i,1),activity_pcs(i,2),activity_pcs(i,3),'.k')
%    end
%     title(i);
%     drawnow;
% end
% grid on; xlabel('PC1');ylabel('PC2');zlabel('PC3');
% legend('all','above threshold');
% title('activity - PCs');

%% Center of mass trajectory - 2D covarage

% electrode name to 2D index
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
labeled_im = bwlabel(true([a b]));

T_ = 30*60/dt;
% center of mass trajectory - 30min
activity_2D = zeros(a,b,T_);
CoM_trajectory = zeros(2,T_);
for e = 1:N
    [r,c] = ind2sub([a b],find(contains(electrode_names_ordered,electrodes_names(e))));
    activity_2D(r,c,:) = squeeze(activity_smoothed(e,end-T_+1:end));
end
parfor t = 1:T_
    tmp = squeeze(activity_2D(:,:,t));
    props = regionprops(labeled_im, tmp, 'Centroid', 'WeightedCentroid');
    CoM_trajectory(:,t) = props.WeightedCentroid;
end


% plot 
figure; hold on;
title('center of mass trajectory - 30min');
grid on; xlim([1 b]); ylim([1 a]);
cx = CoM_trajectory(1,:);
cy = CoM_trajectory(2,:);
plot(cx,cy,'.k','MarkerSize',.1);
% plot(cx(1),cy(1),'.k','LineWidth',.5);
% for t = 2:T_
%     plot(cx(t),cy(t),'.k','LineWidth',.5);
%     drawnow;
% end

%% save .mat
fprintf('Saving mat...');
save([dir_,'data.mat'], '-v7.3')

disp('Done.');