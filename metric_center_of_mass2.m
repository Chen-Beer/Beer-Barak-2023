%% electrode name to 2D index
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

%% save spatio-temporal response & center of mass trajectory

responses_smoothed_2D = cell(n_patterns,1);
T = t2_ - t1_ + 1;
labeled_im = bwlabel(true([a b]));
parfor p = 1:n_patterns
    disp(p)
    responses_ = squeeze(responses_smoothed(p,non_zero_responses(p,:) == 1,:,t1_:t2_));
    R = size(responses_,1);
    responses_smoothed_2D{p} = cell(R,1);
    center_of_mass{p} = zeros(R,2,T);
    for i = 1:R
        response_2D = zeros(a,b,T);
        response = squeeze(responses_(i,:,:));
        for e = relevant_electrodes_idx{p}
            [r,c] = ind2sub([a b],find(contains(electrode_names_ordered,electrodes_names(e))));
            response_2D(r,c,:) = squeeze(response(e,:));
        end
        responses_smoothed_2D{p}{i} = response_2D;
        for t = 1:T
            tmp = squeeze(response_2D(:,:,t));
            props = regionprops(labeled_im, tmp, 'Centroid', 'WeightedCentroid');
            center_of_mass{p}(i,:,t) = props.WeightedCentroid;
        end
    end
end

%% plot mean center of mass trajectory 

time_color_vec = winter(T);
figure;
sgtitle([num2str(t1_),':',num2str(t2_),' [dt = ',num2str(dt),']']);
for p = 1:n_patterns
    if n_patterns > 5
        subplot(q,q,p);
    else
        subplot(1,n_patterns,p);
    end
    hold on;
    [r,c] = ind2sub([a b],find(contains(electrode_names_ordered,electrodes_names(stimulus_electrods_idx{p}))));
    xlim([1 b]); ylim([1 a]);
    plot(c,r,'*r','LineWidth',3,'MarkerSize',5);
    title([num2str(p),' (',num2str(length(find(non_zero_responses(p,:) == 1))),')']);
    % plot all 
    n = size(center_of_mass{p},1);
    for i = 1:n
        cx = squeeze(center_of_mass{p}(i,1,:));
        cy = squeeze(center_of_mass{p}(i,2,:)); 
        plot(cx,cy,'color',[.7 .7 .7],'LineWidth',.5);
        plot(cx(1),cy(1),'o','color',[.7 .7 .7],'LineWidth',.5);    
    end
    cx = nanmean(squeeze(center_of_mass{p}(:,1,:)));
    cy = nanmean(squeeze(center_of_mass{p}(:,2,:)));
%     plot(cx,cy,'k','LineWidth',2);
%     plot(cx(1),cy(1),'ok','LineWidth',1);
    scatter(cx,cy,10,time_color_vec,'filled')
    grid on; drawnow;
end

%% euclidean dist & 2d corr between trajectories - within & across

mean_euclidean_dist_within = cell(n_patterns,1);
mean_euclidean_dist_across = cell(n_patterns,1);
corr_2d_within = cell(n_patterns,1);
corr_2d_across = cell(n_patterns,1);

% across & within
parfor p1 = 1:n_patterns
    disp(p1)
    mean_euclidean_dist_across{p1} = [];
    corr_2d_across{p1} = [];
    for p2 = 1:n_patterns
        R1 = size(center_of_mass{p1},1);
        R2 = size(center_of_mass{p2},1);
        
        if p2 ~= p1 % across
            corr_ = zeros(R1,R2,2);
            norm_dist = zeros(R1,R2);
            for i = 1:R1 
                for j = 1:R2
                    corr_(i,j,1) = corr(squeeze(center_of_mass{p1}(i,1,:)),squeeze(center_of_mass{p2}(j,1,:)),'Rows','complete');
                    corr_(i,j,2) = corr(squeeze(center_of_mass{p1}(i,2,:)),squeeze(center_of_mass{p2}(j,2,:)),'Rows','complete');
                    norm_dist(i,j) = nanmean(vecnorm(squeeze(center_of_mass{p1}(i,:,:))-squeeze(center_of_mass{p2}(j,:,:))));
                end
            end
            corr_2d_across{p1} = cat(2,corr_2d_across{p1},corr_);
            mean_euclidean_dist_across{p1} = cat(2,mean_euclidean_dist_across{p1},norm_dist);
        
        else % within
            corr_2d_within{p1} = zeros(R1,R1,2);
            corr_2d_within{p1}(:,:,1) = corr(squeeze(center_of_mass{p1}(:,1,:))',squeeze(center_of_mass{p1}(:,1,:))','Rows','complete');
            corr_2d_within{p1}(:,:,2) = corr(squeeze(center_of_mass{p1}(:,2,:))',squeeze(center_of_mass{p1}(:,2,:))','Rows','complete');
            norm_dist = zeros(R1,R1);
            for i = 1:R1 
                for j = 1:R1
                    norm_dist(i,j) = nanmean(vecnorm(squeeze(center_of_mass{p1}(i,:,:))-squeeze(center_of_mass{p1}(j,:,:))));
                end
            end
            mean_euclidean_dist_within{p1} = norm_dist;
        end
    end
end

%% plot

% mean euclidean distance
figure;
sgtitle(['mean dist | ',num2str(t1_),':',num2str(t2_),' [dt = ',num2str(dt),']']);
for p = 1:n_patterns
    if n_patterns > 5
        subplot(q,q,p);
    else
        subplot(1,n_patterns,p);
    end
    hold on;
    
    % across
    across = mean_euclidean_dist_across{p};
    histogram(across(:),n_repetitions,'normalization','probability')
    
    % within
    within = mean_euclidean_dist_within{p};
    histogram(within(:),round(n_repetitions),'normalization','probability')
    if p == 1
        legend({'across','within'},'location','best'); 
    end
    title([num2str(p),' (',num2str(length(find(non_zero_responses(p,:) == 1))),')']);
    
end