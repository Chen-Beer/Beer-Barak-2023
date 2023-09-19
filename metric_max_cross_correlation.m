T = length(t1_:t2_);

% within
cc_score_within = cell(n_patterns,1);
parfor p = 1:n_patterns
    disp(p)
    responses_ = squeeze(responses_smoothed(p,non_zero_responses(p,:) == 1,relevant_electrodes_idx{p},t1_:t2_));
    R = size(responses_,1);
    E = size(responses_,2);
    cc_score_within{p} = zeros(R,R);
    for j = 1:R
               
        for k = 1:R
            cc = zeros(E,2*T-1);
            for e = 1:E
                tmp = squeeze(responses_(:,e,:));
                [c, ~] = xcorr(tmp(j,:),tmp(k,:),'coef');
                cc(e,:) = c;
            end
            cc_score_within{p}(j,k) = max(nanmean(cc));
        end
        
    end
end

% across 
cc_score_across = cell(n_patterns,1);
parfor p1 = 1:n_patterns
    disp(p1)
    electrodes1 = relevant_electrodes_idx{p1};
    cc_score_across{p1} = [];
    for p2 = 1:n_patterns
        if p2 ~= p1
            electrodes2 = relevant_electrodes_idx{p2};
            ee = intersect(electrodes1,electrodes2);
            reps1 = find(non_zero_responses(p1,:) == 1);
            reps2 = find(non_zero_responses(p2,:) == 1);
            responses1 = squeeze(responses_smoothed(p1,reps1,ee,t1_:t2_));
            responses2 = squeeze(responses_smoothed(p2,reps2,ee,t1_:t2_));
            R1 = size(responses1,1);
            R2 = size(responses2,1);
            corr_ = zeros(R1,R2);
            for i = 1:R1 
                for j = 1:R2
                    cc = nan(length(ee),2*T-1);
                    if length(ee) == 1
                        current = squeeze(responses1(i,:));
                        other = squeeze(responses2(j,:));
                        cc = xcorr(current,other,'coef');
                    else 
                        for e = 1:length(ee)
                            current = squeeze(responses1(i,e,:));
                            other = squeeze(responses2(j,e,:));
                            [c, ~] = xcorr(current,other,'coef');
                            cc(e,:) = c;
                        end
                    end
                    corr_(i,j) = max(nanmean(cc));
                end
            end
            cc_score_across{p1} = cat(2,cc_score_across{p1},corr_);
        end
        
    end
end

%% plot xcorr within & across distributions

vals = 0:0.01:1;

figure;
sgtitle([num2str(t1_),':',num2str(t2_),' [dt = ',num2str(dt),']']);
for p = 1:n_patterns
    if n_patterns > 5
        subplot(q,q,p);
    else
        subplot(1,n_patterns,p);
    end
    hold on;
    
    % across
    across = cc_score_across{p};
    histogram(across(:),n_repetitions,'normalization','probability')
    
    % within
    within = cc_score_within{p};
    histogram(within(:),round(n_repetitions),'normalization','probability')
    
    if p == 1
        legend({'across','within'},'location','best'); 
    end
    title([num2str(p),' (',num2str(length(find(non_zero_responses(p,:) == 1))),')']);
    
end