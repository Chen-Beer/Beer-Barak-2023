%% events\NS detection
% input: a countinuous time series of activity
% assumption: input data is taken from ~the same day
% output: time frames [t1 t2] of all global events (\NS) & croped smoothed activity (continuous) & center of mass

function [time_frames,events_smoothed,events_2D,events_CoM] = event_detection(activity_smoothed_,dt,electrode_names_ordered,electrodes_names,chosen_electrodes,plot_flag,clean_spont_events)
    
    fprintf('Extracting events...');
    high_threshold = 4; %[*std]
    low_threshold = .5;
    min_dist = 3/dt;
    min_crossing_low_threshold = .5/dt;
    
    if clean_spont_events
        min_event_length = 0.15/dt;
        max_event_length = 0.5/dt;
        % -- evoked length is 0.35/dt in [dt] units
    end

%     dt_hist = dt/units;
%     figure; 
%     h =  histogram(activity_raw,0:dt_hist:T);
%     NS_threshold = threshold*std(h.BinCounts);
%     time_vec = h.BinEdges(2:end)-h.BinWidth/2;
%     counts = h.BinCounts;
%     clf; hold on;
%     plot(time_vec,counts,'k');
%     plot(time_vec,NS_threshold*ones(size(time_vec)),'r','LineWidth',2);
%     findpeaks(counts,time_vec,'MinPeakHeight',NS_threshold,'MinPeakDistance',2/units,'Annotate','extents');
%     [NS_pks,NS_locs] = findpeaks(counts,time_vec,'MinPeakHeight',NS_threshold,'MinPeakDistance',2/units);
    
    activity = sum(activity_smoothed_);

     
    NS_threshold = high_threshold*std(activity);
    baseline_threshold = low_threshold*std(activity);
    if plot_flag 
        figure;
        plot(activity); hold on;
        plot([0 length(activity)],[NS_threshold NS_threshold],'r','LineWidth',2);
        plot([0 length(activity)],[baseline_threshold baseline_threshold],'k','LineWidth',2);
        findpeaks(activity,'MinPeakHeight',NS_threshold,'MinPeakDistance',min_dist);
    end
    [NS_pks,NS_locs] = findpeaks(activity,'MinPeakHeight',NS_threshold,'MinPeakDistance',min_dist);

    % time frames
    % t1 - up crossing baseline_threshold
    % t2 - down crossing baseline_threshold for at least min_crossing_low_threshold
    time_frames = zeros(length(NS_pks),1);
    for i = 1:length(NS_pks)
        tt = NS_locs(i);
        t1 = tt;
        while t1 > 1 && activity(t1-1) >= baseline_threshold
            t1 = t1 - 1;
        end
        
        t2 = tt;
        while activity(t2+1) >= baseline_threshold
            t2 = t2 + 1;
        end
        while (t2 + min_crossing_low_threshold) < size(activity_smoothed_,2) && activity(t2 + min_crossing_low_threshold) > baseline_threshold
            t2 = t2 + min_crossing_low_threshold;
            while activity(t2+1) >= baseline_threshold
                t2 = t2 + 1;
            end
        end
        if t2 > length(activity)
            t2 = length(activity);
        end
        time_frames(i,1) = t1;
        time_frames(i,2) = t2;
        
        if plot_flag 
            % plot
            plot([t1 t1],[0 max(activity)],'-k');
            plot([t2 t2],[0 max(activity)],'-k');
        end
    end
    
    % remove overlapping events
    i = 2;
    while i < size(time_frames,1)
        if time_frames(i,1) < time_frames(i-1,2) && time_frames(i,2) <= time_frames(i-1,2)
            time_frames(i,:) = [];
        else
            i = i+1;
        end
    end
    
    % clean too short\long events 
    if clean_spont_events
        i = 1;
        while i < size(time_frames,1)
            if (time_frames(i,2)- time_frames(i,1) < min_event_length) ||  (time_frames(i,2)- time_frames(i,1) > max_event_length)
                time_frames(i,:) = [];
            else
                i = i+1;
            end    
        end
    end
    
    if plot_flag 
        % plot final events - overall activity
        for i = 1:length(time_frames)
            t1 = time_frames(i,1);
            t2 = time_frames(i,2);
            plot([t1 t1],[0 max(activity)],'-g');
            plot([t2 t2],[0 max(activity)],'-g');
        end
    end
    
    % crop events 
    events_smoothed = cell(size(time_frames,1),1);
    for i = 1:size(time_frames,1)
        t1 = time_frames(i,1);
        t2 = time_frames(i,2);
        events_smoothed{i} = activity_smoothed_(:,t1:t2);
    end
    
    if plot_flag 
        % plot events - per electrode
        q = ceil(sqrt(size(time_frames,1)));
        figure; 
        idx = 1;
        for i = 1:size(time_frames,1)
            t1 = time_frames(i,1);
            t2 = time_frames(i,2);
            subplot (q,2*q,idx); plot(activity(t1:t2));
            subplot (q,2*q,idx+1); imagesc(events_smoothed{i});
            idx = idx + 2;
        end
        sgtitle('all global events')
    end
    
    % center of mass
    a = 10; b = 12;
    E = length(chosen_electrodes);
    events_2D = cell(size(time_frames,1),1);
    events_CoM = cell(size(time_frames,1),1);
    labeled_im = bwlabel(true([a b]));
    parfor i = 1:size(time_frames,1)    
        event_ = squeeze(events_smoothed{i});
        T = size(events_smoothed{i},2);
        event_2D = zeros(a,b,T);
        for e = 1:E
            [r,c] = ind2sub([a b],find(contains(electrode_names_ordered,electrodes_names(chosen_electrodes(e)))));
            event_2D(r,c,:) = squeeze(event_(e,:));
        end
        events_2D{i} = event_2D;
        events_CoM{i} = zeros(2,T);
        for t = 1:T
            tmp = squeeze(event_2D(:,:,t));
            props = regionprops(labeled_im, tmp, 'Centroid', 'WeightedCentroid');
            events_CoM{i}(:,t) = props.WeightedCentroid;
        end
    end
 
    disp('Done.')
    disp([num2str(size(time_frames,1)),' events in total']);
end