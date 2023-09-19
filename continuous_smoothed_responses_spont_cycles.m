%% create smoothed & binned time series

% Responses

nn = min(size(responses_raw{2},1),size(responses_raw{1},1));
[ii,pp] = groupcounts(sequence');
zero_rep = max(ii);
less_frequent = pp(end);
pp(less_frequent) = [];

n_bins = round((T_response-pre/units)*units/dt);
responses_smoothed = zeros(n_patterns,nn,N,n_bins);
for p = pp'
%     disp(p)
    parfor e = 1:N
        disp([p, e])
        for i = 1:nn
            t1 = double(responses_time_frames{p}(i,1))*units; %first pulse
            t2 = double(responses_time_frames{p}(i,2))*units; %second pulse
            t2 = t2 - pre; % remove pre-pulse delay 
%             t2 = t1 + double((T_response-pre/units)*units);
            edges = t1:dt:t2-dt;
            if length(edges)+1 == n_bins 
                edges = t1:dt:t2;
            end
            if length(edges)-1 == n_bins 
                edges = t1:dt:t2-2*dt;
            end
            edges = edges - t1;
            tmp = double(responses_raw{p}{i,e})*units - t1;
            if ~isempty(tmp)
                binned_response = zeros(length(edges),1);
                for s = 1:length(tmp)
                    mu = tmp(s);
                    if mu < t2
                        fun = @(x) (1/(sqrt(2*pi)*sigma))*exp((-(x-mu).^2)./(2*sigma.^2));
                        tt1 = mu - kernel_frame*sigma;
                        tt2 = mu + kernel_frame*sigma;
                        [~,tt1_idx] = min(abs(edges - tt1));
                        [~,tt2_idx] = min(abs(edges - tt2));
                        for t = tt1_idx:tt2_idx-1
                            a = edges(t);
                            b = edges(t+1);
                            res = integral(fun,a,b);
                            binned_response(t) = binned_response(t) + res;
                        end
                    end
                end
                responses_smoothed(p,i,e,:) = binned_response;  
%                 figure; hold on;plot(edges,binned_response);plot(tmp,ones(size(tmp)),'x');xlabel('Time [sec]');legend('smoothed FR','spikes');title(['electrode ',electrodes_names{e},' pattern #',num2str(p),' repetition #',num2str(i)]);
            end
        end
    end
end

p = less_frequent;
parfor e = 1:N
    disp([p, e])
    ii_ = 0;
    for i = 1:nn
        if mod(i,zero_rep) ~= 0
            ii_ = ii_ + 1;
            t1 = double(responses_time_frames{p}(ii_,1))*units; %first pulse
            t2 = double(responses_time_frames{p}(ii_,2))*units; %second pulse
            t2 = t2 - pre; % remove pre-pulse delay 
    %             t2 = t1 + double((T_response-pre/units)*units);
            edges = t1:dt:t2-dt;
            if length(edges)+1 == n_bins 
                edges = t1:dt:t2;
            end
            if length(edges)-1 == n_bins 
                edges = t1:dt:t2-2*dt;
            end
            edges = edges - t1;
            tmp = double(responses_raw{p}{ii_,e})*units - t1;
            if ~isempty(tmp)
                binned_response = zeros(length(edges),1);
            for s = 1:length(tmp)
                mu = tmp(s);
                if mu < t2
                    fun = @(x) (1/(sqrt(2*pi)*sigma))*exp((-(x-mu).^2)./(2*sigma.^2));
                    tt1 = mu - kernel_frame*sigma;
                    tt2 = mu + kernel_frame*sigma;
                    [~,tt1_idx] = min(abs(edges - tt1));
                    [~,tt2_idx] = min(abs(edges - tt2));
                    for t = tt1_idx:tt2_idx-1
                        a = edges(t);
                        b = edges(t+1);
                        res = integral(fun,a,b);
                        binned_response(t) = binned_response(t) + res;
                    end
                end
            end
            responses_smoothed(p,i,e,:) = binned_response;  
%                 figure; hold on;plot(edges,binned_response);plot(tmp,ones(size(tmp)),'x');xlabel('Time [sec]');legend('smoothed FR','spikes');title(['electrode ',electrodes_names{e},' pattern #',num2str(p),' repetition #',num2str(i)]);
            end
        end
    end
end

% Spontaneous activity
T_spont = double(spont_off(end));
n_bins = round(T_spont*units/dt);
spont_activity_smoothed = zeros(N,n_bins);
edges = 0:dt:(T_spont*units);
if length(edges) == n_bins + 1
    edges = edges(1:end-1);
end

parfor e = 1:N
    disp(e)
    tmp = double(spont_activity_raw{e})*units;
    if ~isempty(tmp)
        binned_response = zeros(length(edges),1);
        for s = 1:length(tmp)
            mu = tmp(s);
            fun = @(x) (1/(sqrt(2*pi)*sigma))*exp((-(x-mu).^2)./(2*sigma.^2));
            tt1 = mu - kernel_frame*sigma;
            tt2 = mu + kernel_frame*sigma;
            [~,tt1_idx] = min(abs(edges - tt1));
            [~,tt2_idx] = min(abs(edges - tt2));
            for t = tt1_idx:tt2_idx-1
                a = edges(t);
                b = edges(t+1);
                res = integral(fun,a,b);
                binned_response(t) = binned_response(t) + res;
            end
        end
        spont_activity_smoothed(e,:) = binned_response;  
%                 figure; hold on;plot(edges,binned_response);plot(tmp,ones(size(tmp)),'x');xlabel('Time [sec]');legend('smoothed FR','spikes');title(['electrode ',electrodes_names{e},' pattern #',num2str(p),' repetition #',num2str(i)]);
    end
end