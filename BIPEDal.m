participants = {'P01','P02','P03','P04','P05','P07','P10','P11','P12','P14'};
sets = {'s1','s2','s3'};
muscles = {'BFl','BFr','GAl','GAr','RFl','RFr','TAl','TAr','VLl','VLr','VMl','VMr'};
fs = 1500;

[b_high, a_high] = butter(2, 30/(fs/2), 'high');
[b_low, a_low] = butter(2, 10/(fs/2), 'low');

target_len = 2900;

for p = 1:length(participants)
    pid = participants{p};
    
    for s = 1:3
        set_name = sets{s};
        dataMat = [];
        
        for m = 1:length(muscles)
            varname = sprintf('%s_%s', muscles{m}, set_name);
            raw = eval(varname);  % e.g., BFl_s1
            
            if isfield(raw, pid)
                emg = raw.(pid);  % get signal for that participant
                emg = double(emg(:))';  % ensure row vector
                
                % Filter and preprocess
                emg = filtfilt(b_high, a_high, emg);
                emg = abs(hilbert(emg));  % rectify
                emg = filtfilt(b_low, a_low, emg);  % envelope
                emg = emg / max(abs(emg));  % max normalization
                emg = emg / std(emg);  % unit variance
                
                emg = resample(emg, target_len, length(emg));  % time normalize
                dataMat = [dataMat; emg];  % append muscle × time
            else
                warning('%s missing in %s', pid, varname);
                dataMat = [dataMat; zeros(1, target_len)];
            end
        end
        
        % NNMF
        max_synergies = 10;
        reps = 50;
        maxIter = 1000;
        RA_thresh = 0.9;
        RA_delta = 0.03;
        RA_list = zeros(1, max_synergies);

        for k = 1:max_synergies
            RAs = zeros(1, reps);
            for r = 1:reps
                [W, H] = nnmf(dataMat, k, 'replicates', 1, ...
                    'options', statset('MaxIter', maxIter, 'Display', 'off'));
                recon = W * H;
                err = norm(dataMat - recon, 'fro') / norm(dataMat, 'fro');
                RAs(r) = 1 - err;
            end
            RA_list(k) = mean(RAs);
            if k > 1 && RA_list(k) > RA_thresh && (RA_list(k)-RA_list(k-1)) < RA_delta
                break;
            end
        end
        
        final_k = k;
        [W_final, H_final] = nnmf(dataMat, final_k, 'replicates', reps, ...
            'options', statset('MaxIter', maxIter));
        
        % Normalize W
        for i = 1:final_k
            normF = norm(W_final(:,i));
            W_final(:,i) = W_final(:,i) / normF;
            H_final(i,:) = H_final(i,:) * normF;
        end
        
        synergyResult(p).participant = pid;
        synergyResult(p).set(s).W = W_final;
        synergyResult(p).set(s).H = H_final;
        synergyResult(p).set(s).RA = RA_list(final_k);
        synergyResult(p).set(s).k = final_k;
    end
end

%% Graphs
colors = lines(3);  % One color per load/set
muscle_labels = {'BFl','BFr','GAl','GAr','RFl','RFr','TAl','TAr','VLl','VLr','VMl','VMr'};

% Choose a participant to plot
pid_idx = 1;
pid = synergyResult(pid_idx).participant;

for set_id = 1:3
    W = synergyResult(pid_idx).set(set_id).W;
    H = synergyResult(pid_idx).set(set_id).H;
    k = synergyResult(pid_idx).set(set_id).k;

    figure('Name', sprintf('%s - Set %d', pid, set_id), 'Color', 'w');
    
    for s = 1:5
        % Plot waveform (H)
        subplot(2, 5, s)
        plot(H(s,:), 'Color', colors(set_id,:), 'LineWidth', 2)
        title(sprintf('Synergy %d - H', s))
        xlabel('Time (normalized)')
        ylabel('Activation')
        xlim([1 size(H,2)])
        grid on

        % Plot weight (W)
        subplot(2, 5, 5 + s)
        bar(W(:,s), 'FaceColor', colors(set_id,:))
        title(sprintf('Synergy %d - W', s))
        ylabel('Weight')
        xticks(1:length(muscle_labels))
        xticklabels(muscle_labels)
        xtickangle(45)
        ylim([0 1])
        grid on
    end
end

%% Temporal activation patterns
% Select participant index (e.g., P01 is index 1)
pid_idx = 1;
pid = synergyResult(pid_idx).participant;

% Number of synergies from any set (assuming same across sets)
k = synergyResult(pid_idx).set(1).k;
colors = [0 0.4470 0.7410;   % blue (Set 1)
          0.8500 0.3250 0.0980;  % red (Set 2)
          0.4660 0.6740 0.1880]; % green (Set 3)
set_labels = {'50%', '62.5%', '75%'};

figure('Name', ['Synergy Waveforms - Participant ' pid], 'Color', 'w');

for s = 1:k
    subplot(ceil(k/2), 2, s);  % adjust layout for clarity

    hold on;
    for set_id = 1:3
        H = synergyResult(pid_idx).set(set_id).H;
        if s <= size(H, 1)
            plot(H(s,:), 'LineWidth', 2, 'Color', colors(set_id, :));
        end
    end
    hold off;

    title(['Synergy ' num2str(s)]);
    xlabel('Normalized Time');
    ylabel('Activation');
    legend(set_labels, 'Location', 'northeast');
    grid on;
    xlim([1 size(H,2)]);
end

%% Weights

participants = 1:length(synergyResult);
muscles = {'BFl','BFr','GAl','GAr','RFl','RFr','TAl','TAr','VLl','VLr','VMl','VMr'};
numMuscles = length(muscles);
numSets = 3;

% Initialize storage: muscle x set x participant
weights = nan(numMuscles, numSets, length(participants));

% Collect Synergy 1 weights for all
for p = 1:length(participants)
    for s = 1:numSets
        W = synergyResult(p).set(s).W;
        if size(W,2) >= 1
            weights(:, s, p) = W(:,1);  % Synergy 1 weights
        end
    end
end

% Compute mean and std across participants
meanWeights = mean(weights, 3, 'omitnan');  % muscle x set
stdWeights = std(weights, 0, 3, 'omitnan');

% Plot
figure('Color','w');
hold on;
b = bar(meanWeights');  % Transpose to get bars grouped by muscle

% Colors matching the paper (blue, red, black)
b(1).FaceColor = [0 0.4470 0.7410];   % 50% - blue
b(2).FaceColor = [0.8500 0.3250 0.0980];  % 62.5% - red
b(3).FaceColor = [0.2 0.2 0.2];       % 75% - dark gray/black

% Add error bars
ngroups = size(meanWeights, 1);
nbars = size(meanWeights, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, meanWeights(:,i), stdWeights(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1);
end

% Labels
set(gca, 'XTick', 1:numMuscles, 'XTickLabel', muscles);
xtickangle(45);
ylabel('Weight contribution');
ylim([0 1]);
legend({'50%', '62.5%', '75%'}, 'Location', 'northoutside', 'Orientation', 'horizontal');
title('Synergy 1 Weight Contributions Across Loads');
grid on;

%%


