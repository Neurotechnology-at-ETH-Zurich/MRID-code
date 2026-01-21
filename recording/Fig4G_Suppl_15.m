%% Magnetic Resonance Identification Tags for Ultra-Flexible Electrodes
%  This script loads spike-sorted single-unit activity from rat rEO_10 and
%  quantifies peri–sharp-wave ripple (SWR) population dynamics. If a
%  previously computed summary structure (Cells) is available
%  (Fig4G_Suppl_15.mat), it is loaded to skip the full preprocessing; 
%  otherwise the script:
%
%   1) Loads spike sorting outputs (spike times, cluster assignments, sites,
%      waveforms, SNR, etc.) together with curated SWR event timestamps.
%   2) Loads probe mapping/metadata from the .prm file to obtain channel
%      geometry, scaling, and sampling rate.
%   3) Assigns each cluster a putative anatomical structure and hemisphere
%      based on external channel-to-structure mapping text files (right and
%      left hemisphere maps).
%   4) For each well-sampled unit (>=100 spikes), computes basic spike-train
%      descriptors (autocorrelogram, ISI histogram) and builds per-ripple
%      peri-event spike histograms in fixed time windows around ripple
%      centers. Peri-event responses are Gaussian-smoothed and converted to
%      z-scored kernel-smoothed firing rate traces.
%   5) Builds a population activity matrix (neurons × time) and visualizes
%      it as a heatmap for the left hemisphere (Fig. 4G / Suppl. Fig. 15A).
%   6) Runs PCA on population responses around SWRs, estimates the optimal
%      k-means cluster number using the Calinski–Harabasz criterion, and
%      sorts the heatmap by cluster assignment.
%   7) Performs cluster-level statistics by comparing mean activity during
%      a ripple-centered window versus a spontaneous baseline window using
%      two-sample t-tests (alpha = 0.01), and visualizes distributions with
%      mean±STD and violin-style plots (Suppl. Fig. 15B).

% Author: Peter Gombkoto, Ph.D.
% Email: pgombkoto@ethz.ch

%% =========================================================================
%  Rat identifier
% =========================================================================
ratID = 'rEO_10';

hasStatSummary = false;   % alapállapot

try
    load('Fig4G_Suppl_15.mat')
    CellIDs = [Cells.CellID_original];
    hasStatSummary = true;
catch
    hasStatSummary = false;
end

%% =========================================================================
%  SPIKE DATA LOADING
% =========================================================================

if  hasStatSummary == 0
    folder = ['G:\Yanik Lab Dropbox\Peter Gombkoto\Localization Manuscript 2024\' ...
        'RAT DATA\rEO_10\recording\05_241022_134111_selected\' ...
        'New_Sorting_for_the_paper\05_241022_134111'];

    cd(folder)

    % Load spike sorting results
    load(fullfile(folder,'amplifier_res.mat'), 'spikeTimes')
    load(fullfile(folder,'amplifier_res.mat'), 'spikesByCluster')
    load(fullfile(folder,'amplifier_res.mat'), 'spikeSites')
    load(fullfile(folder,'amplifier_res.mat'), 'spikeClusters')
    load(fullfile(folder,'amplifier_res.mat'), 'clusterSites')
    load(fullfile(folder,'amplifier_res.mat'), 'clusterNotes')
    load(fullfile(folder,'amplifier_res.mat'), 'unitVppRaw')
    load(fullfile(folder,'amplifier_res.mat'), 'unitSNR')
    load(fullfile(folder,'amplifier_res.mat'), 'unitCount')
    load(fullfile(folder,'amplifier_res.mat'), 'meanWfLocal')
    load(fullfile(folder,'amplifier_res.mat'), 'filtShape')

    % Load curated ripple events
    load(fullfile(folder,'curated_ripples.mat'))

    CellsID=unique((spikeClusters(spikeClusters>0)));

    %% =========================================================================
    %  PROBE MAP AND RECORDING METADATA
    % =========================================================================
    [probePad, shankMap, siteLoc, siteMap, bitScaling, nChans, sample_rate] = ...
        loadProabMap(fullfile(folder,'amplifier.prm'));

    filename = 'amplifier.dat';
    FileInfo = dir(fullfile(folder, filename));
    nChansInRawFile = nChans;
    num_samples = FileInfo.bytes / (nChansInRawFile * 2);   % int16 → 2 bytes

    %% =========================================================================
    %  PUTATIVE ANATOMICAL STRUCTURE ASSIGNMENT
    % =========================================================================
    PutitativeStructure = {};

    % -------------------- Right hemisphere --------------------
    filename = 'rEO_10-CHmap-right_hemisphere-duo.txt';
    [siteMap_Nano, StructureName] = importEminhanTXT(filename);

    for n_cluster_Site = 1:length(clusterSites)
        PutitativeStructure(n_cluster_Site,1) = ...
            {char(StructureName(siteMap_Nano == siteMap(clusterSites(n_cluster_Site)) - 1))};
        Sites_Neuorscope(n_cluster_Site) = siteMap(clusterSites(n_cluster_Site)) - 1;
        PutitativeStructure(n_cluster_Site,2) = {'Right hem.'};
    end

    % -------------------- Left hemisphere --------------------
    filename = 'rEO_10-CHmap-left_hemisphere-penta.txt';
    [siteMap_Nano, StructureName] = importEminhanTXT(filename);

    siteMap_Nano = siteMap_Nano + 64;

    for n_cluster_Site = 1:length(clusterSites(1:max(CellsID)))
        if ~isempty(StructureName(siteMap_Nano == siteMap(clusterSites(n_cluster_Site)) - 1))
            PutitativeStructure(n_cluster_Site,1) = ...
                {char(StructureName(siteMap_Nano == siteMap(clusterSites(n_cluster_Site)) - 1))};
            Sites_Neuorscope(n_cluster_Site) = siteMap(clusterSites(n_cluster_Site)) - 1;
            PutitativeStructure(n_cluster_Site,2) = {'Left hem.'};
        end
    end
end

%% =========================================================================
%  KERNELS FOR SMOOTHING FIRING RATES
% =========================================================================

sigma_for_trial = 0.005;
bin_trial = 0.001;

edges_trial = -0.5 : 0.015 : 0.5;
edges_norm_for_trial = -0.5*sigma_for_trial : bin_trial : 0.5*sigma_for_trial;

kernel_trial = normpdf(edges_norm_for_trial, 0, sigma_for_trial);
kernel_trial = kernel_trial * bin_trial;

sigma = 0.005;
edges_norm = -0.5*sigma : 0.001 : 0.5*sigma;

kernel = normpdf(edges_norm, 0, sigma);
kernel = kernel * 0.001;

edges       = -0.5 : 0.005 : 0.5;   % 5 ms bins
edges_large = -0.5 : 0.025 : 0.5;   % 25 ms bins
bin = 0.01;

%% =========================================================================
%  CONSTANTS AND INITIALIZATION
% =========================================================================
if hasStatSummary == 0
    spikeSites    = double(clusterSites);
    spikeTimes    = double(spikeTimes);
    spikeClusters = spikeClusters;

    a = memmapfile('amplifier_filt.jrc', 'Format', 'int16');
    spikesFilt_unshaped = a.Data;

    waveform = meanWfLocal;
    Recording_CellID = unique(spikeClusters(spikeClusters > 0));

    WfLocalRaw = double(reshape(spikesFilt_unshaped, filtShape)) .* bitScaling;

    halfdeltaT = 0.5;
    y1 = logspace(0, 10, 100);

    CellID_counter = 1;

    %% =========================================================================
    %  MAIN LOOP OVER CELLS
    % =========================================================================
    for CellID = 1:length(Recording_CellID)

        centers  = sum(ripple_timestamps,2) / 2;
        time_on  = centers - halfdeltaT;
        time_off = centers + halfdeltaT;

        fprintf('Processing CellID %d\n', CellID);

        SpikeSecund = spikeTimes(spikeClusters == CellID) ./ sample_rate;

        if length(SpikeSecund) > 100

            % Autocorrelogram
            [tR, ccgR] = spike_autocorrelation(SpikeSecund, 0.4, 0.002);
            binCCG = 0.002 * 1000;

            % Inter-spike interval histogram
            sejt_ido_isi = SpikeSecund * 1000;
            ISI = histcounts(diff(sejt_ido_isi), y1);
            xlimitation = [0 mean(diff(sejt_ido_isi)) + 1000];

            % Allocate matrices
            spikes_trials = zeros(length(centers), length(edges));
            spikes_trials_assembly = zeros(length(centers), length(edges_trial));
            spikes_trials_assembly_Gauss = zeros(length(centers), length(edges_trial));

            trial_ripple_ON = 1;

            for trial = 1:length(centers)

                idx = SpikeSecund >= time_on(trial) & SpikeSecund <= time_off(trial);

                if any(idx)

                    spike(trial_ripple_ON).DuringRippleTime = ...
                        SpikeSecund(idx) - centers(trial);

                    spike(trial_ripple_ON).time = ...
                        SpikeSecund(idx) - centers(trial);

                    spike(trial_ripple_ON).counts_pre  = ...
                        sum(spike(trial_ripple_ON).time >= -0.5 & spike(trial_ripple_ON).time < -0.3);

                    spike(trial_ripple_ON).counts_post = ...
                        sum(spike(trial_ripple_ON).time > -0.1 & spike(trial_ripple_ON).time <= 0.1);

                    spike(trial_ripple_ON).ripple_classes = char(ripple_classes(trial));

                    spikes_trials(trial,:) = histc(spike(trial_ripple_ON).time, edges)';
                    spikes_trials_assembly(trial,:) = histc(spike(trial_ripple_ON).time, edges_trial)';

                    s = conv(spikes_trials_assembly(trial,:), kernel_trial);
                    center_gauss = ceil(length(edges_trial)/2);

                    spikes_trials_assembly_Gauss(trial,:) = ...
                        s(ceil(length(s)/2)-(center_gauss-1) : ceil(length(s)/2)+(center_gauss-1));

                    trial_ripple_ON = trial_ripple_ON + 1;
                end
            end

            % Z-scored kernel-smoothed population response
            temp_conv = conv(zscore(mean(spikes_trials,1)), kernel);

            Cells(CellID_counter).spikes_kernel_mean = temp_conv(3:end-3);
            Cells(CellID_counter).individual_trials = spikes_trials;
            Cells(CellID_counter).individual_trials_count_assembly = spikes_trials_assembly;
            Cells(CellID_counter).individual_trials_gauss_assembly = spikes_trials_assembly_Gauss;
            Cells(CellID_counter).ccgR = ccgR;
            Cells(CellID_counter).tR = tR;
            Cells(CellID_counter).ISI = ISI;
            Cells(CellID_counter).ISI_lim = xlimitation;
            Cells(CellID_counter).Time = edges(1:end-1);
            Cells(CellID_counter).Structre = PutitativeStructure(CellID,1);
            Cells(CellID_counter).Hemisphare = PutitativeStructure(CellID,2);
            Cells(CellID_counter).CellID_original = CellID;
            Cells(CellID_counter).spikeSites = spikeSites(CellID);
            Cells(CellID_counter).wavform = waveform(:,1,CellID);
            Cells(CellID_counter).NeuroscopeID = Sites_Neuorscope(CellID);
            Cells(CellID_counter).ratID = ratID;

            CellID_counter = CellID_counter + 1;
        end
    end

    %% =========================================================================
    %  POPULATION ACTIVITY HEATMAP (LEFT HEMISPHERE)
    % =========================================================================

    % Select neurons recorded in the left hemisphere
    num_cell = find(ismember([Cells.Hemisphare], {'Left hem.'}));
    CellIDs  = num_cell;

    % Extract unique anatomical structure labels
    unique_Structure = unique([Cells.Structre]);

    % Assemble z-scored, kernel-smoothed firing rate matrix
    spikes_kernel_mean = [Cells(CellIDs).spikes_kernel_mean];

    % Neurons × Time matrix for visualization and PCA
    brain_struct = reshape( ...
        spikes_kernel_mean, ...
        length(Cells(CellIDs(1)).spikes_kernel_mean), ...
        length(CellIDs) )';

    %% =========================================================================
    %  HEATMAP OF Z-SCORED FIRING RATES DURING SWRs
    % =========================================================================
end
figure();

ax1 = subplot(1,2,1);
imagesc(brain_struct);
set(gca, 'YDir', 'normal')          % Dorsoventral axis: bottom → top

% Y-axis: neuron indices labeled by anatomical structure
yticks(1:length(CellIDs))
labelsforleft = [Cells.Structre];
yticklabels(labelsforleft(CellIDs))

% Colormap and scaling
colormap(CustomColormap)
clim([-2 2])

% Title and axes
title({'Figure 4G – Z-scored neuronal firing rates', ...
    'along the dorsoventral axis during SWRs'})

xticks([1, ...
    (length(Cells(1).Time)/2)+1, ...
    length(Cells(1).Time)])

xticklabels(num2cell( ...
    round(Cells(1).Time([1, ...
    (length(Cells(1).Time)/2)+1, ...
    length(Cells(1).Time)]), 2) * 1000))

xlabel('Time (ms)')
ylabel('Neurons')

% Colorbar
colorbar(ax1, 'southoutside')
hold off

%% =========================================================================
%  PCA ON Z-SCORED FIRING RATES (POPULATION DYNAMICS)
% =========================================================================
if  hasStatSummary == 0
    % Define time window around ripple for PCA (in ms)
    time_around_stimulus = [ ...
        find(round(Cells(1).Time * 1000) == -150), ...
        find(round(Cells(1).Time * 1000) ==  200) ];

    % Perform PCA on population activity
    [coeff, score, ~, ~, explained, ~] = ...
        pca(brain_struct(:, ...
        time_around_stimulus(1):time_around_stimulus(2)));

    % Retain components explaining <90% cumulative variance
    PCA_best = find(cumsum(explained) < 90);

    opts = statset('Display','final');

    % Estimate optimal cluster number using Calinski–Harabasz criterion
    try
        eva = evalclusters( ...
            score(:, PCA_best), ...
            'kmeans', ...
            'CalinskiHarabasz', ...
            'KList', max(PCA_best));
    catch
        disp('Error evaluating number of clusters for k-means')
    end

    % K-means clustering in PCA space
    [label, C] = kmeans( ...
        score(:, PCA_best), ...
        eva.OptimalK, ...
        'Replicates', 100, ...
        'Options', opts);

    % Store PCA cluster labels per neuron
    for i = 1:length(CellIDs)
        Cells(CellIDs(i)).PCA_Labels = label(i);
    end
end
%% =========================================================================
%  SORT NEURONS BY CLUSTER ASSIGNMENT
% =========================================================================

Sort_index     = [];
cluster_border = [0];
Cluster_kmeans_Label = [];

label=[Cells(CellIDs(:)).PCA_Labels]';

for i = 1:max(label)
    Sort_index = [Sort_index; find(label == i)];
    cluster_border = [cluster_border, ...
        length(find(label == i)) + cluster_border(end)];
end

% Create cluster index vector
for darab_cluster = 1:length(cluster_border)-1
    Cluster_kmeans_Label( ...
        cluster_border(darab_cluster)+1 : ...
        cluster_border(darab_cluster+1)) = darab_cluster;
end

% Store sorted indices and cluster IDs
for darab_sort_index = 1:length(Sort_index)
    Cells(CellIDs(darab_sort_index)).label_kmeans = Sort_index(darab_sort_index);
    Cells(CellIDs(darab_sort_index)).clusterID_kmeans = label(darab_sort_index);
end

cluster_border = cluster_border + 0.5;

%% =========================================================================
%  CLUSTER-SORTED HEATMAP
% =========================================================================
ax2 = subplot(1,2,2);
imagesc(brain_struct(Sort_index,:));
hold on

% Zero-time reference line
xline(find(Cells(CellIDs(1)).Time == 0), 'LineWidth', 2)

% Horizontal cluster boundaries
line( ...
    repmat(1:size(brain_struct,2), length(cluster_border), 1)', ...
    repmat(cluster_border, size(brain_struct,2), 1), ...
    'LineWidth', 2);

% Time window for display
time_stamp = [ ...
    find(round(Cells(1).Time * 1000) == -400), ...
    find(round(Cells(1).Time * 1000) ==  400) ];

xlim(time_stamp)
xticks(time_stamp)

xticklabels(num2cell( ...
    round(Cells(CellIDs(1)).Time(time_stamp), 2) * 1000))

% Y-axis labels (sorted neurons)
yticks(1:length(CellIDs))
yticklabels(char([Cells(CellIDs(Sort_index)).Structre]))

xlabel('Time (ms)')
title('Suppl.Fig.15A','Clustered Z-scored firing rate during SWR')

axis square
set(ax2, 'CLim', [min(brain_struct(:)), max(brain_struct(:))])
set(ax2, 'FontName', 'Arial', 'FontSize', 8)

colormap(ax2, 'jet')
colorbar(ax2, 'southoutside')

fig = gcf;

fig.PaperUnits = 'points';
fig.Renderer='painters'
fig.PaperPosition = [0 0 2000 2000];
fig.PaperSize = [2000 2000];
saveas(fig,['Fig_4_G_Suppl_15_A'],'fig')
saveas(fig,['Fig_4_G_Suppl_15_A'],'svg')
saveas(fig,['Fig_4_G_Suppl_15_A'],'tiff')

%% =========================================================================
%  STATISTICAL ANALYSIS OF CLUSTERED POPULATION ACTIVITY
% =========================================================================
figure()

% Extract k-means cluster labels
Labels_extracted = [Cells.clusterID_kmeans];

% Define colormap for clusters
color_for_stem = colormap("lines");

% Define time windows (ms) relative to ripple center
time_around_stimulus = [ ...
    find(round(Cells(1).Time * 1000) == -25), ...
    find(round(Cells(1).Time * 1000) ==  25) ];

time_around_spont = [ ...
    find(round(Cells(1).Time * 1000) == -250), ...
    find(round(Cells(1).Time * 1000) == -200) ];

% Remove empty Cells entries
empty_not_Cells = find(~cellfun(@isempty, {Cells.spikes_kernel_mean}));
Cells = Cells(empty_not_Cells);



%% =========================================================================
%  LOOP OVER CLUSTERS
% =========================================================================
subplot_counter = 0;
for cluster_ID = unique(Labels_extracted(~isnan(Labels_extracted)))

    % ---------------------------------------------------------
    % Mean ± STD population response per cluster
    % ---------------------------------------------------------
    subplot(length(unique(Labels_extracted)), 2, cluster_ID + subplot_counter)

    for num_State = 1:length(unique([Cells.Hemisphare]))

        State = [Cells.Hemisphare];

        % Select neurons belonging to current cluster
        ID = find([Cells.clusterID_kmeans] == cluster_ID);

        % Mean firing rate across neurons
        patch_data = mean( ...
            reshape( ...
            [Cells(ID).spikes_kernel_mean], ...
            length(Cells(CellIDs(1)).spikes_kernel_mean), ...
            length(ID)), ...
            2)';

        firing_rate_matrix(cluster_ID)={reshape( ...
            [Cells(ID).spikes_kernel_mean], ...
            length(Cells(CellIDs(1)).spikes_kernel_mean), ...
            length(ID))};

        % Standard deviation across neurons
        errBar_sig = std( ...
            reshape( ...
            [Cells(ID).spikes_kernel_mean], ...
            length(Cells(CellIDs(1)).spikes_kernel_mean), ...
            length(ID))');

        hold on

        % Plot mean trace
        plot(edges, patch_data, ...
            'LineWidth', 1, ...
            'Color', color_for_stem(cluster_ID,:))

        % Shaded error region
        p1 = patch( ...
            [edges fliplr(edges)], ...
            [patch_data + errBar_sig, fliplr(patch_data - errBar_sig)], ...
            color_for_stem(cluster_ID,:));

        p1.EdgeColor = color_for_stem(cluster_ID,:);
        p1.FaceColor = color_for_stem(cluster_ID,:);
        p1.FaceAlpha = 0.3;
        p1.LineStyle = '--';

        title(['Cluster: ' num2str(cluster_ID) ...
            ' | Neurons n: ' num2str(length(ID))])

        % -----------------------------------------------------
        % Summary statistics (stimulus vs spontaneous)
        % -----------------------------------------------------
        Mean_Firing(cluster_ID) =mean(firing_rate_matrix{cluster_ID}(time_around_stimulus(1):time_around_stimulus(2),:),"all"); %avg. across time window and neurons / group
           % mean(patch_data(time_around_stimulus(1):time_around_stimulus(2)));

        STD_Firing(cluster_ID) =std(mean(firing_rate_matrix{cluster_ID}(time_around_stimulus(1):time_around_stimulus(2),:)));
            %std(patch_data(time_around_stimulus(1):time_around_stimulus(2)));

        Mean_Firing_spont(cluster_ID) =mean(firing_rate_matrix{cluster_ID}(time_around_spont(1):time_around_spont(2),:),"all");
            %mean(patch_data(time_around_spont(1):time_around_spont(2)));

        STD_Firing_spont(cluster_ID) =std(mean(firing_rate_matrix{cluster_ID}(time_around_spont(1):time_around_spont(2),:))); 
            %std(patch_data(time_around_spont(1):time_around_spont(2)));

        % Two-sample t-test
        [h,p,ci,stats] = ttest2((mean(firing_rate_matrix{cluster_ID}(time_around_stimulus(1):time_around_stimulus(2),:))), ...
        (mean(firing_rate_matrix{cluster_ID}(time_around_spont(1):time_around_spont(2),:))),'Alpha', 0.01);

        STAT.p(cluster_ID) = p;
        STAT.h(cluster_ID) = h;
        STAT.ci(cluster_ID,1:2) = ci;
        STAT.stats(cluster_ID) = {stats};

        % Axis formatting
        axis square
        ylim([-2 2])
        xline(0)
        xline(0.5)

        clear ID
    end

    axis square

    % ---------------------------------------------------------
    % Bar plot: spontaneous vs ripple-period activity
    % ---------------------------------------------------------
    % subplot(length(unique(Labels_extracted)), 2, cluster_ID + subplot_counter + 1)
    %
    % try
    %     b = bar([ ...
    %         Mean_Firing_spont(cluster_ID); ...
    %         Mean_Firing(cluster_ID) ]);
    %
    %     b.BarWidth = 0.5;
    %     b.FaceColor = color_for_stem(cluster_ID,:);
    %
    %     hold on
    %     errorbar( ...
    %         [Mean_Firing_spont(cluster_ID); Mean_Firing(cluster_ID)], ...
    %         [STD_Firing_spont(cluster_ID); STD_Firing(cluster_ID)])
    %
    %     xlim([0.5 2.5])
    %     ylim([-2 2])
    %     title(round(STAT.p(cluster_ID), 4))
    %     axis square
    % catch
    %     axis square
    % end

    % ---------------------------------------------------------
    % Violin plot: spontaneous vs duirng SWR period activity
    % ---------------------------------------------------------

    subplot(length(unique(Labels_extracted)), 2, cluster_ID + subplot_counter + 1)

    try
        hold on

        % Distributions derived from patch_data

        data_spont = (mean(firing_rate_matrix{cluster_ID}(time_around_spont(1):time_around_spont(2),:)));%patch_data(time_around_spont(1):time_around_spont(2));
        data_stim  = (mean(firing_rate_matrix{cluster_ID}(time_around_stimulus(1):time_around_stimulus(2),:)));%patch_data(time_around_stimulus(1):time_around_stimulus(2));

        data_all = {data_spont, data_stim};
        xpos = [1 2];

        for v = 1:2
            data = data_all{v};

            % Kernel density estimate for violin
            [f, xi] = ksdensity(data);
            f = f ./ max(f) * 0.3;   % violin width scaling

            % Violin shape
            patch( ...
                [xpos(v)+f, xpos(v)-fliplr(f)], ...
                [xi, fliplr(xi)], ...
                color_for_stem(cluster_ID,:), ...
                'FaceAlpha', 0.3, ...
                'EdgeColor', 'none');

            % Jittered data points
            scatter( ...
                xpos(v) + 0.05*randn(size(data)), ...
                data, ...
                15, ...
                color_for_stem(cluster_ID,:), ...
                'filled', ...
                'MarkerFaceAlpha', 0.6);
        end

        % Overlay mean values (same variables as before)
        % plot([1 2], ...
        %     [Mean_Firing_spont(cluster_ID), Mean_Firing(cluster_ID)], ...
        %     'k-', 'LineWidth', 1.5)

        hold on
        errorbar( ...
            [Mean_Firing_spont(cluster_ID); Mean_Firing(cluster_ID)], ...
            [STD_Firing_spont(cluster_ID); STD_Firing(cluster_ID)], 'k', ...
            'LineStyle','none', ...
            'LineWidth',1, ...
            'CapSize',10)
        %

        scatter([1 2], ...
            [Mean_Firing_spont(cluster_ID), Mean_Firing(cluster_ID)], ...
            40, 'k', 'filled')

        % Axes formatting
        xlim([0.5 2.5]);
        ylim([-2 2]);
        xticks([1 2]);
        xticklabels({'Spont', 'SWR'})
        title(['p-vaule (sig. if alpha<0.01): ' num2str(round(STAT.p(cluster_ID), 2))])
        axis square

    catch
        axis square
    end



    subplot_counter = subplot_counter + 1;
end

sgtitle({'Suppl. Fig. 15B', 'Avg. firing rate during SWR'})

fig=gcf;

fig.PaperUnits = 'points';
fig.Renderer='painters'
fig.PaperPosition = [0 0 2000 2000];
fig.PaperSize = [2000 2000];
saveas(fig,['Suppl_15_B'],'fig')
saveas(fig,['Suppl_15_B'],'svg')
saveas(fig,['Suppl_15_B'],'tiff')
