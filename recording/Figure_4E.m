%% Magnetic Resonance Identification Tags for Ultra-Flexible Electrodes
%% Figure 4E – Spike amplitude vs. depth and mean spike waveform
%
% This script generates Figure 4E, showing:
%   (1) Average spike amplitude as a function of electrode depth, and
%   (2) Mean spike waveforms for high-amplitude units.
%
% The script supports cached execution: if a precomputed Figure_4E.mat
% file is found in the current directory, the figure is generated
% directly without recomputing intermediate variables.
%
% Author: Peter Gombkoto, Ph.D.


%% FAST PATH: load cached Figure_4E if it exists (current folder)
cacheFile = 'Figure_4E.mat';
d = dir(cacheFile);

if isempty(d)
    fprintf('Figure_4E.mat NOT found in current folder.\n');
    goto_plot = false;  
else
    fprintf('Figure_4E.mat found:\n');
    disp(d)
    load(cacheFile);
    goto_plot = true;    
end

%% =======================================================================
%% COMPUTE ONLY IF CACHE IS MISSING; if you need the files ask: pgombkoto@ethz.ch
if ~goto_plot

    folder = 'G:\Yanik Lab Dropbox\Peter Gombkoto\Localization Manuscript 2024\RAT DATA\rEO_10\recording\05_241022_134111_selected\New_Sorting_for_the_paper\05_241022_134111';

    % --- Load spike sorting outputs
    load(fullfile(folder,'amplifier_res.mat'), 'spikeTimes');
    load(fullfile(folder,'amplifier_res.mat'), 'spikesByCluster');
    load(fullfile(folder,'amplifier_res.mat'), 'spikeSites');
    load(fullfile(folder,'amplifier_res.mat'), 'spikeClusters');
    load(fullfile(folder,'amplifier_res.mat'), 'clusterSites');
    load(fullfile(folder,'amplifier_res.mat'), 'clusterNotes');
    load(fullfile(folder,'amplifier_res.mat'), 'meanWfGlobal');
    load(fullfile(folder,'amplifier_res.mat'), 'spikeAmps');
    load(fullfile(folder,'amplifier_res.mat'), 'unitVpp');
    load(fullfile(folder,'05_241022_134111_Preproc.mat'), 'CWT_plot')

    % --- Load ripple/theta info
    load(fullfile(folder,'curated_ripples.mat'));
    load(fullfile(folder,'05_241022_134111.theta_info.mat'));

    % --- Probe map
    [probePad, shankMap, siteLoc, siteMap, bitScaling, nChans, sampleRate] = ...
        loadProabMap(fullfile(folder,'amplifier.prm')); 

    % --- Raw file length
    filename = 'amplifier.dat';
    FileInfo = dir(fullfile(folder, filename));
    nChansInRawFile = nChans;
    num_samples = FileInfo.bytes / (nChansInRawFile * 2); % int16 = 2 bytes

    % --- Load xml metadata (Neuroscope format) to get channel groups and skips
    [xml, ~] = LoadXml([filename(1:end-3) 'xml']);

    %% --- Initialize LFP struct (small)
    lfp.samplerate   = xml.lfpSampleRate;
    lfp.num_channels = xml.nChannels;
    lfp.allChannels  = [xml.AnatGrps(1).Channels, xml.AnatGrps(2).Channels]' + 1;

    lfp.LFP_structures = categorical( ...
        [repmat(1,1,length(xml.AnatGrps(1).Channels)), repmat(2,1,length(xml.AnatGrps(2).Channels))]', ...
        [1 2], {'dHP left','dHP right'} );

    lfp.skipped_channels = 1 - [xml.AnatGrps(1).Skip, xml.AnatGrps(2).Skip]';
    lfp.siteLoc = siteLoc;

    distance_left_all  = lfp.siteLoc(ismember(lfp.LFP_structures,'dHP left'),  2);
    distance_left      = distance_left_all(logical(lfp.skipped_channels(ismember(lfp.LFP_structures,'dHP left'))));

    distance_right_all = lfp.siteLoc(ismember(lfp.LFP_structures,'dHP right'), 2);
    distance_right     = distance_right_all(logical(lfp.skipped_channels(ismember(lfp.LFP_structures,'dHP right'))));

    CWT_plot.distance_left  = distance_left;
    CWT_plot.distance_right = distance_right;

    %% --- Spike arrays (sorted)
    spikeTime_Real = double(spikeTimes(spikeClusters > 0));
    spikeind_Real  = spikeClusters(spikeClusters > 0);
    spikeTime_Real = spikeTime_Real(:);
    spikeind_Real  = spikeind_Real(:);
    spikeSites     = spikeSites(spikeClusters > 0);

    clear spikeTimes

    [spikeTime_Real, sortindex] = sort(spikeTime_Real);
    spikeind_Real = spikeind_Real(sortindex);
    spikeSites    = spikeSites(sortindex);

    %% --- Putative structure mapping
    PutitativeStructure = {};

    filename = 'rEO_10-CHmap-right_hemisphere-duo.txt';
    [siteMap_Nano, StructureName] = importEminhanTXT(filename);
    siteMap_Nano  = flipud(siteMap_Nano);
    StructureName = flipud(StructureName);

    for n_cluster_Site = 1:length(clusterSites)
        PutitativeStructure(n_cluster_Site) = { ...
            char(StructureName(siteMap_Nano == siteMap(clusterSites(n_cluster_Site)) - 1)) ...
            };
    end

    filename = 'rEO_10-CHmap-left_hemisphere-penta.txt';
    [siteMap_Nano, StructureName] = importEminhanTXT(filename);
    siteMap_Nano = siteMap_Nano + 64;

    for n_cluster_Site = 1:length(clusterSites)
        tmp = StructureName(siteMap_Nano == siteMap(clusterSites(n_cluster_Site)) - 1);
        if ~isempty(tmp)
            PutitativeStructure(n_cluster_Site) = {char(tmp)};
        end
    end

    %% --- Format cell labels
    cellID = unique(spikeind_Real);

    formattedCells = arrayfun(@(x, y) sprintf('channelMap: %d Neuron ID: %d', x, y), ...
        siteMap(clusterSites(cellID))' - 1, cellID, 'UniformOutput', false);

    if length(formattedCells) == length(PutitativeStructure)
        formattedCells = arrayfun(@(idx) sprintf('%s Structure: %s', ...
            formattedCells{idx}, PutitativeStructure{idx}), ...
            1:length(formattedCells), 'UniformOutput', false);
    else
        error('The lengths of formattedCells and PutitativeStructure must match.');
    end

end 

%% =======================================================================
%% PLOT AVG. SPIKE AMPLITUDE
figure();

fig1 = subplot(1,2,2);
temp_color = colormap("copper");
temp_color = temp_color(1:floor(length(temp_color) ./ length(cellID(clusterSites <= 64))):length(temp_color), :);

markerSize = 200;

scatter(distance_left_all(clusterSites(cellID(clusterSites <= 64))), ...
        unitVpp(cellID(clusterSites <= 64)), ...
        markerSize, 'filled', ...
        'MarkerFaceColor', [0.5 0.5 0.5], ...
        'MarkerEdgeColor', [0.5 0.5 0.5]);
hold on;

FigureE_Data_Excel.neuron         = cellID(clusterSites <= 64);
FigureE_Data_Excel.clusterSites   = clusterSites(cellID(clusterSites <= 64))';
FigureE_Data_Excel.unitVpp        = unitVpp(cellID(clusterSites <= 64));
FigureE_Data_Excel.formattedCells = formattedCells(cellID(clusterSites <= 64))';

SelectedUnits = (double(unitVpp) >= 200)';

for colord_neurons = cellID(and(SelectedUnits, clusterSites <= 64))'
    s = char(PutitativeStructure{colord_neurons});

    if contains(s, 'Primary motor')
        temp_color = [0.6 0.9 0.6];
    elseif contains(s, 'Cornu ammonis 1')
        temp_color = [0.27 0.51 0.71];
    elseif contains(s, 'Cornu ammonis 3')
        temp_color = [0.27 0.51 0.71];
    elseif contains(s, 'Dentate gyrus')
        temp_color = [1.0 0.7 0.85];
    elseif contains(s, 'corpus callosum')
        temp_color = [0.7 1.0 0.7];
    end

    scatter(distance_left_all(clusterSites(colord_neurons)), unitVpp(colord_neurons), ...
        markerSize, 'filled', 'MarkerFaceColor', temp_color);

    FigureE_Data_Excel.temp_color(1:3, colord_neurons) = temp_color(:);
end

xticks(fig1,distance_left_all(unique(clusterSites(cellID(clusterSites<=64)))));
[~,wh]=(unique(clusterSites(cellID(clusterSites<=64))));
formattedCells(wh);
xticklabels(fig1,formattedCells(wh));

set(fig1, 'XAxisLocation','top')
box on

title('Avg. Spike Amplitude');
xlim(fig1, [min(distance_left_all) max(distance_left_all)]);
set(fig1, 'XDir', 'reverse');
ylabel('Spike Waveform Amplitude (µm)');
xlabel('Channel number');
view(90,90);

%% Mean waveform
fig2 = subplot(1,2,1);

SelectedUnits = (double(unitVpp) >= 200)';
cellID = unique(spikeind_Real);

SelectedWavforms     = cellID(and(SelectedUnits, clusterSites <= 64));
SelectedClusterSties = clusterSites(SelectedWavforms);

selectedWaveforms = meanWfGlobal(:, SelectedClusterSties, SelectedWavforms);
time_waveform = linspace(-0.5, 1.5, size(selectedWaveforms,1));

timeshift = 0;
for i = 1:max(size(selectedWaveforms,2))
    colord_neurons_for_waveform = cellID(and(SelectedUnits, clusterSites <= 64))';
    s = char(PutitativeStructure{colord_neurons_for_waveform(i)});

    if contains(s, 'Primary motor')
        temp_color = [0.6 0.9 0.6];
    elseif contains(s, 'Cornu ammonis 1')
        temp_color = [0.27 0.51 0.71];
    elseif contains(s, 'Cornu ammonis 3')
        temp_color = [0.27 0.51 0.71];
    elseif contains(s, 'Dentate gyrus')
        temp_color = [1.0 0.7 0.85];
    elseif contains(s, 'corpus callosum')
        temp_color = [0.7 1.0 0.7];
    end

    plot(time_waveform + timeshift, ...
        squeeze(selectedWaveforms(:,i,i))' + distance_left_all(SelectedClusterSties(i)), ...
        'LineWidth', 2, 'Color', temp_color);
    hold on;
  accumulated_timeshift(i)=   timeshift;
    timeshift = timeshift + 2;
  
end

title('Avg. Spike Waveforms (>200 µV)');
yticks(fig2,distance_left_all);
yticklabels(fig2,flipud(distance_left_all - CWT_plot.left_pyr_ch(2)));
ylabel('Dorso/Ventral axis (µm) from PyL dCA1');
xlabel('Selected neuron ID#');

xticks(fig2,accumulated_timeshift+0.5);
xticklabels(fig2, formattedCells(  colord_neurons_for_waveform));

fig1.XLim = fig2.YLim;
sgtitle('Figure 4E');

fig=gcf;
fig.PaperUnits = 'points';
fig.PaperPosition = [0 0 2000 1800];
fig.PaperSize = [2000 1800];
set(fig,'Renderer','painters');
saveas(fig,['Figure_4_E'],'svg');
saveas(fig,['Figure_4_E'],'jpg');
saveas(fig,['Figure_4_E'],'fig');

%% =======================================================================
%% SAVE cache
if ~goto_plot
    cacheFile = fullfile(folder, 'Figure_4E.mat');

    save(cacheFile, ...
        'spikeind_Real', ...
        'distance_left_all', ...
        'clusterSites', ...
        'distance_left', ...
        'distance_left_all', ...
        'unitVpp', ...
        'formattedCells', ...
        'meanWfGlobal', ...
        'CWT_plot', ...
        'cellID', ...
        'PutitativeStructure', ...
        '-v7.3');

    fprintf('Saved cache to:\n  %s\n', cacheFile);
end
