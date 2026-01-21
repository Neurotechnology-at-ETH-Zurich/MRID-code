%% Magnetic Resonance Identification Tags for Ultra-Flexible Electrodes
%% -------------------- Hierarchical Clustering Data --------------------
% Author: Peter Gombkoto, Ph.D.
%
% Description:
%   Loads JRC spike sorting results, assigns putative anatomical structures,
%   computes binned spike-count correlations, and visualizes the resulting
%   correlation matrix with anatomical labels. Also performs hierarchical clustering of the correlation matrix.
%
%   Related figures: Fig. 4F, Supplementary Figs. 13 and 14.
%
% Requirements:
%   - JRC output (_res.mat)
%   - loadProabMap.m
%   - importEminhanTXT.m
%   - CustomColormap.mat

%% -------------------- Check for Pre-ActivityMatrix saved earlier --------------------

checkpointFile = 'rEO_06_preActivityMatrix.mat'; %Suppl.14
%checkpointFile = 'rEO_10_preActivityMatrix.mat'; %Figure 4, Suppl.13

hasCheckpoint = false;

if exist(checkpointFile, 'file') == 2
    load(checkpointFile);   
    hasCheckpoint = true;
    fprintf('Checkpoint loaded: %s\n', checkpointFile);
else
    fprintf('Checkpoint NOT found: %s\n', checkpointFile);
end

if  hasCheckpoint==false

%% -------------------- Parameters --------------------
bin = 25; % Bin size [ms]
load('CustomColormap.mat', 'CustomColormap3')

%% -------------------- Load JRC Spike Sorting Output --------------------
[filename, raw_data_dir] = uigetfile('*.dat');
cd(raw_data_dir);
folder = raw_data_dir;

[filepath, name, ~] = fileparts(filename);

[probePad, shankMap, siteLoc, siteMap, bitScaling, nChans, sample_rate] = ...
    loadProabMap(fullfile(folder, [name '.prm']));

FileInfo   = dir(fullfile(folder, filename));
num_samples = FileInfo.bytes / (nChans * 2); % int16 = 2 bytes

time = linspace( ...
    0, ...
    (num_samples ./ sample_rate), ...
    ceil((num_samples ./ sample_rate) / (bin / 1000)) ...
);


%% -------------------- Extract Session Name --------------------
temp_slash = findstr(raw_data_dir, '\');
try
    session_name = erase(raw_data_dir(temp_slash(end-1):end), '\');
catch
    temp_slash = findstr(raw_data_dir, '/');
    session_name = erase(raw_data_dir(temp_slash(end-1):end), '\');
end


%% -------------------- Extract Animal ID --------------------
% Regular expression to match 'rEO_' followed by one or more digits
pattern = 'rEO_\d+';
match = regexp(raw_data_dir, pattern, 'match');

if ~isempty(match)
    animal = match{1};
    disp(['Animal ID# ', animal]);
else
    disp('Animal ID not found');
end


%% -------------------- Load JRC Result File --------------------
% Create the filename with '_res.mat'
res_filename = fullfile(filepath, [name '_res.mat']);

load(fullfile(folder, char(res_filename)), 'spikeTimes')
load(fullfile(folder, char(res_filename)), 'spikesByCluster')
load(fullfile(folder, char(res_filename)), 'spikeSites')
load(fullfile(folder, char(res_filename)), 'spikeClusters')
load(fullfile(folder, char(res_filename)), 'clusterSites')
load(fullfile(folder, char(res_filename)), 'clusterNotes')


%% -------------------- Preprocess Spike Times --------------------
% IDs less than or equal to 0 are noise.
spikeTime_Real = double(spikeTimes(spikeClusters > 0)); % must be double
spikeind_Real  = spikeClusters(spikeClusters > 0);     % keep as int32

spikeTime_Real = spikeTime_Real(:); % must be column vector
spikeind_Real  = spikeind_Real(:);
spikeSites     = spikeSites(spikeClusters > 0);

clear spikeTimes


%% -------------------- Sort Spikes by Time --------------------
% Spike times must be in ascending order
[spikeTime_Real, sortindex] = sort(spikeTime_Real);

spikeind_Real = spikeind_Real(sortindex);
spikeSites    = spikeSites(sortindex);

CellsID = unique(spikeind_Real);

%% -------------------- Assign Putative Anatomical Structures --------------------
PutitativeStructure = {};

%% -------------------- Animal: rEO_10 --------------------
if contains(animal, 'rEO_10')

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

 
    %% -------------------- Format Cell Labels (Left Hemisphere) --------------------

    cellIDs = unique(spikeind_Real);
    clusterSites_selected = (clusterSites(ismember(PutitativeStructure(:,2),'Left hem.')));

    cellID = cellIDs (ismember(clusterSites,clusterSites_selected ));

    formattedCells = arrayfun( ...
        @(x, y) sprintf('channelMap: %d Neuron ID: %d', x, y), ...
        siteMap(clusterSites_selected(cellID))' - 1, ...
        cellID, ...
        'UniformOutput', false);

    % % Preallocate cell array
    % outputCell = cell(length( clusterSites_selected), 1);
    % 
    % % Fill with formatted strings
    % for i = 1:length( clusterSites_selected)
    %     formattedCells{i} = ...
    %         sprintf('Channel: %d Neuron ID: %d',  clusterSites_selected(i), cellID(i));
    % end
    %% -------------------- Build Activity Matrix --------------------
    clear Activitymatrix

  %  cellID = find(64 >= clusterSites); % select valid clusters
    Activitymatrix = zeros(length(cellID), length(time));

    for Cell_ID = 1:length(cellID)
        Activitymatrix(Cell_ID, :) = ...
            histc((spikeTime_Real(spikeind_Real == Cell_ID)) ./ sample_rate, time);
    end

%% -------------------- Animal: rEO_06 --------------------
elseif contains(animal, 'rEO_06')

    % Left hemisphere only; 
    siteMap_right = siteMap(shankMap == 1); 
   
    [Site, Structure] = importEminhanTXT('channel_atlas_coordinates 2.txt');
    siteMap_Nano  = flip(Site);
    StructureName = flip(Structure);

    PutitativeStructure = {};

    for n_cluster_Site = 1:length(clusterSites(ismember(clusterSites, siteMap_right)))
        PutitativeStructure(n_cluster_Site,1) = ...
            {char(StructureName(siteMap_Nano == siteMap_right(clusterSites(n_cluster_Site)) - 1))};
        PutitativeStructure(n_cluster_Site,2) = {'Left hem.'};
    end


    %% -------------------- Format Cell Labels (Right Hemisphere) --------------------
    clear Activitymatrix

    cellIDs = unique(spikeind_Real);
    clusterSites_right = clusterSites(ismember(clusterSites, siteMap_right));

    cellID = cellIDs(ismember(clusterSites, siteMap_right));

    formattedCells = arrayfun( ...
        @(x, y) sprintf('channelMap: %d Neuron ID: %d', x, y), ...
        siteMap_right(clusterSites_right(cellID))' - 1, ...
        cellID, ...
        'UniformOutput', false); % D->V direction

    % Preallocate cell array
    %outputCell = cell(length(clusterSites_right), 1);

    % % Fill with formatted strings
    % for i = 1:length(clusterSites_right)
    %     formattedCells{i} = ...
    %         sprintf('Channel: %d Neuron ID: %d', clusterSites_right(i), cellID(i));
    % end

    %% -------------------- Build Activity Matrix --------------------
    clear Activitymatrix

    Activitymatrix = zeros(length(cellID), length(time));

    for Cell_ID = 1:length(cellID)
        Activitymatrix(Cell_ID, :) = ...
            histc((spikeTime_Real(spikeind_Real == Cell_ID)) ./ sample_rate, time);
    end

end


%% -------------------- Save Data for plotting --------------------

checkpointFile = [animal '_preActivityMatrix.mat'];

if exist(checkpointFile, 'file')
    fprintf('Loading checkpoint: %s\n', checkpointFile);
    load(checkpointFile, ...
        'cellID',...
        'Activitymatrix', ...
        'spikeTime_Real', ...
        'spikeind_Real', ...
        'clusterSites', ...
        'sample_rate', ...
        'time', ...
        'formattedCells', ...
        'PutitativeStructure', ...
        'CustomColormap3', ...
        'animal');

    % Jump directly to Activity Matrix
else
    fprintf('Saving checkpoint: %s\n', checkpointFile);

    save(checkpointFile, ...
        'cellID',...
        'Activitymatrix', ...
        'spikeTime_Real', ...
        'spikeind_Real', ...
        'clusterSites', ...
        'sample_rate', ...
        'time', ...
        'formattedCells', ...
        'PutitativeStructure', ...
        'animal', ...
        'CustomColormap3', ...
        '-v7.3');
end

end


%% -------------------- Correlation Matrix Computation --------------------
Activitymatrix = flipud(Activitymatrix);
zSpikeCount = zscore(Activitymatrix');
CorrMatrix  = corr(zSpikeCount);


%% -------------------- Correlation Matrix Visualization --------------------
figure()

if strcmp(animal, 'rEO_06')
    CorrMatrixFlipped = flipud(fliplr(CorrMatrix - diag(diag(CorrMatrix))));
    labels_for = flip(PutitativeStructure(64 >= clusterSites));
elseif strcmp(animal, 'rEO_10')
    CorrMatrixFlipped =(CorrMatrix - diag(diag(CorrMatrix)));
    labels_for = flip(PutitativeStructure(ismember(PutitativeStructure(:,2), 'Left hem.')));
end

imagesc(CorrMatrixFlipped)

xlim([0.5 max(cellID) + 0.5])
ylim([0.5 max(cellID) + 0.5])

axis square

if strcmp(animal, 'rEO_06')
   labels_for = (PutitativeStructure(ismember(PutitativeStructure(:,2), 'Left hem.')));

title(['Pairwise Neuronal Spike-Count Correlation ' strrep(animal, '_', '-')],'Suppl.14A')
elseif strcmp(animal, 'rEO_10')
    labels_for = flip(PutitativeStructure(ismember(PutitativeStructure(:,2), 'Left hem.')));

title(['Pairwise Neuronal Spike-Count Correlation ' strrep(animal, '_', '-')],'Fig.4F')
end

colorbar
colormap(CustomColormap3)
clim([-0.4 0.4])


%% -------------------- Axis Labels --------------------
xticks(cellID)
xticklabels(labels_for)
xtickangle(45)

yticks(cellID)
yticklabels(labels_for)


%% -------------------- Save Figure --------------------
fig = gcf;

fig.PaperUnits    = 'points';
fig.PaperPosition = [0 0 3840 2160];
fig.PaperSize     = [3840 2160];

saveas(fig, [animal '_Correlataion_Matrix'], 'svg');
saveas(fig, [animal '_Correlataion_Matrix'], 'fig');


%% -------------------- Hierarchical Clustering --------------------
if strcmp(animal, 'rEO_06')
   labels_for = flip(labels_for);
end

[Labels] = hierachical_clastering_PG(CorrMatrix, labels_for', CustomColormap3);


fig=gcf;

if strcmp(animal, 'rEO_06')
sgtitle(fig,['animalID: ' strrep(animal, '_', '-') ' Suppl.14B-C'])
elseif strcmp(animal, 'rEO_10')
sgtitle(fig,['animalID: ' strrep(animal, '_', '-') ' Suppl.13A-B'])
end

fig.PaperUnits = 'points';
fig.PaperPosition = [0 0 3840 2160];
fig.PaperSize = [3840 2160];
saveas(fig,[animal '_Hierarchichal_Clustering'],'svg');

saveas(fig,[animal '_Hierarchichal_Clustering'],'fig');




% UniqueStructure=unique(labels_for)
% load('HierarchicalClustering.mat')
%
% CA1_percent_cluster=(sum(ismember(CA1.RowNodeNames,UniqueStructure(1)))./length(CA1.RowNodeNames)).*100
% M1_percent_cluster=(sum(ismember(M1.RowNodeNames,UniqueStructure(3)))./length(M1.RowNodeNames)).*100
% DG_percent_cluster=(sum(ismember(DG.RowNodeNames,UniqueStructure(2)))./length(DG.RowNodeNames)).*100
%
%
%
%
%
%
%
% fig=gcf;
%
% fig.PaperUnits = 'points';
% fig.PaperPosition = [0 0 3840 2160];
% fig.PaperSize = [3840 2160];
% saveas(fig,['hierachical_clastering_PG'],'svg');
% saveas(fig,['hierachical_clastering_PG'],'fig');
