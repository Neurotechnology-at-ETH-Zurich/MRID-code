%% this code for plotting Figure 4C for Magnetic resonance identification tags for ultra-flexible electrodes
% -------------------------------------------------------------------------
% STRUCTURE
%   1) Load cached data needed for Fig 4B/4C and Theta panels
%   2) Figure 4C: Left SWR row-data traces (stacked)
%   3) Figure 4C: Left CWT magnitude heatmap
%   4) Figure 4C: In paper it is superimposed (ripple + SW cut)
%   5) Figure 4B: MRI overlay/topography
%   6) Figure 4C: THETA panels: quiver field + traces + phase/amplitude summary
%   7) Helper: plotWithCI()
% -------------------------------------------------------------------------

%% Figure 4C first column: Plot row data SWR segment for left hemisphere

%% 1) Load caches and open figure canvas
f = figure('Visible','on','Units','normalized','Position',[0 0 1 1]);

load('Figure_4BC.mat');

% load('05_241022_134111_Fig4C_cache.mat');

%% 2) Prepare left-side raw-data variables (from cache)
cacheLR.left.yticklabels_structure = StructureName(logical(skipped_left_channels)); %Raw data SWR, original sampling rate

t    = cacheLR.left.time;
X    = cacheLR.left.traces;
x0   = cacheLR.left.firstSample;
dist = cacheLR.left.distance;
g    = cacheLR.left.gain;

%% 3) Figure 4C (Column 1): Left raw-data traces (stacked by depth)
figRow = subplot(1,4,2);
hold on;

% Stack signals by subtracting first sample and adding depth offsets
plot(figRow, t, X - repmat(x0,1,size(X,2)) + (-dist .* g), 'b');

yticks(figRow, cacheLR.left.yticks);
yticklabels(figRow, cacheLR.left.yticklabels);

figRow.YLim = [(-max(cacheFig4C.distance_left)*g)-100, 0+100];

xlabel('Time (ms)');
title('Figure 4C 1st');

%% Figure 4C second column: Visualize the wavelet transformation magnitude using a heatmap across channels.

%% 4) Unpack CWT cache variables for left/right plot only LEFT
freq = cacheFig4C.frequency;
dL   = cacheFig4C.distance_left;
dR   = cacheFig4C.distance_right;

CWTL   = cacheFig4C.Left.CWT;
CWTLsd = cacheFig4C.Left.CWT_std;

CWTR   = cacheFig4C.Right.CWT;
CWTRsd = cacheFig4C.Right.CWT_std;

leftRef  = cacheFig4C.Left.left_pyr_ch(1,2);
rightRef = cacheFig4C.Right.right_pyr_ch(1,2);

%% 5) Figure 4C (Column 2): Left heatmap (CWT magnitude)
fig1 = subplot(1,4,3);
fig1.FontName = 'Arial';
fig1.FontSize = 12;

surf(dL', freq, CWTL, ...
    'FaceColor','interp', 'EdgeColor','none', 'FaceLighting','gouraud');

clim([0 0.8]);
view(90,90);
title('Figure 4C 2nd');
colormap('jet');

xlim([min(dL), max(dL)]);
ylim([min(freq), max(freq)]);

xlabel('Raw Channel ID - DV-axis (µm) from PyL');
ylabel('Frequency [Hz]');

xticks(fig1, dL);
xticklabels(fig1, flip(cacheLR.left.yticklabels));

c = colorbar(fig1,'eastoutside');
c.Label.String = 'Norm. Mean |WT|';

% Optional MRI overlay (kept commented as original)
% if isfield(cacheFig4C,'hasMRI') && cacheFig4C.hasMRI
%     hAx = axes('Position', fig1.Position, 'XAxisLocation','top', ...
%         'YAxisLocation','right', 'color','none');
%     hold(hAx,'on');
%     plot(hAx, cacheFig4C.Left.MRI.x, cacheFig4C.Left.MRI.y, '--w', 'LineWidth', 2);
%     view(hAx, 90, 90);
%     xlim(hAx, [min(dL), max(dL)]);
% end

%% 6) Figure 4C (Column 2 superimposed): Left line plots (ripple + SW) with CI
fig2 = subplot(1,4,4);
fig2.FontName = 'Arial';
fig2.FontSize = 12;
grid on;

% Ripple line at max ripple freq (normalized)
ixRipL = (freq == double(cacheFig4C.Left.max_ripple_frequency));
mean_y = CWTL(ixRipL, :);
mean_y = (mean_y - min(mean_y)) / (max(mean_y) - min(mean_y));
sem_y  = CWTLsd(ixRipL, :);
plotWithCI(fig2, dL, mean_y, sem_y);

% SW line (original used fixed index 37)
mean_y = CWTL(37, :);
sem_y  = CWTLsd(37, :);
plotWithCI(fig2, dL, mean_y, sem_y);

view(fig2, 90, 90);
xlim(fig2, [min(dL), max(dL)]);
xticks(fig2, dL);

labelsL = arrayfun(@(a,b) sprintf('%d %d', a, b), ...
    cacheFig4C.Left.channels_Neuroscope, ...
    dL - leftRef, 'UniformOutput', false);
xticklabels(fig2, labelsL);

legend1 = legend('Avg. R-Power @ 180Hz','±SD','Avg. SW-Power @ 20Hz','±SD');
set(legend1, 'Position', [0.792298788888441 0.839241166161896 0.0953124979510903 0.0701663182455884]);

title('Figure 4C 2nd (superimposed)');
ylabel('Amplitude (a.u)');

%% 7) Figure 4B MRI topography (left)
if isfield(cacheFig4C,'hasMRI') && cacheFig4C.hasMRI
    hAx = subplot(1,4,1);

    plot(hAx, cacheFig4C.Left.MRI.x, cacheFig4C.Left.MRI.y, '--k', 'LineWidth', 2);
    view(hAx, 90, 90);

    xlim(hAx, [min(dL), max(dL)]);
    ylabel('1D MRI Topography');

    xticks(hAx, dL);
    xticklabels(hAx, flip(cacheLR.left.yticklabels_structure));

    title('Figure 4B');
end

%% 8) Figure 4C Theta Phase (Column 3)
f2 = figure('Visible','on','Units','normalized','Position',[0 0 1 1]);

%% 8.1) Theta quiver: phase vector field (left hemisphere)
fig_quiv = subplot(1,3,1);

x = ThetaPhasePlotCache.left.quiver.x;
y = ThetaPhasePlotCache.left.quiver.y;
u = ThetaPhasePlotCache.left.quiver.u;
v = ThetaPhasePlotCache.left.quiver.v;

quiver(x, y, u, v, 0, '-r','LineWidth',1);

yticks(fig_quiv, ThetaPhasePlotCache.left.quiver.yticks);
yticklabels(fig_quiv, ThetaPhasePlotCache.left.quiver.yticks_labels);
ylim(fig_quiv, ThetaPhasePlotCache.left.quiver.ylim);

xticks(fig_quiv, ThetaPhasePlotCache.left.quiver.xticks);
xlim(fig_quiv, ThetaPhasePlotCache.left.quiver.xlim);
xticklabels(fig_quiv, ThetaPhasePlotCache.left.quiver.xticklabels);

view(90,90);
xlabel('DV-axis (µm) 0µm -> PyL dHP');
ylabel('Time Period of Theta (ms)');
title(['Figure 4C 3rd - Left hemi. Phase Vector Field']);
grid on;

%% 8.2) Figure 4C (Column 3) TOP: Theta traces: PyL / MaxTheta / HFTheta
subplot(1,3,2);

data(1,:) = ThetaPhasePlotCache.left.PyLTheta.data(1,:); % filtered theta from PyL dHP
data(2,:) = ThetaPhasePlotCache.left.MaxTheta.data(1,:); % filtered theta at maxTheta amplitude
data(3,:) = ThetaPhasePlotCache.left.HFTheta.data(1,:); % filtered theta around HF dHP

extract_first_sample = 0; 
plot(ThetaPhasePlotCache.left.Theta.time, data(1,:), 'r');
hold on;

plot(ThetaPhasePlotCache.left.Theta.time, data(2,:), 'k');
hold on;

plot(ThetaPhasePlotCache.left.Theta.time, data(3,:), 'b');

legend1 = legend('PyL Θ','max Θ','HF Θ');
title('Figure 4C 3rd TOP','filtered Θ')
xlim(ThetaPhasePlotCache.left.Theta.xlim);
ylim(ThetaPhasePlotCache.left.Theta.ylim);

yticks(ThetaPhasePlotCache.left.Theta.yticks);
yticklabels(ThetaPhasePlotCache.left.Theta.yticklabels);

xticks(ThetaPhasePlotCache.left.Theta.xticks);
xticklabels(ThetaPhasePlotCache.left.Theta.xticklabels);

%% 8.3)  Figure 4C (Column 4) Avg. Phase + amplitude summary
fig_phase = subplot(1,3,3);
fig_phase.FontName = 'Arial';
fig_phase.FontSize = 12;

% Phase summary (mean +/- std)
fig_info = plot(ThetaPhasePlotCache.left.PhaseAmplitude.xticks, ...
    mean(ThetaPhasePlotCache.left.PhaseAmplitude.ThetaPhase,1), ...
    'DisplayName', 'Phase Data');
hold on;

xticks(ThetaPhasePlotCache.left.PhaseAmplitude.xticks);
xlim([min(ThetaPhasePlotCache.left.PhaseAmplitude.xticks)-1, ...
    max(ThetaPhasePlotCache.left.PhaseAmplitude.xticks)+1]);
xticklabels(ThetaPhasePlotCache.left.PhaseAmplitude.xticklabels);

ylabel('Phase from Pyramidal Layer (degree)');

upper_bound = mean(ThetaPhasePlotCache.left.PhaseAmplitude.ThetaPhase, 1) + ...
    std(ThetaPhasePlotCache.left.PhaseAmplitude.ThetaPhase);
lower_bound = mean(ThetaPhasePlotCache.left.PhaseAmplitude.ThetaPhase, 1) - ...
    std(ThetaPhasePlotCache.left.PhaseAmplitude.ThetaPhase);

v_upper = [ThetaPhasePlotCache.left.PhaseAmplitude.xticks, upper_bound'];
v_lower = [flip(ThetaPhasePlotCache.left.PhaseAmplitude.xticks), flip(lower_bound')];
f = 1:length([v_upper; v_lower]);

patch('Faces', f, 'Vertices', [v_upper; v_lower], ...
    'FaceColor', fig_info.Color, 'EdgeColor', fig_info.Color, ...
    'FaceAlpha', 0.1, 'LineStyle', '--', ...
    'DisplayName', 'Standard Deviation (Data1)');

grid on;
view(90,90);
hold on;

% Overlay axis for normalized amplitude (top/right axis)
hAx = axes('Position', fig_phase.Position, ...
    'XAxisLocation', 'top', 'YAxisLocation', 'right', 'color', 'none');

hAx.YColor = [0.850980401039124 0.325490206480026 0.0980392172932625];
hold on;

plot(hAx, ThetaPhasePlotCache.left.PhaseAmplitude.xticks, ...
    mean(ThetaPhasePlotCache.left.PhaseAmplitude.scale_factor, 1), ...
    'DisplayName', 'Norm.Ampl.', 'Color', hAx.YColor);

plot(hAx, ThetaPhasePlotCache.left.PhaseAmplitude.xticks(ThetaPhasePlotCache.left.PhaseAmplitude.maxPeakIndex), ...
    ThetaPhasePlotCache.left.PhaseAmplitude.avg_norm_amplitude(ThetaPhasePlotCache.left.PhaseAmplitude.maxPeakIndex), ...
    'ro', 'MarkerSize', 10, 'DisplayName', 'Last Peak');

plot(hAx, ThetaPhasePlotCache.left.PhaseAmplitude.xticks(ThetaPhasePlotCache.left.PhaseAmplitude.decline_start_idx), ...
    ThetaPhasePlotCache.left.PhaseAmplitude.avg_norm_amplitude(ThetaPhasePlotCache.left.PhaseAmplitude.decline_start_idx), ...
    'rx', 'MarkerSize', 10, 'DisplayName', 'Decline Start');

ylabel(hAx,'Norm. Ampl.');

% Norm amplitude CI patch (mean +/- std) already prepared in cache
v_upper_scale = [ThetaPhasePlotCache.left.PhaseAmplitude.xticks, ...
    ThetaPhasePlotCache.left.PhaseAmplitude.upper_bound_scale'];
v_lower_scale = [flip(ThetaPhasePlotCache.left.PhaseAmplitude.xticks), ...
    flip(ThetaPhasePlotCache.left.PhaseAmplitude.lower_bound_scale')];

f_scale = 1:length([v_upper_scale; v_lower_scale]);

patch(hAx, 'Faces', f_scale, 'Vertices', [v_upper_scale; v_lower_scale], ...
    'FaceColor', hAx.YColor, 'EdgeColor', hAx.YColor, ...
    'FaceAlpha', 0.1, 'LineStyle', '--', ...
    'DisplayName', 'Standard Deviation');

grid on;
view(90,90);
title(['Figure 4C 4rd - Avg. Phase/Amplitude of Theta']);
xlim(hAx, [min(ThetaPhasePlotCache.left.PhaseAmplitude.xticks), ...
    max(ThetaPhasePlotCache.left.PhaseAmplitude.xticks)]);

%% Save (kept commented as original)
% fig = gcf;
% fig.PaperUnits = 'points';
% fig.PaperPosition = [0 0 1920 1080];
% fig.PaperSize = [1920 1080];
% set(fig, 'Renderer', 'painters');
%
% saveas(fig, [cacheFig4C.session_name '_SWR_CWT_all_STD_fromCache'], 'svg');
% saveas(fig, [cacheFig4C.session_name '_SWR_CWT_all_STD_fromCache'], 'fig');

%% helper
function plotWithCI(ax, x, mean_y, sem_y)
% plotWithCI  Plot mean line with a confidence-interval band (mean +/- sem_y)

upper = mean_y + sem_y;
lower = mean_y - sem_y;

fig_info = plot(ax, x, mean_y', 'LineWidth', 1.5);
hold(ax, 'on');

verts = [x(:), upper(:); flipud(x(:)), flipud(lower(:))];
faces = 1:size(verts,1);

patch(ax, 'Faces', faces, 'Vertices', verts, ...
    'FaceColor', fig_info.Color, 'EdgeColor', fig_info.Color, 'LineStyle','--');

alpha(ax, 0.1);
end
