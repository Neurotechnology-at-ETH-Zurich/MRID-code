%% Magnetic Resonance Identification Tags for Ultra-Flexible Electrodes
% This script generates Figure 4D and Supplementary Figure 12
% from the manuscript:
% "Magnetic Resonance Identification Tags for Ultra-Flexible Electrodes"
%
% Description:
% This script processes electrophysiological and MRI-related data to
% reproduce the figures associated with ultra-flexible electrode
% localization and validation.
%
% If the raw data required for full processing are not available
% (e.g., due to missing files or restricted access), the script
% automatically falls back to loading precomputed variables from a
% cached file. These cached variables are sufficient to reproduce
% the figures without rerunning the full analysis pipeline.
%
% For access to the raw data or questions regarding data availability,
% please contact the author.
%
% Author: Peter Gombkoto, Ph.D.
% Email: pgombkoto@ethz.ch

try

    %% =========================
    % ===== Session Paths ======
    % =========================

    % Animal rEO_10 (Bundles 1 and 2)
    sessions_bundle_1_2 = [
        "G:\Yanik Lab Dropbox\Peter Gombkoto\Localization Manuscript 2024\RAT DATA\rEO_10\recording\02_240925_144727_selected"
        "G:\Yanik Lab Dropbox\Peter Gombkoto\Localization Manuscript 2024\RAT DATA\rEO_10\recording\05_241022_134111_selected"
        "G:\Yanik Lab Dropbox\Peter Gombkoto\Localization Manuscript 2024\RAT DATA\rEO_10\recording\07_241118_134451_selected"
        "G:\Yanik Lab Dropbox\Peter Gombkoto\Localization Manuscript 2024\RAT DATA\rEO_10\recording\10_241213_120506_selected"
        "G:\Yanik Lab Dropbox\Peter Gombkoto\Localization Manuscript 2024\RAT DATA\rEO_10\recording\12_250110_065258_selected"
        "G:\Yanik Lab Dropbox\Peter Gombkoto\Localization Manuscript 2024\RAT DATA\rEO_10\recording\16_250212_151238_selected"
        ];

    % Animal rEO_06 (Bundle 3)
    sessions_bundle_3 = [
        "G:\Yanik Lab Dropbox\Peter Gombkoto\Localization Manuscript 2024\RAT DATA\rEO_06\recording\01_231218_153006"
        "G:\Yanik Lab Dropbox\Peter Gombkoto\Localization Manuscript 2024\RAT DATA\rEO_06\recording\06_240205_150311"
        "G:\Yanik Lab Dropbox\Peter Gombkoto\Localization Manuscript 2024\RAT DATA\rEO_06\recording\18_240513_151304"
        ];

    % Combine all sessions
    sessions = [sessions_bundle_1_2; sessions_bundle_3];

    %% =========================
    % ===== Animal IDs ========
    % =========================

    AnimalID = regexp(sessions, 'rEO_\d+', 'match', 'once');
    AnimalID = string(AnimalID);
    AnimalID_unique = unique(AnimalID);

    clear sessions_bundle_1_2 sessions_bundle_3

    %% =========================
    % ===== Main Loop ==========
    % =========================

    for i = 1:length(sessions)

        currentID = AnimalID(i);

        cd(sessions(i));
        file = dir('*_Theta_Results.mat');
        load(fullfile(file.folder, file.name));

        %% -------- Right Hemisphere --------
        try
            % Absolute difference between Pyr layer and HF layer
            PyrL_Hf_diff_Right(i,:) = abs(diff(CWT_plot.Right_Pyr_HF_MRI_EPHYS'));

            % Detect first session or animal switch
            if i == 1
                fprintf('Selecting electrode start position (0 µm) for %s\n', currentID);
            elseif currentID ~= prevID
                fprintf('Animal changed from %s to %s at row %d — resetting reference\n', ...
                    prevID, currentID, i);
            end

            if i == 1 || currentID ~= prevID
                % Reference positions (Ephys)
                Start_Electrod_position_PyrL_Ephys = CWT_plot.Right_Pyr_HF_MRI_EPHYS(1,2);
                Start_Electrod_position_HF_Ephys   = CWT_plot.Right_Pyr_HF_MRI_EPHYS(2,2);

                % Reference positions (MRI)
                Start_Electrod_position_PyrL_MRI = CWT_plot.Right_Pyr_HF_MRI_EPHYS(1,1);
                Start_Electrod_position_HF_MRI   = CWT_plot.Right_Pyr_HF_MRI_EPHYS(2,1);
            end

            % Absolute displacement relative to first session
            PyrL_Hf_diff_Right_abs_EPHY(:,i) = ...
                CWT_plot.Right_Pyr_HF_MRI_EPHYS(:,2) - ...
                [Start_Electrod_position_PyrL_Ephys; Start_Electrod_position_HF_Ephys];

            PyrL_Hf_diff_Right_abs_MRI(:,i) = ...
                CWT_plot.Right_Pyr_HF_MRI_EPHYS(:,1) - ...
                [Start_Electrod_position_PyrL_MRI; Start_Electrod_position_HF_MRI];

        catch
            fprintf('\rRight hemisphere data not available in %s', char(AnimalID(i)));
        end

        %% -------- Left Hemisphere --------
        try
            PyrL_Hf_diff_Left(i,:) = abs(diff(CWT_plot.Left_Pyr_HF_MRI_EPHYS'));

            if i == 1 || currentID ~= prevID
                % Reference positions (Ephys)
                Start_Electrod_position_PyrL_Ephys_left = CWT_plot.Left_Pyr_HF_MRI_EPHYS(1,2);
                Start_Electrod_position_HF_Ephys_left   = CWT_plot.Left_Pyr_HF_MRI_EPHYS(2,2);

                % Reference positions (MRI)
                Start_Electrod_position_PyrL_MRI_left = CWT_plot.Left_Pyr_HF_MRI_EPHYS(1,1);
                Start_Electrod_position_HF_MRI_left   = CWT_plot.Left_Pyr_HF_MRI_EPHYS(2,1);
            end

            PyrL_Hf_diff_Left_abs_EPHY(:,i) = ...
                CWT_plot.Left_Pyr_HF_MRI_EPHYS(:,2) - ...
                [Start_Electrod_position_PyrL_Ephys_left; Start_Electrod_position_HF_Ephys_left];

            PyrL_Hf_diff_Left_abs_MRI(:,i) = ...
                CWT_plot.Left_Pyr_HF_MRI_EPHYS(:,1) - ...
                [Start_Electrod_position_PyrL_MRI_left; Start_Electrod_position_HF_MRI_left];

        catch
            fprintf('\rLeft hemisphere data not available in %s', char(AnimalID(i)));
        end

        % Store previous animal ID
        prevID = currentID;
    end

    clear prevID currentID CWT_plot CWT_Theta file i

catch
    % Load precomputed variables if processing fails
    load('Fig4D_Suppl_12.mat')
    disp('Loaded precomputed variables for plotting.')
end

%% Figure 4.D.
figure()


data = [mean([PyrL_Hf_diff_Right(:,1);PyrL_Hf_diff_Left(:,1)]) mean([PyrL_Hf_diff_Right(:,2);PyrL_Hf_diff_Left(:,2)])];
errhigh = [std([PyrL_Hf_diff_Right(:,1);PyrL_Hf_diff_Left(:,1)]) std([PyrL_Hf_diff_Right(:,2);PyrL_Hf_diff_Left(:,2)])];
errlow  =[std([PyrL_Hf_diff_Right(:,1);PyrL_Hf_diff_Left(:,1)]) std([PyrL_Hf_diff_Right(:,2);PyrL_Hf_diff_Left(:,2)])];

% Combine data for two conditions
data1 = [PyrL_Hf_diff_Right(:,1); PyrL_Hf_diff_Left(:,1)];
data2 = [PyrL_Hf_diff_Right(:,2); PyrL_Hf_diff_Left(:,2)];
% Prepare data for violin plot
data_all = [data1, data2];

v = violinplot(data_all);

v(1).FaceColor = [1 0 0];
v(2).FaceColor = [0 0 1];

hold on;

% --- SETTINGS ---
jitterAmount = 0.4;
N = length(data1);

% Generate jittered x positions
x1 = 1 + (rand(N,1) - 0.5) * jitterAmount;
x2 = 2 + (rand(N,1) - 0.5) * jitterAmount;

% --- DRAW PAIRED LINES ---
% for i = 1:N
%     plot([x1(i), x2(i)], [data1(i), data2(i)], '-', ...
%         'Color', [0 0 0 0.25], 'LineWidth', 0.8);
% end

% --- SCATTER POINTS ---
scatter(x1, data1, 20, 'red', 'filled', 'MarkerFaceAlpha', 0.5,'MarkerFaceColor','r','MarkerEdgeColor','r','Marker','x');
scatter(x2, data2, 20, 'blue', 'filled', 'MarkerFaceAlpha', 0.5,'MarkerFaceColor','b','MarkerEdgeColor','b','Marker','*');


ylim([-10 400]);

[h,p,ci,stats] = ttest2([PyrL_Hf_diff_Right(:,1)] ,[PyrL_Hf_diff_Right(:,2)]);

% Create ylabel
ylabel({'Localization difference between ','electrophysiology and MRI (µm)'});

% Create title
title(['Fig.4D t-test, p-value:' num2str(round(p,2))])

% Set the remaining axes properties
set(gca,'FontName','Arial','FontSize',12,'XTickLabel',{'PyL','HF'});


% bar(1:2,data)
set(gca,'FontName','Arial','FontSize',12);
hold on

er = errorbar(1:2,data,errlow,errhigh);
er.Color = [0 0 0];
er.LineStyle = 'none';
ylim([-10 400])
fig=gcf;
fig.PaperUnits = 'points';
fig.PaperPosition = [0 0 150 250];
fig.PaperSize = [150 250];
set(fig,'Renderer','painters');
saveas(fig,['Figure_4D_violin'],'svg')
saveas(fig,['Figure_4D_violin'],'jpg')
saveas(fig,['Figure_4D_violin'],'fig')

%% Suppllementary figure 12.B

figure()
subplot(1,2,1)

AnimalID=categorical(AnimalID);
idx =AnimalID == "rEO_10";

plot(PyrL_Hf_diff_Right(idx,1),'r','LineWidth',2,'Marker','o')
hold on
plot(PyrL_Hf_diff_Right(idx,2),'b','LineWidth',2,'Marker','o')
hold on
plot(PyrL_Hf_diff_Left(idx,1),'--r','LineWidth',2,'Marker','o')
hold on
plot(PyrL_Hf_diff_Left(idx,2),'--b','LineWidth',2,'Marker','o')
ylim([-10 500])
xlim([0.9 6.1])
grid on
xlabel('Consecutive months')
ylabel(["Localization difference between"; ...
    "electrophysiology and MRI (µm)"]);

subplot(1,2,2)

idx =AnimalID == "rEO_06";

plot(PyrL_Hf_diff_Right(idx,1),'r','LineWidth',2,'Marker','o')
hold on
plot(PyrL_Hf_diff_Right(idx,2),'b','LineWidth',2,'Marker','o')
ylim([-10 500])
grid on
xlim([0.9 3.1])
xticks([1:0.5:length(PyrL_Hf_diff_Right(idx,1))])
xticklabels([1 2 3 4 5])

xlabel('Consecutive months')
% ylabel(["Localization difference between"; ...
%     "electrophysiology and MRI (µm)"]);
sgtitle('Suppl. 12B')
fig=gcf;
fig.PaperUnits = 'points';
fig.PaperPosition = [0 0 600 200];
fig.PaperSize = [600 100];
set(fig,'Renderer','painters');
saveas(fig,['Suplementary_12_B'],'svg')
saveas(fig,['Suplementary_12_B'],'jpg')
saveas(fig,['Suplementary_12_B'],'fig')


%% Suppllementary figure 12.A


figure('Color','w');

AnimalID = categorical(AnimalID);

% Common plotting parameters
ylims = [-400 400];
lw    = 2;
mk    = 'o';

%% =========================
% ===== EPHYS (Row 1) ======
% =========================

% ---- rEO_10 : Right ----
subplot(2,3,1)
idx = AnimalID == "rEO_10";

plot(PyrL_Hf_diff_Right_abs_EPHY(1,idx), 'r', ...
    'LineWidth', lw, 'Marker', mk); hold on
plot(PyrL_Hf_diff_Right_abs_EPHY(2,idx), '--b', ...
    'LineWidth', lw, 'Marker', mk);

ylabel(["Variation of electrophysiological localization"; ...
    "relative to the first session (µm)"])
legend('PyrL bundle 1','HF bundle 1')

ylim(ylims)
xlim([0.9 6.1])
axis square
grid on


% ---- rEO_10 : Left ----
subplot(2,3,2)

plot(PyrL_Hf_diff_Left_abs_EPHY(1,idx), 'r', ...
    'LineWidth', lw, 'Marker', mk); hold on
plot(PyrL_Hf_diff_Left_abs_EPHY(2,idx), '--b', ...
    'LineWidth', lw, 'Marker', mk);

legend('PyrL bundle 2','HF bundle 2')
xlabel('Consecutive months')

ylim(ylims)
xlim([0.9 6.1])
axis square
grid on


% ---- rEO_06 : Right ----
subplot(2,3,3)
idx = AnimalID == "rEO_06";

plot(PyrL_Hf_diff_Right_abs_EPHY(1,idx), 'r', ...
    'LineWidth', lw, 'Marker', mk); hold on
plot(PyrL_Hf_diff_Right_abs_EPHY(2,idx), '--b', ...
    'LineWidth', lw, 'Marker', mk);

legend('PyrL bundle 3','HF bundle 3')

ylim(ylims)
xlim([0.9 3.1])
axis square
grid on

xticks(1:0.5:length(PyrL_Hf_diff_Right(idx,1)))
xticklabels(1:5)


%% ======================
% ===== MRI (Row 2) =====
% ======================

% ---- rEO_10 : Right ----
subplot(2,3,4)
idx = AnimalID == "rEO_10";

plot(PyrL_Hf_diff_Right_abs_MRI(1,idx), 'r', ...
    'LineWidth', lw, 'Marker', mk); hold on
plot(PyrL_Hf_diff_Right_abs_MRI(2,idx), '--b', ...
    'LineWidth', lw, 'Marker', mk);

ylabel(["Variation of MRI-based loaclization"; ...
    "relative to the first session (µm)"])

legend('PyrL bundle 1','HF bundle 1')

ylim(ylims)
xlim([0.9 6.1])
axis square
grid on


% ---- rEO_10 : Left ----
subplot(2,3,5)

plot(PyrL_Hf_diff_Left_abs_MRI(1,idx), 'r', ...
    'LineWidth', lw, 'Marker', mk); hold on
plot(PyrL_Hf_diff_Left_abs_MRI(2,idx), '--b', ...
    'LineWidth', lw, 'Marker', mk);

legend('PyrL bundle 2','HF bundle 2')
xlabel('Consecutive months')

ylim(ylims)
xlim([0.9 6.1])
axis square
grid on


% ---- rEO_06 : Right ----
subplot(2,3,6)
idx = AnimalID == "rEO_06";

plot(PyrL_Hf_diff_Right_abs_MRI(1,idx), 'r', ...
    'LineWidth', lw, 'Marker', mk); hold on
plot(PyrL_Hf_diff_Right_abs_MRI(2,idx), '--b', ...
    'LineWidth', lw, 'Marker', mk);

legend('PyrL bundle 3','HF bundle 3')

ylim(ylims)
xlim([0.9 3.1])
axis square
grid on

xticks(1:0.5:length(PyrL_Hf_diff_Right(idx,1)))
xticklabels(1:5)
sgtitle('Suppl. 12A')

fig=gcf;
fig.PaperUnits = 'points';
fig.PaperPosition = [0 0 1800 1600];
fig.PaperSize = [1800 1600];
set(fig,'Renderer','painters');
saveas(fig,['Suplementary_12_A'],'svg');
saveas(fig,['Suplementary_12_A'],'jpg');
saveas(fig,['Suplementary_12_A'],'fig');



