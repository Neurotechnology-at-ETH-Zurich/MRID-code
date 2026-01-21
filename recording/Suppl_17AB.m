%% Supplementary Figure 17A–B (Bundles 1–4)
% This script loads bundle .mat files and plots:
%   - Suppl. 17A: number of neurons over recording days (Left/Right)
%   - Suppl. 17B: percentage of active channels over recording days (Left/Right)


BundleFileName={'Suppl_17AB_Bundle_1_2.mat','Suppl_17AB_Bundle_3_4.mat'}

%% Color definitions (RGB) for bundles 1–4
% Row order corresponds to bundle 1,2,3,4 respectively

subplot(1,2,2)
colors = [
    0.000000000000000  0.447058826684952  0.741176486015320  % bundle 1 Left MRID
    0.466666668653488  0.674509823322296  0.188235297799110  % bundle 2 Right MRID
    0.929411768913269  0.694117665290833  0.125490203499794  % bundle 3 Left MRID
    1.000000000000000  0.000000000000000  0.000000000000000  % bundel 4 Right no-MRID
    ];

% Color index mapping: each row corresponds to a file, columns are [Left, Right]
color_code=[1,2;3,4];

%% Figure setup
bundle_counter=1;
fig1=figure(1);

%% Loop over bundle files
for num_file=1:length(BundleFileName)
    %% ---------- Panel B (Suppl. 17B): Active channels (%) ----------
    subp1=subplot(1,2,2)

    load(char(BundleFileName(num_file)));
    %% Plot partial results
    session_counter=1 % this beacouse we do not have impednace for all and it jump to zero.
    for i=1:length(sortedSessions)
        idx_session = contains(AllData.Session, sortedSessions(i), 'IgnoreCase', true);  % case-insensitive match
        hemi=AllData.Hemisphere(idx_session);
        Active=AllData.Functional(idx_session);
        %split Left Right

        if sum(Active(contains(hemi,'Left')))~=0 && sum(Active(contains(hemi,'Right')))~=0
            Active_per_session(session_counter).left=(sum(Active(contains(hemi,'Left'))));
            Active_per_session(session_counter).right=(sum(Active(contains(hemi,'Right'))));
            session_counter_x_axis(session_counter)=i;
            session_counter=session_counter+1;
        end

    end

    plot1=plot(subp1,session_counter_x_axis,([Active_per_session.left]./FunctionalChannels_Ref_left).*100,'LineWidth',2,'Color',colors(color_code(num_file,1),:))
    hold on
    plot2=plot(subp1,session_counter_x_axis,([Active_per_session.right]./FunctionalChannels_Ref_right).*100,'LineWidth',2,'Color',colors(color_code(num_file,2),:))

    ylim(subp1,[80,104]);
    xlim(subp1,[0.5,16.5]);
    set(plot1,'DisplayName',['bundle (' num2str(bundle_counter) ')']);

    bundle_counter=  bundle_counter+1;
    set(plot2,'DisplayName',['bundle (' num2str(bundle_counter) ')']);

    bundle_counter=  bundle_counter-1;

    % Create ylabel
    ylabel({'Percentage of active channels', '(25 kΩ < Z < 5 MΩ)'});

    % Create xlabel
    xlabel('Recording day (relative to start, day 0)');
    % xticks([2 4 6 8 10 12 14 16])
    % xticklabels([1 26 44 62 79 107 120 140])
    legend(gca,'show');
    xticks([2 4 6 8 10 12 14 16])
    xticklabels([1 26 44 62 79 107 120 140])
    title('Suppl. 17B')
    %% ---------- Panel A (Suppl. 17A): Number of neurons ----------
    subp2=subplot(1,2,1)

    title(subp2,'Suppl. 17A')

    hold on
    plot3=plot(subp2,counter_calid_segment_x_plot,[Neurons_active_per_session.left],'LineWidth',2,'Color',colors(color_code(num_file,1),:))
    hold on
    plot4=plot(subp2,counter_calid_segment_x_plot,[Neurons_active_per_session.right],'LineWidth',2,'Color',colors(color_code(num_file,2),:))

    set(plot3,'DisplayName',['bundle (' num2str(bundle_counter) ')']);
    bundle_counter=  bundle_counter+1;
    set(plot4,'DisplayName',['bundle (' num2str(bundle_counter) ')']);
    bundle_counter=  bundle_counter+1;
    % Create y-axis label
    ylabel('Number of neurons');

    % Create x-axis label
    xlabel('Recording day (relative to start, day 0)');
    legend(gca,'show');
    %    ylim([0 max([[Neurons_active_per_session.left] [Neurons_active_per_session.right]]+5)])
    ylim([0 65])
    xticks([2 4 6 8 10 12 14 16])
    xticklabels([1 26 44 62 79 107 120 140])

    box on

end

%% Save (vector output for publication)
fig = gcf;
fig.PaperUnits = 'points';
fig.PaperPosition = [0 0 1500 500];
fig.PaperSize = [1500 500];
set(fig, 'Renderer', 'painters');

saveas(fig, ['SuppleFig_17AB'], 'svg');
saveas(fig, ['SuppleFig_17AB'], 'fig');



