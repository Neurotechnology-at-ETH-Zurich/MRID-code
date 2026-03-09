%% Chronic Impedances

mainPath= ".\";


% Surgery date for rEO_06
surgeryDate = datetime(2023, 12, 08);

files=dir(mainPath);
numSessions=0;
for i=1:length(files)
    filename=files(i).name;
    if not(contains(filename, ".")) && not(contains(filename, "128ch")) && not(contains(filename, "ignore"))
        splitstr=strsplit(filename, '_');
        weekID=splitstr{1};
        weekID=str2double(weekID);
        numSessions = numSessions + 1;
    end
end

datesArr = datetime(zeros(numSessions,1), 0, 0);
impedanceSelection = zeros(numSessions,1);
impFiles=[];

sessionIdx = 1;
for i=1:length(files)
    filename=files(i).name;
    if not(contains(filename, ".")) && not(contains(filename, "128ch")) && not(contains(filename, "ignore"))
        splitstr=strsplit(filename, '_');
        dateChr=splitstr{2};
        year = str2double(dateChr(1:2));
        year = 2000 + year;
        month = str2double(dateChr(3:4));
        day = str2double(dateChr(5:6));
        date = datetime(year, month, day);
        datesArr(sessionIdx)=date;
        
        sessionDir = strjoin([mainPath, filename], '');
        sessionFiles = dir(sessionDir);
        
        for k=1:length(sessionFiles)
            filename2=sessionFiles(k).name;
            if contains(filename2, "impedances")
                impedanceSelection(sessionIdx)=1;
                impFiles = [impFiles; strjoin([sessionDir, filename2], '\')];
                
            end
        end
        sessionIdx = sessionIdx + 1;
    end
end
impedanceSelection = logical(impedanceSelection);

numDays = daysact(datesArr, surgeryDate);
numDays = abs(numDays);

allWorkingImp=zeros(128,length(impFiles));
refWorkingChIdx = zeros(128,1);

numWorkingChs=zeros(length(impFiles), 1);
workingChArr = zeros(64, length(impFiles));

meanImps=zeros(length(impFiles), 1);
stdImps=zeros(length(impFiles), 1);

threshold=5e6;
thresholdLower = 5e4;

for k=1:length(impFiles)
   filename=impFiles(k);
   
   if contains(filename, ".csv")
    splitstr=strsplit(filename, '_');
    impCSV=readtable(filename);

    portAidx=string(impCSV.Port) == 'Port A';
    portBidx=string(impCSV.Port) == 'Port B';
    portSelected = logical(portAidx);
   
    allImp=impCSV.ImpedanceMagnitudeAt1000Hz_ohms_(portSelected);
    
    if k==1
        workingIdx1 = allImp>thresholdLower;
        workingIdx2 = allImp<threshold;
        workingChsIdx = logical(workingIdx1 .* workingIdx2);
        refWorkingChIdx = workingChsIdx;
    else
        workingIdx1 = allImp>thresholdLower;
        workingIdx2 = allImp<threshold;
        workingChsIdx = logical(workingIdx1 .* workingIdx2);
    end
    %workingChsIdx = logical(workingChsIdx.*refWorkingChIdx);
    workingChsIdx = logical(workingChsIdx);
    workingChArr(:,k) = allImp/1e3;

    allWorkingImp(workingChsIdx, k)=allImp(workingChsIdx);
    numWorkingChs(k)=sum(workingChsIdx);
    meanImps(k)=mean(allImp(workingChsIdx));
    stdImps(k)=std(allImp(workingChsIdx));
   end
end


meanImps=nonzeros(meanImps);
meanImps=meanImps(sortIdx);
meanImps=meanImps';

stdImps=nonzeros(stdImps);
stdImps=stdImps(sortIdx);
stdImps=stdImps';

numWorkingChs=nonzeros(numWorkingChs);
numWorkingChs=numWorkingChs(sortIdx);
numWorkingChsT = numWorkingChs';
figure(5)

errorbar(numDays(logical(impedanceSelection)), meanImps/1e3, stdImps/1e3./sqrt(numWorkingChsT))
xlabel("Weeks")
ylabel("Impedance")

%% Loading the data
clc
foldername = "128ch_concatenated_sessions";
filename = "128ch_concat_data.dat";
fullpath = fullfile(mainPath, foldername, filename);
a=memmapfile(filename, 'Format','int16');
length(a.Data)./128;
rawData=a.Data;
rawData=reshape(rawData,[128,length(a.Data)./128]);

%% 
clc
JRC_filename = "128ch_concat_data_filt.jrc";
jrc_fullpath = fullfile(mainPath, foldername, JRC_filename);
fid=fopen(jrc_fullpath,'r');
waveform_filt=fread(fid,Inf,'int16');                                                           % spikewavform filtered by JRClust RAW data from binary
fclose(fid);

matname = "128ch_concat_data_res.mat";
matfullpath = fullfile(mainPath, foldername, matname);
load(matfullpath);
spikesFilt=double((squeeze(single( ...
    reshape(waveform_filt,rawShape(1),rawShape(2),rawShape(3)) ...
    ) ...
    ) ...
    ).*(37.4./192)); %uV/bit / gain 

%% Loading the time window onset/offsets where data is concatenated from
clc 
foldername = "128ch_concatenated_sessions";
concat_times_name ="rEO_06_concat_rec_times.csv";
excel_path = fullfile(mainPath, foldername, concat_times_name);
timepoints = readtable(excel_path);
sessions= unique(timepoints.SessionName);

%% Find the session bounds in terms of samplepoint
clc
fs=2e4;
tmp_idx=0;
ses_bounds=[1];


for i=1:length(sessions)
    session=sessions{i};
    splitstr=strsplit(session, '_');
    dateChr=splitstr{2};
    snrDatesArrTmp{i}=dateChr;
end
snrDatesArrTmp = unique(snrDatesArrTmp);
snrDatesArr = datetime(zeros(length(snrDatesArrTmp),1), 0, 0);

for i=1:length(snrDatesArrTmp)
    dateChr=snrDatesArrTmp{i};
    year = str2double(dateChr(1:2));
    year = 2000 + year;
    month = str2double(dateChr(3:4));
    day = str2double(dateChr(5:6));
    date = datetime(year, month, day);
    snrDatesArr(i)=date;

    idx=find(contains(timepoints.SessionName, dateChr));
    for j=1:length(idx)
        time_onset=round(timepoints.Onset_secs_(idx(j))*fs);
        time_end=round(timepoints.End_secs_(idx(j))*fs);
        len=time_end-time_onset;
        tmp_idx=tmp_idx+len;
    end

    ses_bounds=[ses_bounds, tmp_idx];

end

snrNumDays = daysact(snrDatesArr, surgeryDate);
snrNumDays = abs(snrNumDays);

num_sessions = length(ses_bounds)-1;
%% 
clc
siteMap =[20 21 41 18 40 19 43 16 42 23 45 22 44 25 38 24 39 27 36 26 37 29 34 28 35 31 32 17 33 15 30 14 46 13 47 12 49 11 48 10 51 9 50 8 53 7 52 6 55 5 54 4 57 3 56 2 59 1 58 0 61 63 60 62]+1;

selected_channels=siteMap;
lastSelectedCluster=sum(clusterSites <= max(siteMap));
unqClusters=1:1:lastSelectedCluster;

clstChrAmp=zeros([length(ses_bounds)-1, length(unqClusters)]);
chrMeanWf=zeros(length(sessions),length(unqClusters), 41); 
sesVrms=zeros([length(ses_bounds)-1, length(unqClusters)]);
singleUnitsPerSession = zeros([length(ses_bounds)-1, 1]);
spikeIDs=int8(1:1:nSpikes);

ses_wf=[];
ses_std=[];


singleUnits=spikeClusters>0;

for i=2:length(ses_bounds)
    idx_start=ses_bounds(i-1);
    idx_end=ses_bounds(i)-1;

    dataWin=double(rawData(selected_channels,idx_start:idx_end)).*0.195;
    dataWin=bandpass(dataWin.', [300 5000], fs).';

    sesIdx=(spikeTimes>=idx_start & spikeTimes<=idx_end);
    sesSites=spikeSites(spikeTimes>=idx_start & spikeTimes<=idx_end);
    sesClustrs=spikeClusters(spikeTimes>=idx_start & spikeTimes<=idx_end);
    sesIdx=logical(sesIdx.*singleUnits);

    singleUnitsPerSession(i-1) = length(unique(spikeClusters(sesIdx)));

    for j=1:length(unqClusters)
        clstIdx=logical(sesIdx.*(spikeClusters==j));
        wftmp=squeeze(spikesFilt(:,1,clstIdx));
        wftmp=wftmp';
        chrMeanWf(i-1,j,:)=mean(wftmp);
        sesVrms(i-1,j)=rms(dataWin(clusterSites(j), :));
        [cell_metrics]=implicitTimeSpike(mean(wftmp),fs);
        if ~isempty(cell_metrics.Amp)
            % clstChrVpp(i-1,j)=cell_metrics.Amp(:,3)-cell_metrics.Amp(:,2);
            clstChrAmp(i-1,j)=abs(cell_metrics.Amp(:,2));
        end
    end    

    %Uncomment below to run on all spike instances
    neuron=squeeze(spikesFilt(:,1,sesIdx));
    neuron=neuron';
    [cell_metrics]=implicitTimeSpike(neuron,fs);
    meanVpp=mean(cell_metrics.Amp(:,3)-cell_metrics.Amp(:,2));
    stdVpp=std(cell_metrics.Amp(:,3)-cell_metrics.Amp(:,2));
    ses_wf = [ses_wf meanVpp];
    ses_std=[ses_std stdVpp];
    %Uncomment above to run on all spike instances

end

clear dataWin

%%
selectedClusterNotes=clusterNotes(unqClusters);
singleUnitsManual=zeros(length(selectedClusterNotes),1);
singleUnitsManualJRC=zeros(length(selectedClusterNotes),1);
for i=1:length(selectedClusterNotes)
    if ~isempty(selectedClusterNotes{i})
        if selectedClusterNotes{i}=="single"
            singleUnitsManual(i)=1;
        end
    else
        if singleUnits(i)==1
            singleUnitsManualJRC(i)=1;
        end
    end
end
singleUnitsManual=logical(singleUnitsManual);
singleUnitsManualJRC=logical(singleUnitsManualJRC);

%%
% If manual single units labels available replace singleUnits with
% singleUnitsManual below.
clc

snrAmp=clstChrAmp(:,singleUnitsManual)./sesVrms(:,singleUnitsManual);
meanSnr=zeros(num_sessions,1);
stdSnr=zeros(num_sessions,1);
nonzeroUnits=zeros(num_sessions,1);

for i=1:(length(ses_bounds)-1)
    tmpSnr=snrAmp(i,:);
    nonzeroUnits(i)=sum(tmpSnr>0);
    meanSnr(i)=mean(nonzeros(tmpSnr));
    stdSnr(i)=std(nonzeros(tmpSnr));
end

%% Chronic spikes plot
clc
close all

ch_offset=2500;
imp_scaler=5;
ch_scaler=1000;
snr_scale=100;
week_offset=10;
w=3;

[C, ia] = unique(numDays);
daysPlotImpedances = numDays(ia);
meanImpsPlot = meanImps(ia);
stdImpsPlot = stdImps(ia);
numWorkingChsPlot = numWorkingChsT(ia);

zMeanSnr = (meanSnr-mean(meanSnr))/std(meanSnr);
plotSnr=meanSnr(abs(zMeanSnr)<3);
plotSnrStd=stdSnr(abs(zMeanSnr)<3);
plotSnrNumDays = snrNumDays(abs(zMeanSnr)<3);
nonzeroUnitsPlot = nonzeroUnits(abs(zMeanSnr)<3);

selectedSingleUnits=[13, 23, 25, 31, 42];     % Selected units for rEO_06
cmap=colormap("lines");

c_offset = 10;


figure();
f3=gcf;
f3.PaperUnits = 'points';
f3.PaperPosition = [0 0 1000 500];
f3.PaperSize = [1200 1920];
set(f3,'Renderer','painters');

grid on

yyaxis left
errorbar(plotSnrNumDays*week_offset, plotSnr*snr_scale, plotSnrStd*snr_scale./sqrt(nonzeroUnitsPlot), 'LineWidth', 2)
yticks((0:5:max(plotSnr)*2)*snr_scale);
yticklabels(0:5:max(plotSnr)*2)
ylim([0, length(selectedSingleUnits)*ch_scaler + 1.5*ch_offset])
hold on
for i=1:length(plotSnrNumDays)
    for j=1:length(selectedSingleUnits)
        singleUnit=selectedSingleUnits(j);
        plot(plotSnrNumDays(i)*week_offset-w:(2*w/40):plotSnrNumDays(i)*week_offset+w, ...
            ch_offset+squeeze(chrMeanWf(i,singleUnit,:))+j*ch_scaler,'-','Color', cmap(j*c_offset,:), 'LineWidth', 2);
        hold on
    end
end

xlabel("weeks")
ylabel("SNR")

yyaxis right

errorbar(daysPlotImpedances*week_offset, ...
    meanImpsPlot/1e3*imp_scaler, ...
    (stdImpsPlot/1e3*imp_scaler)./sqrt(numWorkingChsPlot), ...
    'LineWidth', 2)


yticks((0:100*imp_scaler:max(meanImpsPlot/1e3*imp_scaler)*2));
yticklabels(0:100:max(meanImpsPlot/1e3)*2)
xlabel("Weeks")
ylabel("Impedance")

xticks(daysPlotImpedances*week_offset)
xticklabels(daysPlotImpedances)
ylim([0, length(selectedSingleUnits)*ch_scaler + 1.5*ch_offset])

path=".\";
name="";

%saveas(f3, fullfile(path, name), 'svg')

%%
clc
filename = "Figure5.xlsx";
for j=1:length(selectedSingleUnits)
    singleUnit=selectedSingleUnits(j);
    singleUnitChrWf = squeeze(chrMeanWf(:,singleUnit,:))';
    singleUnitChrWf(:,end) = [];
    xlswrite(filename, singleUnitChrWf', sprintf("selectedUnit%i",singleUnit) , 'B2');
    xlswrite(filename, cellstr(snrDatesArr), sprintf("selectedUnit%i",singleUnit) , 'A2');
    xlswrite(filename, {"Sessions"}, sprintf("selectedUnit%i",singleUnit) , 'A1');
    spike_sampletimes = round(linspace(-0.5, 1.5, 41),2);
    spike_sampletimes_cellArr = {};
    for k=1:41
        spike_sampletimes_cellArr(k) = {sprintf("Time %.2f ms", spike_sampletimes(k))};
    end
    xlswrite(filename, spike_sampletimes_cellArr, sprintf("selectedUnit%i",singleUnit) , 'B1');
    % xlswrite(filename, samples, sprintf("selectedUnit%i",singleUnit) , 'B1');
end
sheetName="Chronic Impedances";
xlswrite(filename, meanImpsPlot'/1e3, sheetName , 'B2');
xlswrite(filename, {"Mean Impedances (kOhms)"}, sheetName , 'B1');
xlswrite(filename, stdImpsPlot'/1e3, sheetName , 'C2');
xlswrite(filename, {"S.t.d Impedances (kOhms)"}, sheetName , 'C1');
xlswrite(filename, numWorkingChsPlot', sheetName , 'D2');
xlswrite(filename, {"Number of Working Channels per Session (Z <5 MOhm)"}, sheetName , 'D1');
xlswrite(filename, {"Sessions"}, sheetName , 'A1');
xlswrite(filename, cellstr(datesArr(ia)), sheetName , 'A2');
