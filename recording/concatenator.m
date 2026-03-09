%% Reading .csv file
clc
clear all
close all

path="E:\rEO_10\";
excel_path=strcat(path,"rEO_10_concat_rec_times.csv");

timepoints = readtable(excel_path);
sessions = unique(timepoints.SessionName);

% Order sessions correctly 
session_numbers = cellfun(@(x) str2double(regexp(x, '^\d+', 'match', 'once')), sessions);
[~, sortIdx] = sort(session_numbers);
sessions = sessions(sortIdx);

%% Writing the concatenated data into .dat file
path_concat_data=strcat(path, '128ch_concatenated_sessions_NEW\128ch_concat_data.dat');

[fid, msg]=fopen(path_concat_data, 'W');
assert(fid>0, msg);

%% Iterating over sessions
totalSampleCount=0;
selected_channels=1:1:128;
for i=1:length(sessions)
    concat_data=[];
    session=sessions{i};
    session_path=strcat(path, session,"\amplifier.dat");
   
    idx=find(contains(timepoints.SessionName, session));
    
    % Reading the data
    a=memmapfile(session_path, 'Format','int16');
    RawData=a.Data;
    RawData=reshape(RawData,[128,length(a.Data)./128]);
    
    % Iterating over predefined timewindows
    for j=1:length(idx)
        time_onset=timepoints.Onset_secs_(idx(j));
        time_end=timepoints.End_secs_(idx(j));
        selected_time=round([time_onset*20000,time_end*20000-1]);
        if selected_time(2)>0
            data_windowed=[];
            if selected_time(1)==0
                selected_time(1)=1;
            end
            data_windowed=double(RawData(selected_channels,selected_time(1):selected_time(2)));

            concat_data=cat(2,concat_data,data_windowed);
        end
    end

    % Reshaping and converting to Int16
    reshaped_concat_data = reshape(concat_data, [1, length(selected_channels)*length(concat_data)]);
    int_concat_data=int16(reshaped_concat_data);

    % Appending to the file
    fwrite(fid, int_concat_data, 'int16');

    % Keeping the total sample count
    totalSampleCount=totalSampleCount+length(concat_data);
    % clear a % keep commented, because it causes Matlab to crash
    clear RawData
end

fclose(fid);

%% Creating time.dat file
path_time_data=strcat(path, '128ch_concatenated_sessions_NEW\time.dat');
[fid, msg]=fopen(path_time_data, 'W');
assert(fid>0, msg);

samples=0:1:(length(totalSampleCount)-1);
fwrite(fid, samples, 'int32');

fclose(fid);


