function [cell_metrics]=implicitTimeSpike(rawWaveforms,sample_rate)
adj=0;
%time vectors, and limitations
time_vector = (linspace(0,((size(rawWaveforms,2))-1)./sample_rate,(size(rawWaveforms,2).*2))).*1000;
peak_zero_limit=find(round(time_vector,1)==0.6,1,'first');

%declare memory for variables.
PeakOne=zeros(size(rawWaveforms,1),3);
PeakMinusOne=zeros(size(rawWaveforms,1),3);
PeakZero=zeros(size(rawWaveforms,1),3);
PeakMinusToPeakZero=zeros(size(rawWaveforms,1),3);
for_PCA=zeros(size(rawWaveforms,1),floor(0.0002.*(2*sample_rate)));
wave=zeros(1,size(rawWaveforms,2).*2);
PeakMinusToPeakZero=zeros(size(rawWaveforms,1),1);
PeakZero_PeakOne=zeros(size(rawWaveforms,1),1);
norm_waveform=zeros(1,size(rawWaveforms,2).*2);
stop=[];
peakA=zeros(size(rawWaveforms,1),1);
peakB=zeros(size(rawWaveforms,1),1);

%loop -> waveform by waveform calculate peaks, and implicit times
for m = 1:size(rawWaveforms,1)
    try
        %   disp(m)
        eucldist=[];
        secund_minimum_temp=[];
        secund_minimum=[];


        % wave = interp1([1:size(rawWaveforms,2)],zscore(rawWaveforms(m,:)),[1:0.5:size(rawWaveforms,2),size(rawWaveforms,2)],'spline');
        wave = interp1([1:size(rawWaveforms,2)],rawWaveforms(m,:),[1:0.5:size(rawWaveforms,2),size(rawWaveforms,2)],'spline');
        [~,loc]=min(wave);
        if time_vector(loc)>0.5
            [~,location_deley]=min(abs(time_vector-abs(diff([time_vector(loc),0.5]))));
            adj(m,:)=abs(diff([time_vector(loc),0.5]));
            hossza=length( wave(location_deley:end));
            temp_wave=wave(location_deley:end);
            wave=zeros(size(wave,1),size(wave,2));
            wave(1:hossza)= temp_wave;

        end

        PeakZero(m,1:3)=[(time_vector(find(min(wave(1:peak_zero_limit))==wave,1,'first'))),
            wave((find(min(wave(1:peak_zero_limit))==wave,1,'first'))),
            (find(min(wave(1:peak_zero_limit))==wave,1,'first'))];
        %         eucldist=diag(squareform(pdist([time_vector;wave]','euclidean')),+1);
        %         secund_minimum_temp=sort(eucldist(1:PeakZero(m,3)));
        %         secund_minimum=find(eucldist==secund_minimum_temp(2));
        %         if length(secund_minimum)~=1
        %             secund_minimum=secund_minimum(find(secund_minimum>1,1,'first'))
        %         end


        % findpeaks(wave,time_vector,'Annotate','extents','WidthReference','halfheight');

        Ypk=[];
        Xpk=[];
        [Ypk,Xpk,~,~,~,~,~,~,~] = findpeaks_peti(wave,time_vector,'WidthReference','halfheight');
        [~,idx]=min(abs(Xpk-PeakZero(m,1)));
        minVal=Xpk(idx);

        % plot(time_vector,wave );
        % hold on; scatter( PeakZero(m,1), PeakZero(m,2));
        % hold on; scatter( PeakMinusOne(m,1), PeakMinusOne(m,2));

        PeakMinusOne(m,1:3)=[minVal,wave(find(time_vector==minVal,1,'first')), (find(time_vector==minVal,1,'first'))];
        PeakOne(m,1:3)=[(time_vector(find(max(wave(PeakZero(m,3):end))==wave,1,'first'))),
            (wave(find(max(wave(PeakZero(m,3):end))==wave,1,'first'))),
            find(max(wave(PeakZero(m,3):end))==wave,1,'first')];
        %         hold on; scatter(PeakMinusOne(m,1),PeakMinusOne(m,2),'*')
        %      hold on; plot(time_vector,wave)
        %

        %          hold on; plot(time_vector,[0; diag(squareform(eucldist),+1)])
        %         hold on; plot([time_vector(2:end)],diag(eucldist,-1))
        %
        % wave_invert=wave.*-1;
        % [M,I]=min(wave_invert(1:(find(time_vector.*1000==PeakZero(m,:)))));
        % PeakMinusOne(m,:)=time_vector(I).*1000;
        PeakMinusToPeakZero(m,:)=diff([PeakMinusOne(m,1), PeakZero(m,1)]);
        %     HalfHigh_times_diff(m,:)=diff(wxPk(1,:));

        %
        %
        % % findpeaks(wave,time_vector,'Annotate','extents','WidthReference','halfheight');
        %
        %   Ypk=[];
        %   Xpk=[];
        % [Ypk,Xpk,~,~,~,~,~,~,~] = findpeaks_peti(wave,time_vector.*1000,'WidthReference','halfheight');
        peakA(m,:)=abs(diff([wave(PeakMinusOne(m,3)),wave(PeakZero(m,3))]));
        peakB(m,:)=abs(diff([wave(PeakZero(m,3)),wave(PeakOne(m,3))]));
        % [MB,IB]=max(Ypk);
        % PeakOne(m,:)=(Xpk(IB));
        PeakZero_PeakOne(m,:)=diff([PeakZero(m,1),PeakOne(m,1)]);
        %
        norm_waveform=wave./abs(min(wave(1:peak_zero_limit)));
        stop=PeakZero(m,3)+floor(0.0002.*(2*sample_rate));
        for_PCA(m,:)=diff(norm_waveform(PeakZero(m,3):stop));
        %
Wave(m,:)=wave;
    catch
        disp('error!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')%        figure(2);  hold on; plot(time_vector,wave,'-'); title(['Error in spikewavformID: ' num2str(m)])
    end
end
%clear Ypk Xpk Wpk Ppk iPk bPk bxPk byPk wxPk ;

%cleaning up the results: waveform (P-1) satrt with 0 -> out of analysis,
%extreme values out of analysis
mostOccuringVal = mode(PeakZero(:,1));
Peak_zero_STD=std(PeakZero(:,1));
loct_correct_peak_zero=find(mostOccuringVal-0.05< PeakZero & PeakZero < mostOccuringVal+0.05);
temop_loc_minus_one=(find(0<PeakMinusOne(:,1) & PeakMinusOne(:,1)< mostOccuringVal-Peak_zero_STD)); %it was -
mostOccuringPeakMinusOne = mode(PeakMinusOne(temop_loc_minus_one,1));
loct_correct_peak_minusone=find(mostOccuringPeakMinusOne-0.001< PeakMinusOne(:,1) & PeakMinusOne(:,1) < mostOccuringVal-Peak_zero_STD); %mostOccuringPeakMinusOne+0.1)

[sharedvals,idx] = intersect(loct_correct_peak_zero,loct_correct_peak_minusone,'stable');

mostOccuringPeakPlusOne = mode(PeakOne(sharedvals));
loct_correct_peak_one=find(mostOccuringPeakPlusOne-0.05< PeakOne(:,1) & PeakOne(:,1) < mostOccuringPeakPlusOne+0.1);
[sharedvals_final,idx] = intersect(loct_correct_peak_one,sharedvals,'stable');


cell_metrics.PeaktoTrough = PeakMinusToPeakZero(sharedvals_final);
cell_metrics.TroughtoPeak = PeakZero_PeakOne(sharedvals_final);
cell_metrics.AB_ratio = ((peakB(sharedvals_final)-peakA(sharedvals_final))./(peakA(sharedvals_final)+peakB(sharedvals_final)))';
%cell_metrics.HalfHigh = HalfHigh_times_diff(sharedvals_final);
cell_metrics.Peaks = [PeakMinusOne(sharedvals_final,1) PeakZero(sharedvals_final,1) PeakOne(sharedvals_final,1)];
cell_metrics.Amp =[PeakMinusOne(sharedvals_final,2) PeakZero(sharedvals_final,2) PeakOne(sharedvals_final,2)];
cell_metrics.for_PCA = for_PCA(sharedvals_final,:);
cell_metrics.filtWaveform = rawWaveforms(sharedvals_final,:);
cell_metrics.Time=time_vector;

%x = linspace(0,((size(rawWaveforms,2))-1)./sample_rate,(size(rawWaveforms,2).*2))

x_old = ([0:1:size(rawWaveforms,2)-1]./(sample_rate));
cell_metrics.Time_original=x_old;

% uncomment below
% figure(2)
% subplot(3,3,[1 2 4 5 7 8])

% findpeaks(mean(rawWaveforms(:,:),1),x_old.*1000,'Annotate','extents','WidthReference','halfheight')
% plot( x_old.*1000,mean(rawWaveforms(sharedvals_final,:)).*6,'LineWidth',1)

% uncomment below
% plot(time_vector,mean(Wave,1),'LineWidth',1)
% hold on; scatter(cell_metrics.Peaks(:,1),cell_metrics.Amp(:,1),'*');
% hold on; scatter(cell_metrics.Peaks(:,2),cell_metrics.Amp(:,2),'*');
% hold on; scatter(cell_metrics.Peaks(:,3),cell_metrics.Amp(:,3),'*');

% vline(PeakMinusOne(sharedvals_final))
% vline(PeakZero(sharedvals_final),'k')
% vline(PeakOne(sharedvals_final))
% xlim([0 1])

% uncomment below
% subplot(3,3,[3])
% histogram(PeakMinusOne(sharedvals_final),0:0.01:1)
% %vline(mode(PeakMinusOne(sharedvals_final)));
% xline(mode(PeakMinusOne(sharedvals_final)),'LineWidth',2);
% title('P-1')
% xlabel('Time[msec]')
% subplot(3,3,[6])
% histogram(PeakZero(sharedvals_final),0:0.01:1)
% %vline(mode(PeakZero(sharedvals_final)));
% xline(mode(PeakZero(sharedvals_final)),'LineWidth',2);
% title('P0')
% xlabel('Time[msec]')
% subplot(3,3,[9])
% histogram(PeakOne(sharedvals_final),0:0.01:1)
% %vline(mode(PeakOne(sharedvals_final)));
% xline(mode(PeakOne(sharedvals_final)),'LineWidth',2);
% title('P1')
% xlabel('Time[msec]')

clear Wave adj eucldist x_old wave time_vector temop_loc_minus_one stop sharedvals_final sharedvals secund_minimum_temp secund_minimum rawWaveforms PeakZero_PeakOne PeakZero PeakOne PeakMinusToPeakZero PeakMinusOne peakB peakA peak_zero_limit norm_waveform mostOccuringVal mostOccuringPeakPlusOne mostOccuringPeakMinusOne m loct_correct_peak_zero loct_correct_peak_one loct_correct_peak_minusone idx for_PCA



