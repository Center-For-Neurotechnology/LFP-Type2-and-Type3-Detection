%% Detecting Type 2 and 3 events from LFP data
% This is code for the detection of Type 2 and 3 events for the manuscript
% titled "Microscale physiological events on the human cortical surface"
% by Paulk et al., located at https://www.biorxiv.org/content/10.1101/770743v1

% The criteria for detecting Type 2 and 3 events are as follows:
% 1) absolute waveforms were detected which were >25 µV in amplitude with a
% series of steps to clean the waveforms and remove noise
% 2) detected waveforms had a correlation above 0.8 with a template event
% waveform
% 3) the second derivative at the onset of the recording was greater
% than 2
% 4) the voltage in the 100 ms preceding the event onset is less
% than 25 µV to reduce the chances of capturing oscillations

%% Input data:
% LFP data should be in the form of channel x sample,
% sampled at 1000 Hz, with the voltage scale in microVolts
% ChannelNum is the selected channels in the recording to be used in the
% detection (generally channels within good impedance ranges). In this
% example, only a single channel is used.

%% Thresholding to detect peaks in the LFP
% First, we detect peaks larger than 25 ÂµV and take -250 ms and 500 ms
% snippets of data around each event per channel (at ch).
MainDirectory='\YourDirectory\LFP-Type2-and-Type3-Detection-master\';

%Example data from a single channel 
load([MainDirectory,'Examp'],'datalfp','DecFS')

ChannelNum=1;
fs=DecFS; %Sampling rate in Hz.
% datalfp=DataFiltBank1(ChannelNum,1:12*10^5) ;

WaveClassAll=[];
for Deflection=1:2
    if Deflection==1 % Detecting the negatively-deflecting waveforms
        [N ,out1] = spike_times(-(datalfp),25);
    elseif Deflection==2 % Detecting the positively-deflecting waveforms
        [N ,out1] = spike_times((datalfp),25);
    end
    out1(find(diff(out1)<50))=[]; % Remove waveforms which are closer than 50 ms apart
    
    for msn=1:length(out1)
        % A check to capture the full waveform
        if out1(msn)-250>0 && out1(msn)+500<size(datalfp,2)
            
            % A step to capture the time span surrounding the detected peaks in the LFP
            TI=(out1(msn)-250:out1(msn)+500);
            % A step to capture waves surrounding detected peaks in the LFP
            plotD3=datalfp(:,TI);
            % Useful for measuring features of the the waveform snippets
            Rangecheck=250-50:250+50;
            
            %Producing a matrix of the waveform snippets along with
            %measures of the amplitudes and timing of the detected waveforms
            WaveClassAll=[WaveClassAll; Deflection ChannelNum max(plotD3(Rangecheck))...
                min(plotD3(Rangecheck)) max(plotD3(Rangecheck))-min(plotD3(Rangecheck))...
                max(diff(plotD3(Rangecheck))) min(diff(plotD3(Rangecheck))) out1(msn) plotD3];
        end
        msn
    end
end

%% Detecting Type 2 and 3 events from LFP data
% Among those events, we then re-align events relative to the largest
% rising or falling edge (for easier comparisons across waveforms).

%Sort the detected waveforms by time per channel
B = sortrows(WaveClassAll,8);

AllWaveforms=(1:size(B,1))';
recenteredAll=[];
TimingrecenterAll=[];
for RC=1:size(AllWaveforms)
    rng=200:300;
    
    % As the original detection approach can result in repeated waveforms
    % and snippets which overlap in time or are multiples of the same
    % event, this step finds the maximum peak in the snippet recording and
    % re-centers the waveform to -100 to +300 ms
    [mx,id]=max(diff(abs(B(AllWaveforms(RC),rng+9))));
    reCenter=id+rng(1)+9;
    TimingrecenterAll=[TimingrecenterAll;B(AllWaveforms(RC),8)+reCenter-250 reCenter-250];
    recenteredAll=[recenteredAll; B(AllWaveforms(RC),reCenter-100:reCenter+300)];
end

% Removing events arriving at the same time after re-aligning the waveforms
NonrepeatsAll=find(diff(TimingrecenterAll(:,1))>0);

% This matrix will now have the waveform snippets along with timing and
% channel information
Snippets=[TimingrecenterAll(NonrepeatsAll,1:2)...
    repmat([1 ChannelNum],size(NonrepeatsAll,1),1)...
    B(AllWaveforms(NonrepeatsAll),1:8) recenteredAll(NonrepeatsAll,:)];
%%
diff1=1;
diff2=20;
 stepvla=13;
% Capturing the differential value at the image onset
diffStep=diff(Snippets(:,stepvla+100:stepvla+101)');

% Measuring the amplitude of each waveform relative to the time
% period before the detected waveform peak
Amplitude=max((abs(Snippets(:,stepvla+95:stepvla+300)-repmat(nanmean(Snippets(:,stepvla+50:stepvla+95),2),1,length(stepvla+95:stepvla+300))))')';

% This step is to remove large amplitude (>2 milliVolt) and small
% amplitude (<50 microV). Also gets rid of noise where the differential
% waveforms are too large (steps of 20 microVolts in 1 ms).
FI=find(Amplitude<2000 &...
    Amplitude>50 &...
    diffStep'>diff1 &...
    diffStep'<diff2);

WavePositive= Snippets(FI,1:end);
WavePositive(:,stepvla:end)=Snippets(FI,stepvla:end)-repmat(nanmean(Snippets(FI,stepvla+50:stepvla+95),2),1,size(Snippets(FI,stepvla:end),2));

% Some large amplitude events with huge amounts of variation can
% filter through. This step is to remove those waveforms:
FG=find(Amplitude<2000 &...
    Amplitude>50 &...
    diffStep'>-diff2 &...
    diffStep'<-diff1);

WaveNegative= Snippets(FG,1:end);
WaveNegative(:,stepvla:end)=Snippets(FG,stepvla:end)-repmat(nanmean(Snippets(FG,stepvla+50:stepvla+95),2),1,size(Snippets(FG,stepvla:end),2));

WaveAll=sortrows([WavePositive;WaveNegative],1);

% Some large amplitude events with huge amounts of variation can
% filter through. This step is to remove those waveforms:
Variance=var(WaveAll(:,14+10:14+75)');
VarG=WaveAll(Variance<30000,1:end);

% Creating the matrix of data with the cleaned up, aligned waveforms.
WaveAll=[ 1*ones(size(VarG,1),1) VarG];

%% Correlation step 1: Correlating the template waveforms to the waveform snippets
% Loading the waves determined from PCA and kmeans approaches
load([MainDirectory,'\Type2Type3WaveformTemplates'],'Type2WaveA','Type2WaveB','Type3Wave')

% Multiply by 100 microvolts to get the template waveforms to a similar scale
multTimes=[100];

% Normalize the waveforms to the maximal value
Type2WaveA=Type2WaveA/max(Type2WaveA); %Similar waveform as Type2WaveB, merged in the manuscript
Type2WaveB=Type2WaveB/max(Type2WaveB);
Type3Wave=Type3Wave/max(Type3Wave);

% Correlate each waveform snippet with the template waves
HR=[];PR=[];HR2=[];PR2=[];HR3=[];PR3=[];
for smsn=1:size(WaveAll,1)
    [rho,p]=corr(multTimes*Type2WaveA',WaveAll(smsn,14:end)');
    HR(smsn)=rho;
    PR(smsn)=p;
    [rho,p]=corr(multTimes*Type2WaveB',WaveAll(smsn,14:end)');
    HR2(smsn)=rho;
    PR2(smsn)=p;
    [rho,p]=corr(multTimes*Type3Wave',WaveAll(smsn,14:end)');
    HR3(smsn)=rho;
    PR3(smsn)=p;
    smsn
end

%% Correlation step 2: Identifying the highest correlation values between snippets and specific templates
% Since there can be correlation values between the waveform to each
% template, in addition to identifying waveforms > 0.8 correlation values,
% this bit of code is to find the the highest correlation value between 
% each snippet and the slow and fast template waveforms.

% Since the waveform snippets can be varyingly be correlated with
% the Type 2 and Type 3 waveforms, this step is to identify which
% template waveforms had the highest correlation to each snippet

AllHRClassCorr1=[HR' HR2' HR3'];
AllWave=WaveAll(abs(AllHRClassCorr1(:,1))>.8 | abs(AllHRClassCorr1(:,2))>.8 | abs(AllHRClassCorr1(:,3))>.8,:);
AllHRClassCorr=AllHRClassCorr1(abs(AllHRClassCorr1(:,1))>.8 | abs(AllHRClassCorr1(:,2))>.8 | abs(AllHRClassCorr1(:,3))>.8,:);
AllWavePerFind=(1:size(AllWave,1))';

ValMxAll=[repmat([1 1],size(AllHRClassCorr,1),1) AllWavePerFind AllHRClassCorr(:,1) AllHRClassCorr(:,2) AllHRClassCorr(:,3) ];

% Matrix the same size as the waveform snippet with the first column
% indicating which type of waveform the snippet is correlated with for that
% snippet.
Categories=[];
for cs=1:size(ValMxAll,1)
    NF=abs(ValMxAll(cs,4:end));
    NFInd=abs(ValMxAll(cs,3));

    NF(isnan(NF)==1)=0;
    [mx,PerInd]=max(NF);
    DiffPerInd=0;
    if PerInd==3
        DiffPerInd=4;
    elseif PerInd==2
        DiffPerInd=3;
    elseif PerInd==1
        DiffPerInd=2;
    end
    
    vl=[4:6];
    Categories(cs,:)=[PerInd DiffPerInd mx ValMxAll(cs,vl(PerInd))];
end

%% Detecting fast deflections at the onset of the waveforms and removing slow onset changes
% Finally, as a step to remove the possibility of including ongoing
% oscillatory waveforms, we keep only waveforms which:
% 1. Have a large second derivative at ± 50 ms around the
% onset of the waveform
% 2. Keep waveforms with average absolute voltages < 25 µV in the
% preceding 100 ms before the event onset.

% Correlation threshold variable where the correlation with the template
% Type 2 and Type 3 waveforms is above 0.8
thresh=0.8;

% The period of the waveform to look at the second derivative
rangeWaveform=90:110;

WaveOnset=14;

Wave2=[];
Wave3=[];

% Identifies the Type 2 and Type 3 events across the recording per channel
% where the onset
for SC=1:2
    
    if SC==1
        % Selecting the waveforms with correlations to the Type 2 template
        % waveform and including the correlation values for later
        % assessment
        ValBit=ValMxAll(Categories(:,1)<3 &...
            Categories(:,3)>thresh ,3);
        VarG=[AllWave(ValBit,1:end) Categories(Categories(:,1)<3 & Categories(:,3)>thresh,:)];
        
        % This step keeps only Type 2 waves where the first 90 sample points are,
        % on average, less than 25 microVolts and absolute values less than 500 microVolts
        VarG=VarG(abs(VarG(:,WaveOnset+97))<500 & nanmean(abs(VarG(:,WaveOnset+1:WaveOnset+90)),2)<25,:);
        
        % This step keeps only Type 3 waves where the onset of the
        % waveforms have a second differential which is greater than 2 to
        % only keep fast onset waveforms to differentiate from rhythmic
        % waveforms based on the saline tests
        PW=VarG(:,13:end);
        % Second differential calculation
        VM=nanmean(abs(diff(diff((PW(:,rangeWaveform))'))));
        % Thresholding for the higher second differential values
        VarG=VarG(VM>2,:);
    elseif SC==2
        % Selecting the waveforms with correlations to the Type 3 template
        % waveform and including the correlation values for later
        % assessment
        ValBit=ValMxAll(Categories(:,1)==3 &...
            Categories(:,3)>thresh,3);
        VarG=[AllWave(ValBit,1:end) Categories(Categories(:,1)==3 & Categories(:,3)>thresh,:)];
        
        % This step keeps only Type 3 waves where the first 90 sample points are,
        % on average, less than 25 microVolts
        VarG=VarG(nanmean(abs(VarG(:,WaveOnset+1:WaveOnset+90)),2)<25,:);
        
        % This step keeps only Type 3 waves where the onset of the
        % waveforms have a second differential which is greater than 2 to
        % only keep fast onset waveforms to differentiate from rhythmic
        % waveforms based on the saline tests
        PW=VarG(:,13:end);
        % Second differential calculation
        VM=nanmean(abs(diff(diff((PW(:,rangeWaveform))'))));
        % Thresholding for the higher second differential values
        VarG=VarG(VM>2,:);
    end
    
    if SC==1
    Wave2=VarG(:,1:end); %Type 2 waveforms with the correlations to the templates and other amplitude details
    elseif SC==2
    Wave3=VarG(:,1:end); %Type 2 waveforms with the correlations to the templates and other amplitude details
    end        
    
end


%% Plot the results with waveforms and raster plots

if isempty(Wave2)==0
    subplot(2,2,1)
    for WV=1:size(Wave2,1)
        plot((1:size(Wave2(WV,WaveOnset:414),2))/fs,Wave2(WV,WaveOnset:414),'color',[.5 .5 .5 0.5],'linewidth',1)
        hold on
    end
    title('Type 2 detected waveforms')
    xlabel('time (sec)')
    
    subplot(2,2,3)
    stem(Wave2(:,2)/fs,ones(size(Wave2(:,2)/fs)),'color',[.5 .5 .5 0.5],'marker','.')
    xlim([0 size(datalfp,2)/fs])
    xlabel('time (sec)')
    title('raster plot of detected times')
end

if isempty(Wave3)==0
    subplot(2,2,2)
    for WV=1:size(Wave3,1)
        plot((1:size(Wave3(WV,WaveOnset:414),2))/fs,Wave3(WV,WaveOnset:414),'color',[0 0 0 .5],'linewidth',1)
        hold on
    end
    title('Type 3 detected waveforms')
    xlabel('time (sec)')
    
    subplot(2,2,4)
    stem(Wave3(:,2)/fs,ones(size(Wave3(:,2)/fs)),'color','k','marker','.')
    xlim([0 size(datalfp,2)/fs])
    xlabel('time (sec)')
    title('raster plot of detected times, Type 3')
end

