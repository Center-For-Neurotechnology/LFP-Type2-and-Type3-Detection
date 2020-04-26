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
% detection (generally channels within good impedance ranges)

%% Thresholding to detect peaks in the LFP
% First, we detect peaks larger than 25 µV and take -250 ms and 500 ms
% snippets of data around each event per channel (at ch).

datalfp=ChanData(ChannelNum(ch),:) ;

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
        if out1(msn)-250>0 && out1(msn)+500<size(ChanData,2)
            
            % A step to capture the time span surrounding the detected peaks in the LFP
            TI=(out1(msn)-250:out1(msn)+500);
            % A step to capture waves surrounding detected peaks in the LFP
            plotD3=datalfp(:,TI);
            % Useful for measuring features of the the waveform snippets
            Rangecheck=250-50:250+50;
            
            %Producing a matrix of the waveform snippets along with
            %measures of the amplitudes and timing of the detected waveforms
            WaveClassAll=[WaveClassAll; Deflection ChannelNum(ch) max(plotD3(Rangecheck))...
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
    % event, this step finds the maximum peak in the recording and
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
TimingWaveAll=[TimingWaveAll; TimingrecenterAll(NonrepeatsAll,1:2)...
    repmat([Patient ChannelNum(hs)],size(NonrepeatsAll,1),1)...
    B(AllWaveforms(NonrepeatsAll),1:8) recenteredAll(NonrepeatsAll,:)];

diff1=1;
diff2=20;

% Capturing the differential value at the image onset
diffStep=diff(TimingWaveAll1(:,stepvla+100:stepvla+101)');

% Measuring the amplitude of each waveform relative to the time
% period before the detected waveform peak
Amplitude=max((abs(TimingWaveAll1(:,stepvla+95:stepvla+300)-repmat(nanmean(TimingWaveAll1(:,stepvla+50:stepvla+95),2),1,length(stepvla+95:stepvla+300))))')';

% This step is to remove large amplitude (>2 microV) and small
% amplitude (<50 microV).
FI=find(ValAll<3  & Amplitude<2000 &...
    Amplitude>50 &...
    diffStep'>diff1 &...
    diffStep'<diff2);

Wave1= TimingWaveAll1(FI,1:end);
Wave1(:,stepvla:end)=TimingWaveAll1(FI,stepvla:end)-repmat(nanmean(TimingWaveAll1(FI,stepvla+50:stepvla+95),2),1,size(TimingWaveAll1(FI,stepvla:end),2));
Val=1:size(Wave1,1);

% Some large amplitude events with huge amounts of variation can
% filter through. This step is to remove those waveforms:
FG=find(ValAll<3 & Amplitude<2000 &...
    Amplitude>50 &...
    diffStep'>-diff2 &...
    diffStep'<-diff1);

Wave2= TimingWaveAll1(FG,1:end);
Wave2(:,stepvla:end)=TimingWaveAll1(FG,stepvla:end)-repmat(nanmean(TimingWaveAll1(FG,stepvla+50:stepvla+95),2),1,size(TimingWaveAll1(FG,stepvla:end),2));

% Some large amplitude events with huge amounts of variation can
% filter through. This step is to remove those waveforms:
Variance=var(WaveAll(:,14+10:14+75)');
VarG=WaveAll(Variance<30000,1:end);

% Creating the matrix of data with the cleaned up waveforms.
WaveAllPat=[ fi*ones(size(VarG,1),1) repmat([TimeBrev(1) TimeMax TimeCold(end,:)],size(VarG,1),1)  ValBaseline VarG];

%% Correlating the template waveforms to the waveform snippets
% Loading the waves determined from PCA and kmeans approaches
load(['N:\PEDOTPtNRData\TimingAlignments\CorValTimeTimingWave',num2str(2)],'MeanVar3','MeanVar2','MeanVar')

Type2WaveA=MeanVar;
Type2WaveB=MeanVar2;
Type3Wave=MeanVar3;

save(['N:\PEDOTPtNRData\TimingAlignments\Type2Type3WaveformTemplates'],'Type2WaveA','Type2WaveB','Type3Wave')

Type2WaveA=MeanVar;
Type2WaveB=MeanVar2;
Type3Wave=MeanVar3;

% Multiply by 100 microvolts to get the waveforms to a similar scale
multTimes=[100];

% Normalize the waveforms to the maximal value
Type2WaveA=Type2WaveA/max(Type2WaveA);
Type2WaveB=Type2WaveB/max(Type2WaveB);
Type3Wave=Type3Wave/max(Type3Wave);

% Correlate each waveform snippet with the template waves
HR=[];PR=[];HR2=[];PR2=[];HR3=[];PR3=[];
for smsn=1:size(WaveAllPat,1)
    [rho,p]=corr(multTimes*Type2WaveA',WaveAllPat(smsn,valstep:end)');
    HR(smsn,ampstep)=rho;
    PR(smsn,ampstep)=p;
    [rho,p]=corr(multTimes*Type2WaveB',WaveAllPat(smsn,valstep:end)');
    HR2(smsn,ampstep)=rho;
    PR2(smsn,ampstep)=p;
    [rho,p]=corr(multTimes*Type3Wave',WaveAllPat(smsn,valstep:end)');
    HR3(smsn,ampstep)=rho;
    PR3(smsn,ampstep)=p;
    smsn
end


dir1=dir([DirSaveAllp,PatName,'_',num2str(fi),'_CorMultipleTimeTimingLatestWave_31*']);
if isempty(dir1)==0
    load([DirSaveAllp,PatName,'_',num2str(fi),'_CorMultipleTimeTimingLatestWave_31.mat'])
    CorrSharp1=HR';
    CorrSharp2=HR2';
    CorrSharp3=HR3';
    SharpWave=WaveAllPat;
    MaxVsharp=MAXV;
    MovementTimesONOFFSharp=MovementTimesONOFF;
    Class2fast=MeanVar;
    Class3fast=MeanVar2;
    Class4fast=MeanVar3;
    % Class2faststd=StdVar;
    % Class3faststd=StdVar2;
    % Class4faststd=StdVar3;
end
dir2=dir([DirSaveAllp,PatName,'_',num2str(fi),'_CorMultipleTimeTimingLatestWave_32*']);
if isempty(dir2)==0
    load([DirSaveAllp,PatName,'_',num2str(fi),'_CorMultipleTimeTimingLatestWave_32.mat'])
    CorrSlow1=HR';
    CorrSlow2=HR2';
    CorrSlow3=HR3';
    SlowWave=WaveAllPat;
    MaxVslow=MAXV;
    MovementTimesONOFFSlow=MovementTimesONOFF;
    Class2slow=MeanVar;
    Class3slow=MeanVar2; %same as Class 3
    Class4slow=MeanVar3;
    % Class2slowstd=StdVar;
    % Class3slowstd=StdVar2;
    % Class4slowstd=StdVar3;
end
dir3=dir([DirSaveAllp,PatName,'_',num2str(fi),'_CorMultipleTimeTimingWaveLatestSub_3',num2str(3),'*']);
if isempty(dir3)==0
    load([DirSaveAllp,PatName,'_',num2str(fi),'_CorMultipleTimeTimingWaveLatestSub_3',num2str(3)],...
        'WAP','HR3P','HR2P','HRP','findls','findns','MAXP',...
        'WAPfast2','HR3P2','HR2P2','HRP2','findls2','findns2','MAXP2',...
        'WAPfast1','HR3P1','HR2P1','HRP1','findls1','findns1','MAXP1')
    CorrCL31=HRP';
    CorrCL32=HR2P';
    CorrCL33=HR3P';
    CL3Wave=WAP;
    MaxVCL3=MAXP';
    MovementTimesONOFFSlow=MovementTimesONOFF;
    Class2slow=MeanVar;
    Class3slow=MeanVar2; %same as Class 3
    Class4slow=MeanVar3;
    % Class2slowstd=StdVar;
    % Class3slowstd=StdVar2;
    % Class4slowstd=StdVar3;
end

AllWave=[zeros(size(SlowWave,1),1) SlowWave;ones(size(SharpWave,1),1) SharpWave;2*ones(size(CL3Wave,1),1) CL3Wave];
AllHRClass2=[zeros(size(SlowWave,1),1) CorrSlow1';ones(size(SharpWave,1),1) CorrSharp1';2*ones(size(CL3Wave,1),1) CorrCL31'];
AllHRClass3=[zeros(size(SlowWave,1),1) CorrSlow2';ones(size(SharpWave,1),1) CorrSharp2';2*ones(size(CL3Wave,1),1) CorrCL32'];
AllHRClass4=[zeros(size(SlowWave,1),1) CorrSlow3';ones(size(SharpWave,1),1) CorrSharp3';2*ones(size(CL3Wave,1),1) CorrCL33'];
AllMaxV=[zeros(size(SlowWave,1),1) MaxVslow';ones(size(SharpWave,1),1) MaxVsharp';2*ones(size(CL3Wave,1),1) MaxVCL3'];

AllHRClassCorr=[AllHRClass2 AllHRClass3 AllHRClass4];
AllHRClassCorr1=AllHRClassCorr;
%             if fi==25
AllWave=AllWave(abs(AllHRClassCorr1(:,2))>.8 | abs(AllHRClassCorr1(:,4))>.8 | abs(AllHRClassCorr1(:,6))>.8,:);
AllMaxV=AllMaxV(abs(AllHRClassCorr1(:,2))>.8 | abs(AllHRClassCorr1(:,4))>.8 | abs(AllHRClassCorr1(:,6))>.8,:);
AllHRClassCorr=AllHRClassCorr(abs(AllHRClassCorr1(:,2))>.8 | abs(AllHRClassCorr1(:,4))>.8 | abs(AllHRClassCorr1(:,6))>.8,:);
%             end

Class2slow=MeanVar;
Class3slow=MeanVar2; %same as Class 3
Class4slow=MeanVar3;


ValMxAll=[];
AllWavePer=AllWave(AllWave(:,2)==fi & AllWave(:,12)==cha,:);
AllWavePerFind=find(AllWave(:,2)==fi & AllWave(:,12)==cha);
AllHRClassCorrPer=AllHRClassCorr(AllWave(:,2)==fi & AllWave(:,12)==cha,:);
AllHRMax=AllMaxV(AllWave(:,2)==fi & AllWave(:,12)==cha,:);
if isempty(AllWavePer)==0
    [GR,ID]=grp2idx(AllWavePer(:,9));
    ID=str2num(char(ID));
    
    for ksn=1:length(ID)
        
        ZI=AllHRClassCorrPer(AllWavePer(:,9)==ID(ksn),:);
        if isempty(ZI)==0
            
            ZIWave=AllWavePer(AllWavePer(:,9)==ID(ksn),:);
            FIZi=AllWavePerFind(AllWavePer(:,9)==ID(ksn));
            ZIAllHRMax=AllHRMax(AllWavePer(:,9)==ID(ksn),:);
            if size(ZI,1)==3
                ValMx=[ZI(1,:) FIZi(1) ZI(2,:) FIZi(2) ZI(3,:) FIZi(3)];
            elseif size(ZI,1)==2
                ValMx=[ZI(1,:) FIZi(1) ZI(2,:) FIZi(2) NaN NaN NaN NaN NaN NaN NaN];
            elseif size(ZI,1)==1
                ValMx=[ZI(1,:) FIZi(1) NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
            else
                pause
            end
            ValMxAll=[ValMxAll;repmat([fi cha],size(ValMx,1),1) ValMx ZIAllHRMax(1)];
            %                 pause
            %             ksn
        end
    end
    cha
end
fi
save([DirSaveAllp,PatName,'_',num2str(fi),'AllDataWavesTogether'],'AllWave',...
    'AllHRClassCorr','AllMaxV','AllHRClass3','AllHRClass4','AllHRClass2')
save([DirSaveAllp,PatName,'_',num2str(fi),'ValMxAll_3'],'ValMxAll')
%%
% load('D:\PEDOT\TimingAlignments\ValMxAll.mat')
load([DirSaveAllp,PatName,'_',num2str(fi),'ValMxAll_3'])
Categories=[];
for cs=1:size(ValMxAll,1)
    NF=abs(ValMxAll(cs,[4:2:8 11:2:15 18:2:22]));
    NFInd=abs(ValMxAll(cs,[3:2:7 10:2:14 17:2:21]));
    %     NF
    %     [mx,PerInd]=max(NF);
    %     PerInd
    %     pause
    NF(isnan(NF)==1)=0;
    [mx,PerInd]=max(NF);
    DiffPerInd=0;
    if NFInd(PerInd)==2
        DiffPerInd=4;
    elseif NFInd(PerInd)==1
        DiffPerInd=3;
    elseif NFInd(PerInd)==0
        DiffPerInd=2;
    end
    
    if mx==1
        pause
    end
    vl=[4:2:8 11:2:15 18:2:22];
    Categories(cs,:)=[PerInd DiffPerInd mx ValMxAll(cs,vl(PerInd))];
end
%%
CategoriesSub=Categories;

save([DirSaveAllp,PatName,'_',num2str(fi),'CategoriesSuball3'],'CategoriesSub','Categories')

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

% Identifies the Type 2 and Type 3 events across the recording per channel
% where the onset
for SC=1:2
    
    if SC==1
        % Selecting the waveforms with correlations to the Type 2 template
        % waveform and including the correlation values for later
        % assessment
        ValBit=ValMxAll(CategoriesSub(:,1)<3 &...
            CategoriesSub(:,3)>thresh ,9);
        VarG=[AllWave(ValBit,1:end) CategoriesSub(CategoriesSub(:,1)<3 & CategoriesSub(:,3)>thresh,:)];
        
        % This step keeps only Type 2 waves where the first 90 sample points are,
        % on average, less than 25 microVolts and absolute values less than 500 microVolts
        VarG=VarG(abs(VarG(:,22+97))<500 & nanmean(abs(VarG(:,22+1:22+90)),2)<25,:);
        
        % This step keeps only Type 3 waves where the onset of the
        % waveforms have a second differential which is greater than 2 to
        % only keep fast onset waveforms to differentiate from rhythmic
        % waveforms based on the saline tests
        PW=VarG(VarG(:,2)>-1,22:421);
        % Second differential calculation
        VM=nanmean(abs(diff(diff((PW(:,rangeview))'))));
        % Thresholding for the higher second differential values
        VarG=VarG(VM>2,:);
    elseif SC==2
        % Selecting the waveforms with correlations to the Type 3 template
        % waveform and including the correlation values for later
        % assessment
        ValBit=ValMxAll(CategoriesSub(:,1)==3 &...
            CategoriesSub(:,3)>thresh,9);
        VarG=[AllWave(ValBit,1:end) CategoriesSub(CategoriesSub(:,1)==3 & CategoriesSub(:,3)>thresh,:)];
        
        % This step keeps only Type 3 waves where the first 90 sample points are,
        % on average, less than 25 microVolts
        VarG=VarG(nanmean(abs(VarG(:,22+1:22+90)),2)<25,:);
        
        % This step keeps only Type 3 waves where the onset of the
        % waveforms have a second differential which is greater than 2 to
        % only keep fast onset waveforms to differentiate from rhythmic
        % waveforms based on the saline tests
        PW=VarG(VarG(:,2)>-1,22:421);
        % Second differential calculation
        VM=nanmean(abs(diff(diff((PW(:,rangeview))'))));
        % Thresholding for the higher second differential values
        VarG=VarG(VM>2,:);
    end
    
    Wave1=VarG(:,1:end); %baseline
    
    WaveCount=[WaveCount;Patient SC size(Wave1,1) size(Wave12,1)];
end



