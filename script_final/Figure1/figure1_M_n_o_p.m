%% Spike trigerred power plots

clear all
%close all
%% Load files =[];
Is_ratio=[];Is_ratio2=[];P_ratio=[];cell_id=[];firing_rate=[];d_id=[];
  nn=0;nn1=0;
%%%%%%%%%%%%%%%
sw_t=3; % switch cell seection, 1==delta, 2==regular, 3=both
%%%%%%%%%%%%%%%%%
for gh=[ 1 2]  % 1= CHAT, 2=SYN, 
    
    %MSN, 1, low SNR
    %MSN, 2, high SNR,regular
    %MSN 3., low SNR, ?
    %MSN, 4, high SNR, regular, fast-spiking
    %MSN 5, low SNR, ?
    %MSN 6, high SNR, regular
    %MSN 7, mid/good, SNR, no structure, regular
    %MSN 8, low SNR,  high motion,,
    %MSN 9, good SNR,  high motion, no delta, theta bursty?
    %MSN 10= theta
    %MSN 11, low SNR, very high motion, no delta,high firing
    %MSN 12, mid SNR, high motion, delta,
    %MSN 13, mid SNR, high motion, delta,
    %MSN 14, god SNR, low motion, no delta, high firing
    %MSN 15, mid SNR, very high motion, delta, burst
    %MSN 16, mid SNR, very high motion, delta, burst
    %MSN 17, severe motion
    %MSN 18, delta, good SNR, sign motion
    %MSN 19, good/mid SNR, high motion, delta, burst
    %MSN 20, low/mid SNR, high motion
    %MSN 21, low SNR
    %MSN 22, delta, good SNR
    %MSN 23, delta, mid SNR
    %MSN 24/25= strong delta
    %MSN 26= no structure, good snr, regular, some theta burst?
    %MSN 27= mid snr, no delta, motion
    %MSN 28= mid snr, no delta, motion
    %MSN 29, mid snr, no delta, high fire periods, lot of motion
    %MSN 30 low snr, motion high
    %MSN 31 mid snr, motion high
    %MSN,32, low SNR, delta?
    %MSN,33, low SNR,?
    %MSN 34. low/mid SNR, no spikes ISI? what happended
    %MSN 35. low/mid SNR, no spikes ISI? what happended
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %CH1, good SNR, delta, single spikes, low trial number,
    %CH2, good SNR, low delta and little bursting
    %CH3, good SNR, strong delta, bursty
    %CH4, good SN, mid delta,
    %CH5 good SNR, strong delta,bursty
    %CH6, mid/good SNR, strong delta,
    %CH7,mid/good SNR, strong delta,
    %CH8, super SNR, mid-strong delta
    %CH9, mid/good SNR, strong delta,
    %CH10, good SNR, strong delta,
    %CH11, good SNR, strong delta,
    %CH12, good SNR, REGULAR,
    %CH13, good SNR, strong delta, low intra-burst freq
    %CH14, mid/good SNR, strong delta,
    %CH15, mid SNR,regular
    %CH16, good SNR, strong delta,
    %CH17, good SNR, strong delta,
    %CH18, low SNR, regular?
    %CH19, good SNR, strong delta,
    %CH20, good SNR, strong delta,
    %CH21, mid SNR, strong delta, motion
    %CH22, low SNR, delta?
    %CH23, good SNR, strong delta
    %CH24, good SNR, strong delta
    %CH25, low/mid SNR, some delta, look more regular
    %CH26, good SNR, no delta,regular
    %CH27, good SNR, strong delta,mixed?
    %CH28, mid/good SNR, regular
    %CH29, mid SNR, delta, motion
    %CH30, mid/good SNR, little delta, more regular
    %CH31, mid SNR, some delta, sparse
    
    
    if gh==1
        %cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Hua-an_Tseng/Data/SomArchon_Striatum/Linear_movement_PinnedBall/good_data_ALL/ChAT/')
        cd('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata2\Hua-an Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\ChAT\')
        
        if sw_t==1 % delta-identified
            chs=[3,4,5,6,7,8,9,10,11,13,14,16,17,19,20,21,23,24,27,29];
        elseif sw_t==2 % regular-identified 
            chs=  [ 2 12 15  26 28  ];
        elseif sw_t==3   % all
            chs=[2 12 15  26 28,3,4,5,6,7,8,9,10,11,13,14,16,17,19,20,21,23,24,27,29];
            
        end
        %%%%%%%%%%%% identification vector, 1=delta, 2=regular/other, 0= low SNR
      dyn_type=   [ 0,2,1,1,1,1,1,1,1,0,1,2,1,1,1,1,1,0,1,1,1,1,1,1,1,2,1,2,0,2,0,1] 
 
        %%%%%%%%%%%%%%
    else
        %cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Hua-an_Tseng/Data/SomArchon_Striatum/Linear_movement_PinnedBall/good_data_ALL/MSN/')
        cd('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata2\Hua-an Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\MSN\')
        if sw_t==1      % delta-identified
            chs=[[ 12 13  15 16  18 19  22 23 24 25]];
        elseif sw_t==2     % regular-identified or theta bursting
            chs=  [  2 4 6 7 14 26 27 28 29 9 10 ]
        elseif sw_t==3      % all
            chs=[[ 2 4 6 7 14 26 27 28 29 9 10 12 13  15 16  18 19  22 23 24 25]];
            
        end
        %%%%%%%%%%%% identification vector, 1=delta, 2=regular/other, 0= low SNR

dyn_type= [2,2,2,2,0,2,2,0,2,2,0,1,1,2,0,1,2,1,1,2,0,1,1,1,1,2,2,2,2,0,2,0  0 0 0 0]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    
    ses= dir('*.mat')
    
    for sesd=[ 1:length(ses)]
        sesd
        load(ses(sesd).name) %LOAD
        
        try
            snr_thres= 5;%
            if length(aligned.trace)>1% Actual - low movement
                
                FS=1000;
                lfp=[];mk=0;
                SN=[];
                
                %% Put in a matrix d - LFP, spikes, Vm, High motion and Low motion
                
                for in =1:length(aligned.trace)
                    windF=ceil(FS*1);
                    clear d A
                    d(1,:)= aligned.lfp{in}(:,1:end); % LFP
                    A=aligned.spike{in}(:,1:end);A(isnan(A))=0;
                    A= spRES{in}.roaster;
                    A=aligned.spike{in}(:,1:end);
                    A(isnan(A))=0;
                    d(2,:)=A;  %% SPIKES
                   % d(3,:)= aligned.trace_nosp{in}(:,1:end);  %Vm
                        LFPtr=   aligned.trace_nosp{in}(:,1:end);
                                  Fn = 1000/2;FB=[ 83 86];
                                   [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                                   LFPg= ((filtfilt(B,A,  LFPtr)));
                                      LFPtr= LFPtr-LFPg;
                    d(3,:)=    LFPtr;%aligned.trace_nosp{in}(:,1:end);  %Vm
                    
                    SN= [SN;spRES{in}.spike_snr{:}];;
                    % d(1,:)= dataM.trial{in}(1,9:end);
                    Nr= floor(length(d)./windF);
                    
                    
                    
                    %% GET LOCMOMOTION INFO
                    A=zeros( 1,length(aligned.mvmt{in}(1:end)));
                    HM= aligned.f_mvmt{in}(1:end);
                    HM(HM>=50)=1000;
                    HM(HM<50)=-1000;
                    d(4,:)=HM;
                    d(5,:)=HM.*-1;
                    A(aligned.vectors{in,8})=1;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    d(6,:)=  aligned.imgmotion{in};  % combined X-Y image motion
                    d(6,[1 2 ])=d(6,3); % remove edge effect
                    d(6,end:-1:end)=d(6,end-2); % remove edge effect
                    
                    d(6,:)= [ 0 fastsmooth(abs(hilbert(diff(fastsmooth(d(6,:),10,1,1)))),100,1,1)];
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    %% Make it suitable for fieldtrip by zscoring
                    if mean(spRES{in}.spike_snr{:}) >snr_thres % & length([spRES{in}.spike_snr{:}])>10
                        for nh= 1:Nr
                            if  sum((d(4,(1  + ((windF)* (nh-1)):   windF * (nh))   ))>0)<100
                            mk=mk+1;
                            lfp.trial{mk}(1:3,:) = zscore(d(1:3,(1  + ((windF)* (nh-1)):   windF * (nh))   ),[],2);
                            lfp.trial{mk}(4:5,:) = (d(4:5,(1  + ((windF)* (nh-1)):   windF * (nh))   )); % locomotion
                            lfp.trial{mk}(6,:) = (d(6,(1  + ((windF)* (nh-1)):   windF * (nh))   ));% image motion
                            lfp.time{mk}= (1:size(lfp.trial{mk},2))./FS;end
                        end
                    end
                end
                n=0;
                for ind=1:size(lfp.trial{1},1)
                    n=n+1;
                    lfp.label(n)= { [ 'CH' num2str(ind)]};
                end
                
                if length(lfp.trial) >1 %  & mean(SN) >  5 & length(SN)>5
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    ff=[];
                    for i=1:length(lfp.trial)
                     ff=[ff, sum(lfp.trial{i}(2,:)>0)];
                    end
                    firing_rate=[firing_rate; mean(ff)];
                    
              
                    cell_id=[cell_id,gh];
                    d_id=[d_id,  dyn_type(sesd) ];
                end
            end
        end
    end
end


cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\git_scripts\fig2\')
savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\git_scripts\fig2\'

pheight=160;




V1=firing_rate(cell_id==1& d_id==1);
V2=firing_rate(cell_id==1& d_id==2);
V3=firing_rate(cell_id==2& d_id==1);
V4=firing_rate(cell_id==2& d_id==2);

%ylim([0.15 0.3])
[h,p,ci,stats] =ttest2(V1,V2)
  figure('COlor','w','Position', [ 300 400 200 150],'Renderer', 'painters')
violinplotSTR(V1,[1.3 ],'ViolinColor', [ 0.9 0 0.0])
hold on,
violinplot2(V2,[1.7 ],'ViolinColor', [ 0. 0. 0.9])
violinplotSTR(V3,[ 2.2],'ViolinColor', [ 0.9 0. 0.])
violinplotSTR(V4,[ 2.6],'ViolinColor', [ 0. 0. 0.9])
line([ 0.8 2.9], [ 0  0],'COlor', [0 0 0 ],'Linewidth',0.5)
 print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'FR_rate_delta_reg_CELLID_barquant.pdf'])

 
  [h,p,ci,stats] =ttest2(V1,V2)
 %df=25,   0.2611
  
 [h,p,ci,stats] =ttest2(V3,V4)
  %df=23,    0.0398

 
   

   
   V1=firing_rate(cell_id==1& d_id>0);
V2=firing_rate(cell_id==2& d_id>0);
V3=firing_rate( d_id==1);
V4=firing_rate(d_id==2);
  figure('COlor','w','Position', [ 300 400 200 150],'Renderer', 'painters')
violinplotSTR(V1,[1.3 ],'ViolinColor', [ 0.9 0 0.0])
hold on,
violinplot2(V2,[1.7 ],'ViolinColor', [ 0. 0. 0.9])
violinplotSTR(V3,[ 2.2],'ViolinColor', [ 0.9 0. 0.])
violinplotSTR(V4,[ 2.6],'ViolinColor', [ 0. 0. 0.9])
