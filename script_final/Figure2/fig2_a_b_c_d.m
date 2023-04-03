%% Spike trigerred power plots

clear all
%close all
%% Load files =[];
Is_ratio=[];Is_ratio2=[];P_ratio=[];cell_id=[];Is_ratio3=[];d_id=[];
  nn=0;nn1=0;
%%%%%%%%%%%%%%%
sw_t=3; % switch cell seection, 1==delta, 2==regular, 3=both
snr_thres= 5;%
i_motion_threshold=0.2; % image motion threshold

%%%%%%%%%%%%%%%%%
for gh=[1 2]  % 1= CHAT, 2=SYN, 
 
 
    
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
dyn_type= [2,2,2,2,0,2,2,0,2,2,0,1,1,2,0,1,2,1,1,2,0,1,1,1,1,2,2,2,2,0,2]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    
    
    ses= dir('*.mat')
    
    for sesd=[ 1:length(ses)]
        sesd
        load(ses(sesd).name) %LOAD
        
        try
            
            if length(aligned.trace)>0% Actual - low movement
                
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
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  d(6,:)=  aligned.imgmotion{in};  % combined X-Y image motion
  d(6,[1 2 ])=d(6,3); % remove edge effect
  d(6,end:-1:end)=d(6,end-2); % remove edge effect
  
  d(6,:)= [ 0 fastsmooth(abs(hilbert(diff(fastsmooth(d(6,:),10,1,1)))),100,1,1)];
  %%%%%%
                    
                    %% Make it suitable for fieldtrip by zscoring
                    if mean(spRES{in}.spike_snr{:}) >snr_thres  %& length([spRES{in}.spike_snr{:}])>10
                        for nh= 1:Nr
                           if  mean(d(6,(1  + ((windF)* (nh-1)):   windF * (nh))))<i_motion_threshold
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
                
                if length(lfp.trial) >0 %  
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    
                    
                    
                    
                               warning off
    cfg = []; %block_type == cfg.blk
    cfg.method ='mtmfft'; %'mvar';
    cfg.output ='fourier';
    cfg.taper='hanning';
    cfg.keeptapers ='yes';
    cfg.keeptrials ='yes';
    cfg.trials='all';cfg.tapsmofrq =5;%
    %cfg.channel= ['all'];
    cfg.channel= [1 2 3]; 
    cfg.foi= [1:100];
   
    freq2 = ft_freqanalysis(cfg, lfp);
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
       wavA = abs(squeeze(freq2.fourierspctrm));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
%        for tr=1:size(wavA,1)
%             imag_mot= lfp.trial{tr}(6,:); 
%           wavA(tr,3,:,imag_mot>i_motion_threshold)=NaN; 
%        end       
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%               
                    
                
            
                    Is_ratio3=[Is_ratio3; (squeeze((nanmean(wavA(:,3,:),1)))')];
                    
              
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



figure('COlor','w','Position', [ 300 300 200 pheight],'Renderer', 'painters')
plot(freq2.freq,nanmean(Is_ratio3(d_id==1,:),1),'r') ,
fill_error_area2(freq2.freq,nanmean(Is_ratio3(d_id==1,:),1),nanstd(Is_ratio3(d_id==1,:))./sqrt(length(find(d_id==1))),[ .5 .5 .5])
plot(freq2.freq,nanmean(Is_ratio3(d_id==2,:),1),'b') ,
fill_error_area2(freq2.freq,nanmean(Is_ratio3(d_id==2,:),1),nanstd(Is_ratio3(d_id==2,:))./sqrt(length(find(d_id==2))),[ .5 .5 .5])
axis tight
xlim([ 1 70])
set(gca,'Xscale','log')
 print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Vm_pow_delta_reg.pdf'])

fsel=find(freq2.freq>1  & freq2.freq<4);
V1=nanmean(Is_ratio3(d_id==1,fsel),2);
V2=nanmean(Is_ratio3(d_id==2,fsel),2);

  figure('COlor','w','Position', [ 300 400 180 pheight],'Renderer', 'painters') 
b1=bar(1,nanmean(V1),'Facecolor',[ 0.8 0 0]);hold on,
set(b1,'FaceAlpha',0.4)
  b1=bar(2,nanmean(V2),'Facecolor',[ 0 0 0.8]);hold on,
set(b1,'FaceAlpha',0.4)
errorbar([1   ],nanmean(V1,1), nanstd(V1)./sqrt(size(V1,1)),'.k');
errorbar([2   ],nanmean(V2,1), nanstd(V2)./sqrt(size(V2,1)),'.k');
ylim([0.15 0.3])
[h,p,ci,stats] =ttest2(V1,V2)
%  0.0011, df=42
 print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Vm_pow_delta_reg_barquant.pdf'])


 
 
 fsel=find(freq2.freq>1  & freq2.freq<4);
V1=nanmean(Is_ratio3(d_id==1&cell_id==1,fsel),2);
V2=nanmean(Is_ratio3(d_id==2&cell_id==1,fsel),2);
V3=nanmean(Is_ratio3(d_id==1&cell_id==2,fsel),2);
V4=nanmean(Is_ratio3(d_id==2&cell_id==2,fsel),2);
  figure('COlor','w','Position', [ 300 400 180 pheight],'Renderer', 'painters') 
b1=bar(1,nanmean(V1),'Facecolor',[ 0.8 0 0]);hold on,
set(b1,'FaceAlpha',0.4)
  b1=bar(2,nanmean(V2),'Facecolor',[ 0 0 0.8]);hold on,
set(b1,'FaceAlpha',0.4)
  b1=bar(4,nanmean(V3),'Facecolor',[ 0.8 0 0]);hold on,
set(b1,'FaceAlpha',0.4)
  b1=bar(5,nanmean(V4),'Facecolor',[ 0 0 0.8]);hold on,
set(b1,'FaceAlpha',0.4)
errorbar([1   ],nanmean(V1,1), nanstd(V1)./sqrt(size(V1,1)),'.k');
errorbar([2   ],nanmean(V2,1), nanstd(V2)./sqrt(size(V2,1)),'.k');
errorbar([4   ],nanmean(V3,1), nanstd(V3)./sqrt(size(V3,1)),'.k');
errorbar([5   ],nanmean(V4,1), nanstd(V4)./sqrt(size(V4,1)),'.k');
  figure('COlor','w','Position', [ 300 400 200 150],'Renderer', 'painters')
violinplotSTR(V1,[1.3 ],'ViolinColor', [ 0.9 0 0.0])
hold on,
violinplot2(V2,[1.7 ],'ViolinColor', [ 0. 0. 0.9])
violinplotSTR(V3,[ 2.2],'ViolinColor', [ 0.9 0. 0.])
violinplotSTR(V4,[ 2.6],'ViolinColor', [ 0. 0. 0.9])
line([ 0.8 2.9], [ 0  0],'COlor', [0 0 0 ],'Linewidth',0.5)
ylim([0.15 0.33])
 print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Vm_pow_delta_reg_CELLIDbarquant.pdf'])

[h,p,ci,stats] =ttest2(V1,V2)
%  0.0475, df=25
[h,p,ci,stats] =ttest2(V3,V4)
%   0.0064, df=23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


selC1=d_id==1&cell_id==1;
   selC2=d_id==2&cell_id==1 ;
  selC3=d_id==1&cell_id==2  ;
    selC4=d_id==2&cell_id==2  ;

figure('COlor','w','Position', [ 300 300 200 pheight],'Renderer', 'painters')
plot(freq2.freq,nanmean(Is_ratio3( selC1,:),1),'r') ,
fill_error_area2(freq2.freq,nanmean(Is_ratio3( selC1,:),1),nanstd(Is_ratio3( selC1,:))./sqrt(length(find( selC1))),[ .5 .5 .5])
plot(freq2.freq,nanmean(Is_ratio3( selC2,:),1),'b') ,
fill_error_area2(freq2.freq,nanmean(Is_ratio3( selC2,:),1),nanstd(Is_ratio3( selC2,:))./sqrt(length(find( selC2))),[ .5 .5 .5])
axis tight
xlim([ 1 70])
set(gca,'Xscale','log')
 print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Vm_pow_CHAT_all_' '.pdf'])


figure('COlor','w','Position', [ 300 300 200 pheight],'Renderer', 'painters')
plot(freq2.freq,nanmean(Is_ratio3( selC3,:),1),'r') ,
fill_error_area2(freq2.freq,nanmean(Is_ratio3( selC3,:),1),nanstd(Is_ratio3( selC3,:))./sqrt(length(find( selC3))),[ .5 .5 .5])
plot(freq2.freq,nanmean(Is_ratio3( selC4,:),1),'b') ,
fill_error_area2(freq2.freq,nanmean(Is_ratio3( selC4,:),1),nanstd(Is_ratio3( selC4,:))./sqrt(length(find( selC4))),[ .5 .5 .5])
axis tight
xlim([ 1 70])
set(gca,'Xscale','log')
 print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Vm_pow_MSN_all_' '.pdf'])

