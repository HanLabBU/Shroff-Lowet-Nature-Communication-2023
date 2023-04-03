%% Spike trigerred power plots
%7274886
%  6478633
% 1700991
% 34.3%
% 10.95
% 26.26%
clear all
%close all
%% Load files =[];
Is_ratio=[];Is_ratio2=[];P_ratio=[];cell_id=[];Is_ratio3=[];d_id=[];
  nn=0;nn1=0;
      amount_of_timepoints=0;  amount_of_timepointsMOV=0; 
    COH_LFP=[]; A_LFP=[]; A_Vm=[];
        COH_Vm=[];
%%%%%%%%%%%%%%%
sw_t=3; % switch cell seection, 1==delta, 2==regular, 3=both
i_motion_threshold=0.2; % image motion threshold  0.0650microns/ms
 snr_thres= 0;%
%%%%%%%%%%%%%%%%%
for gh=[1 2 ]  % 1= CHAT, 2=SYN, 
    
    
    if gh==1
        %cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Hua-an_Tseng/Data/SomArchon_Striatum/Linear_movement_PinnedBall/good_data_ALL/ChAT/')
        cd('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata2\Hua-an Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\ChAT\')
        
   
        %%%%%%%%%%%% identification vector, 1=delta, 2=regular/other, 0= low SNR
      dyn_type=   [ 0,2,1,1,1,1,1,1,1,0,1,2,1,1,1,1,1,0,1,1,1,1,1,1,1,2,1,2,0,2,0,1]     
        %%%%%%%%%%%%%%
    else
        %cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Hua-an_Tseng/Data/SomArchon_Striatum/Linear_movement_PinnedBall/good_data_ALL/MSN/')
        cd('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata2\Hua-an Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\MSN\')
   
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
                     windF=length(aligned.trace{in});%windF=ceil(FS*1);
                    clear d A
                    d(1,:)= aligned.lfp{in}(:,1:end); % LFP
                  
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
  
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    %% Make it suitable for fieldtrip by zscoring
                    if mean(spRES{in}.spike_snr{:}) >=snr_thres%  & length([spRES{in}.spike_snr{:}])>10
                        for nh= 1:Nr
                            if  1 %
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
                
                if length(lfp.trial) >0 %  & 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    

                               warning off
    cfg = []; %block_type == cfg.blk
    cfg.method ='wavelet'; %'mvar';
    cfg.output ='fourier';
    cfg.taper='hanning';
    cfg.keeptapers ='yes';
    cfg.keeptrials ='yes';
    cfg.trials='all';cfg.tapsmofrq =5;%
    %cfg.channel= ['all'];
    cfg.channel= [1 2 3]; 
    cfg.foi= [1:100];
    cfg.toi=lfp.time{1}(1:1:end) ;
    cfg.width =5;
    freq2 = ft_freqanalysis(cfg, lfp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
       wavD = angle(squeeze(freq2.fourierspctrm));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
       for tr=1:size(wavD,1)
            imag_mot= lfp.trial{tr}(6,:); 
          wavD(tr,3,:,imag_mot>i_motion_threshold)=NaN; 
           amount_of_timepoints= amount_of_timepoints+length(imag_mot);
            amount_of_timepointsMOV= amount_of_timepointsMOV+ length(find(imag_mot>i_motion_threshold));
       end       
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 CH=1;sp_shift=0;
     %%%%%%%%%%%%
     [PLV_LFP, LFP_A]= spike_field_ppc_adjCHAT(wavD,aligned.spike,CH,sp_shift) ;    
            
                       CH=3;sp_shift=0;
     %%%%%%%%%%%%
     [PLV_Vm, Vm_A]= spike_field_ppc_adjCHAT(wavD,aligned.spike,CH,sp_shift) ;    
  
    COH_LFP=[COH_LFP;PLV_LFP];
        COH_Vm=[COH_Vm;PLV_Vm];
      A_LFP=[A_LFP,LFP_A];
        A_Vm=[A_Vm,Vm_A];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
              
                    cell_id=[cell_id,gh];
                    d_id=[d_id,  dyn_type(sesd) ];
                end
            end
        end
    end
end


cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\git_scripts\fig3\')
savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\git_scripts\fig3\'

pheight=160;


figure('COlor','w','Position', [ 300 300 200 pheight],'Renderer', 'painters')
plot(nanmean(COH_LFP(d_id==1,:),1),'r');
fill_error_area2(freq2.freq,nanmean(COH_LFP(d_id==1,:),1),nanstd(COH_LFP(d_id==1,:))./sqrt(length(find(d_id==1))),[ .5 .5 .5])
hold on,plot(nanmean(COH_LFP(d_id==2,:),1),'b')
fill_error_area2(freq2.freq,nanmean(COH_LFP(d_id==2,:),1),nanstd(COH_LFP(d_id==2,:))./sqrt(length(find(d_id==2))),[ .5 .5 .5])
axis tight
set(gca,'Xscale','log')
 print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'SP_LFP_coherec.pdf'])

figure('COlor','w','Position', [ 300 300 200 pheight],'Renderer', 'painters')
plot(nanmean(COH_Vm(d_id==1,:),1),'r');
fill_error_area2(freq2.freq,nanmean(COH_Vm(d_id==1,:),1),nanstd(COH_Vm(d_id==1,:))./sqrt(length(find(d_id==1))),[ .5 .5 .5])
hold on,plot(nanmean(COH_Vm(d_id==2,:),1),'b')
fill_error_area2(freq2.freq,nanmean(COH_Vm(d_id==2,:),1),nanstd(COH_Vm(d_id==2,:))./sqrt(length(find(d_id==2))),[ .5 .5 .5])
axis tight
set(gca,'Xscale','log')
 print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'SP_Vm_coherec.pdf'])


% figure('COlor','w','Position', [ 300 300 200 pheight],'Renderer', 'painters')
% plot(nanmean(COH_LFP(cell_id==1,:),1),'r');
% fill_error_area2(freq2.freq,nanmean(COH_LFP(cell_id==1,:),1),nanstd(COH_LFP(cell_id==1,:))./sqrt(length(find(cell_id==1))),[ .5 .5 .5])
% hold on,plot(nanmean(COH_LFP(cell_id==2,:),1),'b')
% fill_error_area2(freq2.freq,nanmean(COH_LFP(cell_id==2,:),1),nanstd(COH_LFP(cell_id==2,:))./sqrt(length(find(cell_id==2))),[ .5 .5 .5])
% axis tight
% set(gca,'Xscale','log')

% 
% figure('COlor','w','Position', [ 300 300 200 pheight],'Renderer', 'painters')
% plot(nanmean(COH_Vm(cell_id==1,:),1),'r');
% fill_error_area2(freq2.freq,nanmean(COH_Vm(cell_id==1,:),1),nanstd(COH_Vm(cell_id==1,:))./sqrt(length(find(cell_id==1))),[ .5 .5 .5])
% hold on,plot(nanmean(COH_Vm(cell_id==2,:),1),'b')
% fill_error_area2(freq2.freq,nanmean(COH_Vm(cell_id==2,:),1),nanstd(COH_Vm(cell_id==2,:))./sqrt(length(find(cell_id==2))),[ .5 .5 .5])
% axis tight
% set(gca,'Xscale','log')



 fsel=find(round(freq2.freq)>=2 & round(freq2.freq)<=3 );
V1=nanmean(COH_Vm(d_id==1&cell_id==1,fsel),2);
V2=nanmean(COH_Vm(d_id==2&cell_id==1,fsel),2);
V3=nanmean(COH_Vm(d_id==1&cell_id==2,fsel),2);
V4=nanmean(COH_Vm(d_id==2&cell_id==2,fsel),2);
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
% % print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Vm_pow_delta_reg_CELLIDbarquant.pdf'])



figure('COlor','w','Position', [ 300 400 180 150],'Renderer', 'painters')
violinplotSTR(V1,[1.3 ],'ViolinColor', [ 0.9 0 0.0])
hold on,
violinplot2(V2,[1.7 ],'ViolinColor', [ 0. 0. 0.9])
violinplotSTR(V3,[ 2.2],'ViolinColor', [ 0.9 0. 0.])
violinplotSTR(V4,[ 2.6],'ViolinColor', [ 0. 0. 0.9])
line([ 0.8 2.9], [ 0  0],'COlor', [0 0 0 ],'Linewidth',0.5)
% % print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Vm_pow_delta_reg_CELLIDbarquant.pdf'])
 print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'SP_Vm_coherec_BARQUANT.pdf'])

   [h,p,ci,stats] =ttest2(V1,V2)
%df=25, tst=2.7,       0.0154    0.0055
   [h,p,ci,stats] =ttest2(V3,V4) 
%df=23, tst=6.6,        1.1719e-07    3.0773e-08

%    [h,p,ci,stats] =ttest2(V1,V4)
% %df=36, tst=5.02,        1.8931e-05
%    [h,p,ci,stats] =ttest2(V3,V2)
% %df=12, tst=3.6,          0.0013
%    [h,p,ci,stats] =ttest2(V1,V3)
% %df=29, tst=-1.7490,          0.2013

 fsel=find(round(freq2.freq)>=2  & round(freq2.freq)<=3);
V1=nanmean(COH_LFP(d_id==1&cell_id==1,fsel),2);
V2=nanmean(COH_LFP(d_id==2&cell_id==1,fsel),2);
V3=nanmean(COH_LFP(d_id==1&cell_id==2,fsel),2);
V4=nanmean(COH_LFP(d_id==2&cell_id==2,fsel),2);
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



figure('COlor','w','Position', [ 300 400 180 150],'Renderer', 'painters')
violinplotSTR(V1,[1.3 ],'ViolinColor', [ 0.9 0 0.0])
hold on,
violinplot2(V2,[1.7 ],'ViolinColor', [ 0. 0. 0.9])
violinplotSTR(V3,[ 2.2],'ViolinColor', [ 0.9 0. 0.])
violinplotSTR(V4,[ 2.6],'ViolinColor', [ 0. 0. 0.9])
line([ 0.8 2.9], [ 0  0],'COlor', [0 0 0 ],'Linewidth',0.5)
% % print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Vm_pow_delta_reg_CELLIDbarquant.pdf'])
 print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'SP_LFP_coherec_BARQUANT.pdf'])

   [h,p,ci,stats] =ttest2(V1,V2)
%df=25, tst=2.17,       0.0393    0.0323
   [h,p,ci,stats] =ttest2(V3,V4)
%df=23, tst=4.87,        6.3905e-05     4.8479e-05
%    [h,p,ci,stats] =ttest2(V1,V4)
% %df=36, tst=3.72,      6.6058e-04
%    [h,p,ci,stats] =ttest2(V3,V2)
% %df=12, tst=2.7,         0.0187
%    [h,p,ci,stats] =ttest2(V1,V3)
% %df=29, tst=-3.0,         0.0054

figure,rose(A_LFP(3,d_id==1&cell_id==1),12);hold on,
rose(A_LFP(3,d_id==1&cell_id==2),12)

figure
polarscatter(A_Vm(2,d_id==1&cell_id==1),COH_Vm(2,d_id==1&cell_id==1));hold on,
polarscatter(A_Vm(2,d_id==1&cell_id==2),COH_Vm(2,d_id==1&cell_id==2));hold on,




figure('COlor','w','Position', [ 300 300 200 pheight],'Renderer', 'painters')
[n1 n2]=hist(circ_mean(A_LFP([ 2  ],d_id==1)),[-pi:0.5:pi])
bar([ n2-pi*2 n2 n2+pi*2],[n1./sum(n1) n1./sum(n1) n1./sum(n1)],1.5)
hold on, plot([-pi.*2:0.1:pi.*2], cos([-pi.*2:0.1:pi.*2])./10+0.1,'k' ,'Linewidth',1)
axis tight
xlim([-pi-2 pi+2])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'LFP_phase_histt.pdf'])


figure('COlor','w','Position', [ 300 300 200 pheight],'Renderer', 'painters')
[n1 n2]=hist(circ_mean(A_Vm([ 2  ],d_id==1)),[-pi:0.15:pi])
bar([ n2-pi*2 n2 n2+pi*2],[n1./sum(n1) n1./sum(n1) n1./sum(n1)],2,'r')
hold on, plot([-pi.*2:0.1:pi.*2], cos([-pi.*2:0.1:pi.*2])./7+(1/7),'k' ,'Linewidth',1)
axis tight
xlim([-pi-2 pi+2])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Vm_phase_histt.pdf'])


% figure
% polarscatter(A_LFP(2,d_id==1&cell_id==1),COH_LFP(2,d_id==1&cell_id==1)+0.002);hold on,
% polarscatter(A_LFP(2,d_id==1&cell_id==2),COH_LFP(2,d_id==1&cell_id==2)+0.002);hold on,


% figure('COlor','w','Position', [ 300 300 200 pheight],'Renderer', 'painters')
% plot(nanmean(COH_Vm(cell_id==1,:),1));
% hold on,plot(nanmean(COH_Vm(cell_id==2,:),1))
% 
% % 
% figure('COlor','w','Position', [ 300 300 200 pheight],'Renderer', 'painters')
% plot(freq2.freq,nanmean(Is_ratio3(d_id==1,:),1),'r') ,
% fill_error_area2(freq2.freq,nanmean(Is_ratio3(d_id==1,:),1),nanstd(Is_ratio3(d_id==1,:))./sqrt(length(find(d_id==1))),[ .5 .5 .5])
% plot(freq2.freq,nanmean(Is_ratio3(d_id==2,:),1),'b') ,
% fill_error_area2(freq2.freq,nanmean(Is_ratio3(d_id==2,:),1),nanstd(Is_ratio3(d_id==2,:))./sqrt(length(find(d_id==2))),[ .5 .5 .5])
% axis tight
% xlim([ 1 70])
% set(gca,'Xscale','log')
%  print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Vm_pow_delta_reg.pdf'])
% 
% fsel=find(freq2.freq>1  & freq2.freq<4);
% V1=nanmean(Is_ratio3(d_id==1,fsel),2);
% V2=nanmean(Is_ratio3(d_id==2,fsel),2);
% 
%   figure('COlor','w','Position', [ 300 400 180 pheight],'Renderer', 'painters') 
% b1=bar(1,nanmean(V1),'Facecolor',[ 0.8 0 0]);hold on,
% set(b1,'FaceAlpha',0.4)
%   b1=bar(2,nanmean(V2),'Facecolor',[ 0 0 0.8]);hold on,
% set(b1,'FaceAlpha',0.4)
% errorbar([1   ],nanmean(V1,1), nanstd(V1)./sqrt(size(V1,1)),'.k');
% errorbar([2   ],nanmean(V2,1), nanstd(V2)./sqrt(size(V2,1)),'.k');
% ylim([0.15 0.3])
% [h,p,ci,stats] =ttest2(V1,V2)
% %  0.0011, df=42
% % print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Vm_pow_delta_reg_barquant.pdf'])
% 
% 
%  
%  
%  fsel=find(freq2.freq>1  & freq2.freq<4);
% V1=nanmean(Is_ratio3(d_id==1&cell_id==1,fsel),2);
% V2=nanmean(Is_ratio3(d_id==2&cell_id==1,fsel),2);
% V3=nanmean(Is_ratio3(d_id==1&cell_id==2,fsel),2);
% V4=nanmean(Is_ratio3(d_id==2&cell_id==2,fsel),2);
%   figure('COlor','w','Position', [ 300 400 180 pheight],'Renderer', 'painters') 
% b1=bar(1,nanmean(V1),'Facecolor',[ 0.8 0 0]);hold on,
% set(b1,'FaceAlpha',0.4)
%   b1=bar(2,nanmean(V2),'Facecolor',[ 0 0 0.8]);hold on,
% set(b1,'FaceAlpha',0.4)
%   b1=bar(4,nanmean(V3),'Facecolor',[ 0.8 0 0]);hold on,
% set(b1,'FaceAlpha',0.4)
%   b1=bar(5,nanmean(V4),'Facecolor',[ 0 0 0.8]);hold on,
% set(b1,'FaceAlpha',0.4)
% errorbar([1   ],nanmean(V1,1), nanstd(V1)./sqrt(size(V1,1)),'.k');
% errorbar([2   ],nanmean(V2,1), nanstd(V2)./sqrt(size(V2,1)),'.k');
% errorbar([4   ],nanmean(V3,1), nanstd(V3)./sqrt(size(V3,1)),'.k');
% errorbar([5   ],nanmean(V4,1), nanstd(V4)./sqrt(size(V4,1)),'.k');
% 
% ylim([0.15 0.33])
% [h,p,ci,stats] =ttest2(V1,V2)
% %  0.0011, df=42
% % print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Vm_pow_delta_reg_CELLIDbarquant.pdf'])
% 
