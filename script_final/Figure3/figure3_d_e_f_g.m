%% Spike trigerred power plots 

clear all
%close all
%% Load files 
cell_id=[];Is_ratio3=[];d_id=[];
pT=[];pT2=[];  nn=0;nn1=0;
CHANNEL=3;sw_t=3;
i_motion_threshold=0.05; % image motion threshold microns/ms
for gh=[1 2]
    
  
    if gh==1
        %cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Hua-an_Tseng/Data/SomArchon_Striatum/Linear_movement_PinnedBall/good_data_ALL/ChAT/')
        cd('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata2\Hua-an Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\ChAT\')
        
%%%% identification vector, 1=delta, 2=regular/other, 0= low SNR
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

for sesd=[ 1:length(ses)] %[1:25 27:32]%length(ses);%[1:15 18:31] %:length(ses)  %22 exclude 16 17
    %sesd =5;
 
    load(ses(sesd).name)
try
%     
%    for id=1:length( aligned.trace_nosp)
%        
%        
%        
%        
%    end
   snr_thres= 0;%6.5;
% Only consider sessions with more than 3 trials 
if length(aligned.trace)>0

   FS=1000;
   lfp=[];mk=0;
   SN=[];

%% Put in a matrix d - LFP, spikes, Vm, High motion and Low motion 

    for in =1:length(aligned.trace)
        windF=length(aligned.trace{in});%ceil(FS*1);
        clear d A
        d(1,:)= aligned.lfp{in}(:,1:end); % LFP
          A=aligned.spike{in}(:,1:end);A(isnan(A))=0;
        A= spRES{in}.roaster;
          A=aligned.spike{in}(:,1:end);
        A(isnan(A))=0;
        d(2,:)=A;  %% SPIKES
        
    %    d(3,:)= aligned.trace_nosp{in}(:,1:end);  %Vm
           LFPtr=   aligned.trace_nosp{in}(:,1:end);
                                  Fn = 1000/2;FB=[ 82 87];
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
    sc= (6.5 ./40).*2;

  d(6,:)=  aligned.imgmotion{in}.*sc;  % combined X-Y image motion
  d(6,[1 2 ])=d(6,3); % remove edge effect
  d(6,end:-1:end)=d(6,end-2); % remove edge effect
  
  d(6,:)= [ 0 fastsmooth(abs((diff(fastsmooth(d(6,:),10,1,1)))),100,1,1)];
   d(7,:)= aligned.trace{in}(:,1:end);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%

 %% Make it suitable for fieldtrip by zscoring
    if mean(spRES{in}.spike_snr{:}) >snr_thres
   for nh= 1:Nr
       mk=mk+1;
       lfp.trial{mk}(1:3,:) = zscore(d(1:3,(1  + ((windF)* (nh-1)):   windF * (nh))   ),[],2);
       
       lfp.trial{mk}(4:5,:) = (d(4:5,(1  + ((windF)* (nh-1)):   windF * (nh))   )); % locomotion
       lfp.trial{mk}(6,:) = (d(6,(1  + ((windF)* (nh-1)):   windF * (nh))   ));% image motion
lfp.trial{mk}(7,:) = (d(7,(1  + ((windF)* (nh-1)):   windF * (nh))   ));% image motion

       lfp.time{mk}= (1:size(lfp.trial{mk},2))./FS;
   end
  end
end
   n=0;
for ind=1:size(lfp.trial{1},1)
    n=n+1;
    lfp.label(n)= { [ 'CH' num2str(ind)]};
end

                
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 % Wavelet spectrum
 
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
  
    cfg.width =6;
   
    cfg.t_ftimwin =ones(1,length(cfg.foi))*0.15;
    freq2 = ft_freqanalysis(cfg, lfp);
    
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 
 
    wind=300;
    clear aM aM2  
    mm=0;   mm2=0; mm3=0;
    for id=1:size(  freq2.fourierspctrm,1)
        % Confirm 3 here -  Vm, no difference between M and MG ?
        M=abs(squeeze(freq2.fourierspctrm(id,CHANNEL,:,:)))      ;
        imag_mot= lfp.trial{id}(6,:); % GET image motion data
        M(:,  imag_mot>i_motion_threshold)=NaN;  % Segements above threshold are NaN
        
        
        MG=abs(squeeze(freq2.fourierspctrm(id,CHANNEL,:,1:end)))      ;
        
        MG(:,  imag_mot>i_motion_threshold)=NaN;% Segements above threshold are NaN
         M1= lfp.trial{id}(7,:);
        
        s= lfp.trial{id}(2,:)>0;s=find(s); % spike times
        
        % High movement and low movement frames 
        HT=find(lfp.trial{id}(4,:)>0);LT=find(lfp.trial{id}(5,:)>0);
        NT=length(s); TL= size(MG,2);
    
        % Actual- high movement 
        for id2=1:length(s)
            if s(id2) >wind  & s(id2) +wind < length(M) & any(s(id2) == HT)
                mm=mm+1;
                aM(:,:,mm)=   M(:,s(id2)-wind:s(id2)+wind);
                %aM_pre(:,:,mm)= M(:,s(id2)-200:s(id2)-100);
            end      
        end
        
        
    % Actual - low movement
    for id2=1:length(s)
        if s(id2) >wind  & s(id2) +wind < length(M) & any(s(id2) == LT)
            mm2=mm2+1;
            aM2(:,:,mm2)=   M(:,s(id2)-wind:s(id2)+wind);
            %aM3_pre(:,:,mm2)= M(:,s(id2)-200:s(id2)-100);
        end
    end
    
      for id2=1:length(s)
        if s(id2) >wind  & s(id2) +wind < length(M) & any(s(id2) == LT)
            mm3=mm3+1;
            aM3(:,mm3)=   M1(s(id2)-wind:s(id2)+wind);
            %aM3_pre(:,:,mm2)= M(:,s(id2)-200:s(id2)-100);
        end
    end

end
    

if size(aM,3)>10 & size(aM2,3)>10
    nn=nn+1;
    % Mean across all windows
    allC(:,:,nn)=  nanmean(abs(aM),3);    
end
if size(aM,3)>10 & size(aM2,3)>10
    nn1=nn1+1;
    allC2(:,:,nn1)=  nanmean(abs(aM2),3);
     allC3(:,nn1)=  nanmean((aM3),2);
 
      cell_id=[cell_id,gh];
                    d_id=[d_id,  dyn_type(sesd) ];end
end
end
end
end
%cd('./Shuffle mat files/Figure 3/Pre_spike_window_norm2_sub_by_union_no_abs/Median')

cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\git_scripts\fig3\')
savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\git_scripts\fig3\'
pheight=160;

tt=-wind:wind;
tsel=find(tt> -201 & tt <-100);
pre_window= allC2(:,tsel,d_id==1);
pre_value2=nanmean(pre_window,2);
figure('COlor','w','Position', [ 300 300 200 pheight],'Renderer', 'painters')
A=bsxfun(@minus, allC2(:,:,d_id==1), pre_value2);B=bsxfun(@plus, allC2(:,:,d_id==1), pre_value2);
imagesc(-wind:wind,freq2.freq,smooth2a(nanmean(bsxfun(@rdivide,A,B),3),1,20));colormap(jet)
axis xy
%hold on,,plot(-wind:wind,zscore(allC3).*10+20,'k','COlor',[0.4 0.4 0.4])
set(gca,'CLim',[-0.1 0.1]) % -+10% modulation
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Vm_sptrig_pow_DELTA.pdf'])

tt=-wind:wind;
tsel=find(tt> -201 & tt <-100);
pre_window= allC2(:,tsel,d_id==2);
pre_value2=nanmean(pre_window,2);
figure('COlor','w','Position', [ 300 300 200 pheight],'Renderer', 'painters')
A=bsxfun(@minus, allC2(:,:,d_id==2), pre_value2);B=bsxfun(@plus, allC2(:,:,d_id==2), pre_value2);
imagesc(-wind:wind,freq2.freq,smooth2a(nanmean(bsxfun(@rdivide,A,B),3),1,20));colormap(jet)
axis xy
%hold on,,plot(-wind:wind,zscore(allC3).*10+20,'k','COlor',[0.4 0.4 0.4])
set(gca,'CLim',[-0.1 0.1]) % -+10% modulation
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Vm_sptrig_pow_REG.pdf'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 fsel=find(freq2.freq>20  & freq2.freq<40);
 selC=cell_id==1&d_id==1;
pre_window= allC2(:,tsel, selC);pre_value2=nanmean(pre_window,2);
A=bsxfun(@minus, allC2(:,:, selC), pre_value2);B=bsxfun(@plus, allC2(:,:, selC), pre_value2);
MM=bsxfun(@rdivide,A,B);
V1=squeeze(nanmean(nanmean(MM(fsel,wind-50:wind+50,:),2),1));
 selC=cell_id==1&d_id==2;
pre_window= allC2(:,tsel, selC);pre_value2=nanmean(pre_window,2);
A=bsxfun(@minus, allC2(:,:, selC), pre_value2);B=bsxfun(@plus, allC2(:,:, selC), pre_value2);
MM=bsxfun(@rdivide,A,B);
V2=squeeze(nanmean(nanmean(MM(fsel,wind-50:wind+50,:),2),1));
 selC=cell_id==2&d_id==1;
pre_window= allC2(:,tsel, selC);pre_value2=nanmean(pre_window,2);
A=bsxfun(@minus, allC2(:,:, selC), pre_value2);B=bsxfun(@plus, allC2(:,:, selC), pre_value2);
MM=bsxfun(@rdivide,A,B);
V3=squeeze(nanmean(nanmean(MM(fsel,wind-50:wind+50,:),2),1));
 selC=cell_id==2&d_id==2;
pre_window= allC2(:,tsel, selC);pre_value2=nanmean(pre_window,2);
A=bsxfun(@minus, allC2(:,:, selC), pre_value2);B=bsxfun(@plus, allC2(:,:, selC), pre_value2);
MM=bsxfun(@rdivide,A,B);
V4=squeeze(nanmean(nanmean(MM(fsel,wind-50:wind+50,:),2),1));
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
%violinplot2(V2,[1.7 ],'ViolinColor', [ 0. 0. 0.9])
violinplotSTR(V3,[ 2.1],'ViolinColor', [ 0.9 0. 0.])
violinplotSTR(V4,[ 2.6],'ViolinColor', [ 0. 0. 0.9])
line([ 0.8 2.9], [ 0  0],'COlor', [0 0 0 ],'Linewidth',0.5)
 print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'spVM_beta_BARQUANT.pdf'])

 [h,p,ci,stats] =ttest2(V3,V4)
 
  [h,p,ci,stats] =ttest2(V2,V4)
 
  
  [h,p,ci,stats] =ttest(V1)
  % df=21, tst=5.5,  1.6941e-05
 [h,p,ci,stats] =ttest(V2)
  % df=3
   [h,p,ci,stats] =ttest(V3)
 
  %df=8; tst=7.8,      5.2610e-05
     [h,p,ci,stats] =ttest(V4)
  % df=15, tst=3.2,     0.0154
 
   [h,p,ci,stats] =ttest2(V1,V4)
%df=36, tst=2.7,     0.0162
   [h,p,ci,stats] =ttest2(V3,V4)
%df=23, tst=3.8,     3.6383e-04
   [h,p,ci,stats] =ttest2(V1,V3)
%df=29, tst=3.8,       0.2789
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for typC=1:4;
% if typC==1
% selC=d_id==1&cell_id==1
% elseif typC==2
%    selC=d_id==2&cell_id==1 
% elseif typC==3
%   selC=d_id==1&cell_id==2  
% elseif typC==4
%     selC=d_id==2&cell_id==2  
% end
% tt=-wind:wind;
% tsel=find(tt> -201 & tt <-100);
% pre_window= allC2(:,tsel,selC);
% pre_value2=nanmean(pre_window,2);
% figure('COlor','w','Position', [ 300 300 200 pheight],'Renderer', 'painters')
% A=bsxfun(@minus, allC2(:,:,selC), pre_value2);B=bsxfun(@plus, allC2(:,:,selC), pre_value2);
% imagesc(-wind:wind,freq2.freq,smooth2a(nanmean(bsxfun(@rdivide,A,B),3),1,20));colormap(jet)
% axis xy
% %hold on,,plot(-wind:wind,zscore(allC3).*10+20,'k','COlor',[0.4 0.4 0.4])
% set(gca,'CLim',[-0.1 0.1]) % -+10% modulation
% print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Vm_sptrig_pow_' num2str(typC) '.pdf'])
% end
