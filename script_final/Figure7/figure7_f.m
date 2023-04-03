%% Spike trigerred power plots 

clear all
%close all
%% Load files 
cell_id=[];Is_ratio3=[];d_id=[];

pT=[];pT2=[];  nn=0;nn1=0;
CHANNEL=1;
i_motion_threshold=300; % image motion threshold
 sw_t=3;
for gh=[ 1 2  ]
    
 if gh==1
        %cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Hua-an_Tseng/Data/SomArchon_Striatum/Linear_movement_PinnedBall/good_data_ALL/ChAT/')
        cd('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata2\Hua-an Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\ChAT\')
          %      cd('Z:\eng_research_handata2\Hua-an_Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\ChAT\')

        %%%%%%%%%%%% identification vector, 1=delta, 2=regular/other, 0= low SNR
      dyn_type=   [ 0,2,1,1,1,1,1,1,1,0,1,2,1,1,1,1,1,0,1,1,1,1,1,1,1,2,1,2,0,2,0,1]     
        %%%%%%%%%%%%%%
    else
        %cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Hua-an_Tseng/Data/SomArchon_Striatum/Linear_movement_PinnedBall/good_data_ALL/MSN/')
       cd('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata2\Hua-an Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\MSN\')
            %   cd('Z:\eng_research_handata2\Hua-an_Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\MSN\')

        %%%%%%%%%%%% identification vector, 1=delta, 2=regular/other, 0= low SNR
dyn_type= [2,2,2,2,0,2,2,0,2,2,0,1,1,2,0,1,2,1,1,2,0,1,1,1,1,2,2,2,2,0,2]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
ses= dir('*.mat')

for sesd=[1:length(ses)] 
    %sesd =5;
 
    load(ses(sesd).name)
try

   snr_thres= 0;%6.5;
% Only consider sessions with more than 3 trials 
if length(aligned.trace)>0% Actual - low movement

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
        
        
        d(3,:)=  locmotion_raw{in}.speed;
        %aligned.trace_nosp{in}(:,1:end);  %Vm
        
        
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
    if mean(spRES{in}.spike_snr{:}) >snr_thres
   for nh= 1:Nr
       mk=mk+1;
       lfp.trial{mk}(1:3,:) = zscore(d(1:3,(1  + ((windF)* (nh-1)):   windF * (nh))   ),[],2);
       
       lfp.trial{mk}(4:5,:) = (d(4:5,(1  + ((windF)* (nh-1)):   windF * (nh))   )); % locomotion
       lfp.trial{mk}(6,:) = (d(6,(1  + ((windF)* (nh-1)):   windF * (nh))   ));% image motion

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
    cfg.channel= [1 2]; 
    cfg.foi= [2:2:140];
    cfg.toi=lfp.time{1}(1:1:end) ;
   
    cfg.width =5;
   
    cfg.t_ftimwin =ones(1,length(cfg.foi))*0.15;
    freq2 = ft_freqanalysis(cfg, lfp);
    
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 
 
    wind=400;
    clear aM aM2  
    mm=0;   mm2=0; mm3=0;
    for id=1:size(  freq2.fourierspctrm,1)
        % Confirm 3 here -  Vm, no difference between M and MG ?
        M=abs(squeeze(freq2.fourierspctrm(id,CHANNEL,:,:)))      ;
        imag_mot= lfp.trial{id}(6,:); % GET image motion data
        M(:,  imag_mot>i_motion_threshold)=NaN;  % Segements above threshold are NaN
        
        
        MG=abs(squeeze(freq2.fourierspctrm(id,CHANNEL,:,1:end)))      ;
        
        MG(:,  imag_mot>i_motion_threshold)=NaN;% Segements above threshold are NaN
       
        
        s4= lfp.trial{id}(2,:)>0;%s4=find(s); % spike times
       vv= lfp.trial{id}(3,:);
              Fn = FS/2;FB=[ 2 3.5];
                                         [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                                         LFPg= ((filtfilt(B,A,  vv)));
                             
                                   Va=  angle(hilbert( LFPg));
                                %  VA= circ_dist(Va,pi);
               deltpeak  =zeros(1,length(Va));  
               timer=50;
                 for x1=1:length(Va)
                       timer=timer+1;
                     if Va(x1) >-0.01 &Va(x1) < 0.1 & timer>50
                    deltpeak(x1)= 1; timer=0; end        
                 end
        s=find(deltpeak);
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
        if s(id2) >wind  & s(id2) +wind < length(M) & any(s(id2) == HT)
            mm2=mm2+1;
            aM2(:,mm2)=   vv(s(id2)-wind:s(id2)+wind);
            %aM3_pre(:,:,mm2)= M(:,s(id2)-200:s(id2)-100);
        end
    end
        % Actual - low movement
    for id2=1:length(s)
        if s(id2) >wind  & s(id2) +wind < length(M) & any(s(id2) == HT)
            mm3=mm3+1;
            aM3(:,mm3)=   s4(s(id2)-wind:s(id2)+wind);
            %aM3_pre(:,:,mm2)= M(:,s(id2)-200:s(id2)-100);
        end
    end
    
end
    

if size(aM,3)>10% & size(aM2,32)>10
    nn=nn+1;
    % Mean across all windows
    allC(:,:,nn)=  nanmean(abs(aM),3);    
end
if size(aM,3)>10 %& size(aM2,2)>10
    nn1=nn1+1;
    allC2(:,nn1)=  nanmean((aM2),2);
   
       allC3(:,nn1)= fastsmooth(nanmean((aM3),2),30,1,1).*1000;
   
      cell_id=[cell_id,gh];
       d_id=[d_id,  dyn_type(sesd) ];
 end
    
end
end
end
end
%cd('./Shuffle mat files/Figure 3/Pre_spike_window_norm2_sub_by_union_no_abs/Median')
savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\git_scripts\fig4\'
pheight=160;

 selC=d_id>0&cell_id>0;%&cell_id==2;
tt=-wind:wind;
tsel=find(tt> -751 & tt <750);
pheight=160;
pre_window= allC(:,tsel, selC);
pre_value2=nanmean(pre_window,2);
figure('COlor','w','Position', [ 300 300 200  pheight],'Renderer', 'painters')
A=bsxfun(@minus, allC(:,:,selC), pre_value2);B=bsxfun(@plus, allC(:,:,selC), pre_value2);
%imagesc(-wind:wind,freq2.freq,smooth2a(nanmean(A,3)./nanmean(B,3),1,30));colormap(jet)
imagesc(-wind:wind,freq2.freq,smooth2a(nanmean(A,3)./nanmean(B,3),1,50));colormap(jet)
%imagesc(-wind:wind,freq2.freq,zscore(smooth2a(nanmean(A,3),1,30),[],2));colormap(jet)
hold on,,plot(-wind:wind,(zscore(nanmean(allC2,2)).*8)+20,'k')
hold on,,plot(-wind:wind,fastsmooth(zscore(nanmean(allC3,2)).*8,10,1,1)+20,'m')

axis xy
set(gca,'CLim',[-0.025 0.025]) % -+10% modulation
ylim([2 100])
% saveas(gcf, 'LFP_power_spec_low_motion.fig')
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Spoow_LFP_CHAT_delta_h' num2str(sw_t) '.pdf'])

%  selC=d_id==1&cell_id==1;
% tt=-wind:wind;
% tsel=find(tt> -201 & tt <-100);
% pheight=160;
% pre_window= allC(:,tsel, selC);
% pre_value2=nanmean(pre_window,2);
% figure('COlor','w','Position', [ 300 300 200  pheight],'Renderer', 'painters')
% A=bsxfun(@minus, allC(:,:,selC), pre_value2);B=bsxfun(@plus, allC(:,:,selC), pre_value2);
% imagesc(-wind:wind,freq2.freq,smooth2a(nanmedian(A,3)./nanmedian(B,3),1,30));colormap(jet)
% axis xy
% set(gca,'CLim',[-0.05 0.05]) % -+10% modulation
% %print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Spoow_LFP_CHAT_delta_l' num2str(sw_t) '.pdf'])
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%
% tt=-wind:wind;
% tsel=find(tt> -201 & tt <-100);
%  selC=d_id==2&cell_id==2;
% tt=-wind:wind;
% tsel=find(tt> -201 & tt <-100);
% pheight=160;
% pre_window= allC(:,tsel, selC);
% pre_value2=nanmedian(pre_window,2);
% A=bsxfun(@minus, allC(:,:,selC), pre_value2);B=bsxfun(@plus, allC(:,:,selC), pre_value2);
% CC1=smooth2a(nanmedian(A,3)./nanmedian(B,3),1,30);
% pre_window= allC2(:,tsel, selC);
% pre_value2=nanmedian(pre_window,2);
% A=bsxfun(@minus, allC2(:,:,selC), pre_value2);B=bsxfun(@plus, allC2(:,:,selC), pre_value2);
% CC2=smooth2a(nanmedian(A,3)./nanmedian(B,3),1,30);
% figure('COlor','w','Position', [ 300 300 200  pheight],'Renderer', 'painters')
% imagesc(-wind:wind,freq2.freq,CC1-CC2);colormap(jet)
% axis xy
% set(gca,'CLim',[-0.07 0.07]) % -+10% modulation
% xlim([-200 200])
% ylim([1 100])
% %print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Sp_LFP_MSN_reg' num2str(sw_t) '.pdf'])
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tt=-wind:wind; fsel=find(freq2.freq>=70  & freq2.freq<=140);
% tsel=find(tt> -201 & tt <-100);
%  selC=d_id==1&cell_id==1;
% pre_window= allC2(fsel,tsel, selC);pre_value2=nanmean(pre_window,2);
% A=squeeze(nanmean(bsxfun(@minus, allC2(fsel,wind+1,selC), pre_value2),1));B=squeeze(nanmean(bsxfun(@plus, allC2(fsel,wind+1,selC), pre_value2),1))
% X1=A./B;pre_window= allC(fsel,tsel, selC);pre_value2=nanmean(pre_window,2);
% A=squeeze(nanmean(bsxfun(@minus, allC(fsel,wind+1,selC), pre_value2),1));B=squeeze(nanmean(bsxfun(@plus, allC(fsel,wind+1,selC), pre_value2),1))
% X2=A./B;
% V1=X2-X1;
%  selC=d_id==2&cell_id==1;
% pre_window= allC2(fsel,tsel, selC);pre_value2=nanmean(pre_window,2);
% A=squeeze(nanmean(bsxfun(@minus, allC2(fsel,wind+1,selC), pre_value2),1));B=squeeze(nanmean(bsxfun(@plus, allC2(fsel,wind+1,selC), pre_value2),1))
% X1=A./B;pre_window= allC(fsel,tsel, selC);pre_value2=nanmean(pre_window,2);
% A=squeeze(nanmean(bsxfun(@minus, allC(fsel,wind+1,selC), pre_value2),1));B=squeeze(nanmean(bsxfun(@plus, allC(fsel,wind+1,selC), pre_value2),1))
% X2=A./B;
% V2=X2-X1;
%  selC=d_id==1&cell_id==2;
% pre_window= allC2(fsel,tsel, selC);pre_value2=nanmean(pre_window,2);
% A=squeeze(nanmean(bsxfun(@minus, allC2(fsel,wind+1,selC), pre_value2),1));B=squeeze(nanmean(bsxfun(@plus, allC2(fsel,wind+1,selC), pre_value2),1))
% X1=A./B;pre_window= allC(fsel,tsel, selC);pre_value2=nanmean(pre_window,2);
% A=squeeze(nanmean(bsxfun(@minus, allC(fsel,wind+1,selC), pre_value2),1));B=squeeze(nanmean(bsxfun(@plus, allC(fsel,wind+1,selC), pre_value2),1))
% X2=A./B;
% V3=X2-X1;
%  selC=d_id==2&cell_id==2;
% pre_window= allC2(fsel,tsel, selC);pre_value2=nanmean(pre_window,2);
% A=squeeze(nanmean(bsxfun(@minus, allC2(fsel,wind+1,selC), pre_value2),1));B=squeeze(nanmean(bsxfun(@plus, allC2(fsel,wind+1,selC), pre_value2),1))
% X1=A./B;pre_window= allC(fsel,tsel, selC);pre_value2=nanmean(pre_window,2);
% A=squeeze(nanmean(bsxfun(@minus, allC(fsel,wind+1,selC), pre_value2),1));B=squeeze(nanmean(bsxfun(@plus, allC(fsel,wind+1,selC), pre_value2),1))
% X2=A./B;
% V4=X2-X1;
% figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters') 
% b1=bar(1,nanmean(V1),'Facecolor',[ 0.8 0 0]);hold on,
% set(b1,'FaceAlpha',0.4)
% %  b1=bar(2,nanmean(V2),'Facecolor',[ 0 0 0.8]);hold on,
% %set(b1,'FaceAlpha',0.4)
%   b1=bar(3,nanmean(V3),'Facecolor',[ 0.8 0 0]);hold on,
% set(b1,'FaceAlpha',0.4)
%   b1=bar(4,nanmean(V4),'Facecolor',[ 0 0 0.8]);hold on,
% set(b1,'FaceAlpha',0.4)
% errorbar([1   ],nanmean(V1,1), nanstd(V1)./sqrt(size(V1,1)),'.k');
% %errorbar([2   ],nanmean(V2,1), nanstd(V2)./sqrt(size(V2,1)),'.k');
% errorbar([3   ],nanmean(V3,1), nanstd(V3)./sqrt(size(V3,1)),'.k');
% errorbar([4   ],nanmean(V4,1), nanstd(V4)./sqrt(size(V4,1)),'.k');
% print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'high_F_BARquant_sppow' num2str(sw_t) '.pdf'])
% ylim([-0.06 0.06])
% 
% 
% 
%  [h,p,ci,stats] = ttest(V1)
%  %df=20, 2.7647, 0.0120
%  [h,p,ci,stats] = ttest(V3)
%  %df=8, -2.2491, 0.054
%   [h,p,ci,stats] = ttest(V4)
%  %df=15, -0.6673, 0.5147
%  
% tt=-wind:wind; fsel=find(freq2.freq>=20  & freq2.freq<=40);
% tsel=find(tt> -201 & tt <-100);
%  selC=d_id==1&cell_id==1;
% pre_window= allC2(fsel,tsel, selC);pre_value2=nanmean(pre_window,2);
% A=squeeze(nanmean(bsxfun(@minus, allC2(fsel,wind+1,selC), pre_value2),1));B=squeeze(nanmean(bsxfun(@plus, allC2(fsel,wind+1,selC), pre_value2),1))
% X1=A./B;pre_window= allC(fsel,tsel, selC);pre_value2=nanmean(pre_window,2);
% A=squeeze(nanmean(bsxfun(@minus, allC(fsel,wind+1,selC), pre_value2),1));B=squeeze(nanmean(bsxfun(@plus, allC(fsel,wind+1,selC), pre_value2),1))
% X2=A./B;
% V1=X2-X1;
%  selC=d_id==2&cell_id==1;
% pre_window= allC2(fsel,tsel, selC);pre_value2=nanmean(pre_window,2);
% A=squeeze(nanmean(bsxfun(@minus, allC2(fsel,wind+1,selC), pre_value2),1));B=squeeze(nanmean(bsxfun(@plus, allC2(fsel,wind+1,selC), pre_value2),1))
% X1=A./B;pre_window= allC(fsel,tsel, selC);pre_value2=nanmean(pre_window,2);
% A=squeeze(nanmean(bsxfun(@minus, allC(fsel,wind+1,selC), pre_value2),1));B=squeeze(nanmean(bsxfun(@plus, allC(fsel,wind+1,selC), pre_value2),1))
% X2=A./B;
% V2=X2-X1;
%  selC=d_id==1&cell_id==2;
% pre_window= allC2(fsel,tsel, selC);pre_value2=nanmean(pre_window,2);
% A=squeeze(nanmean(bsxfun(@minus, allC2(fsel,wind+1,selC), pre_value2),1));B=squeeze(nanmean(bsxfun(@plus, allC2(fsel,wind+1,selC), pre_value2),1))
% X1=A./B;pre_window= allC(fsel,tsel, selC);pre_value2=nanmean(pre_window,2);
% A=squeeze(nanmean(bsxfun(@minus, allC(fsel,wind+1,selC), pre_value2),1));B=squeeze(nanmean(bsxfun(@plus, allC(fsel,wind+1,selC), pre_value2),1))
% X2=A./B;
% V3=X2-X1;
%  selC=d_id==2&cell_id==2;
% pre_window= allC2(fsel,tsel, selC);pre_value2=nanmean(pre_window,2);
% A=squeeze(nanmean(bsxfun(@minus, allC2(fsel,wind+1,selC), pre_value2),1));B=squeeze(nanmean(bsxfun(@plus, allC2(fsel,wind+1,selC), pre_value2),1))
% X1=A./B;pre_window= allC(fsel,tsel, selC);pre_value2=nanmean(pre_window,2);
% A=squeeze(nanmean(bsxfun(@minus, allC(fsel,wind+1,selC), pre_value2),1));B=squeeze(nanmean(bsxfun(@plus, allC(fsel,wind+1,selC), pre_value2),1))
% X2=A./B;
% V4=X2-X1;
% figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters') 
% b1=bar(1,nanmean(V1),'Facecolor',[ 0.8 0 0]);hold on,
% set(b1,'FaceAlpha',0.4)
% %  b1=bar(2,nanmean(V2),'Facecolor',[ 0 0 0.8]);hold on,
% %set(b1,'FaceAlpha',0.4)
%   b1=bar(3,nanmean(V3),'Facecolor',[ 0.8 0 0]);hold on,
% set(b1,'FaceAlpha',0.4)
%   b1=bar(4,nanmean(V4),'Facecolor',[ 0 0 0.8]);hold on,
% set(b1,'FaceAlpha',0.4)
% errorbar([1   ],nanmean(V1,1), nanstd(V1)./sqrt(size(V1,1)),'.k');
% %errorbar([2   ],nanmean(V2,1), nanstd(V2)./sqrt(size(V2,1)),'.k');
% errorbar([3   ],nanmean(V3,1), nanstd(V3)./sqrt(size(V3,1)),'.k');
% errorbar([4   ],nanmean(V4,1), nanstd(V4)./sqrt(size(V4,1)),'.k');
% print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'l_F_BARquant_sppow' num2str(sw_t) '.pdf'])
% ylim([-0.06 0.06])
%  [h,p,ci,stats] = ttest(V1)
%  %df=20, 0.0524, 0.9588
%  [h,p,ci,stats] = ttest(V3)
%  %df=8, --3.4147,0.0092
%   [h,p,ci,stats] = ttest(V4)
%  %df=15, -0.4515,  0.6581
%  
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  
%  %%%%%%%%%%%%%%%%%%%%%%%
%  for typC=1:4;
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
% pheight=160;
% pre_window= allC(:,tsel, selC);
% pre_value2=nanmean(pre_window,2);
% figure('COlor','w','Position', [ 300 300 200  pheight],'Renderer', 'painters')
% A=bsxfun(@minus, allC(:,:,selC), pre_value2);B=bsxfun(@plus, allC(:,:,selC), pre_value2);
% imagesc(-wind:wind,freq2.freq,smooth2a(nanmedian(A,3)./nanmedian(B,3),1,30));colormap(jet)
% axis xy
% set(gca,'CLim',[-0.08 0.08]) % -+10% modulation
% ylim([1 100])
% print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Sp_LFP_suB' num2str(typC) '.pdf'])
%  end
%  
%  
%   
%  %%%%%%%%%%%%%%%%%%%%%%%
%  for typC=1:4;
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
% pheight=160;
% pre_window= allC2(:,tsel, selC);
% pre_value2=nanmean(pre_window,2);
% figure('COlor','w','Position', [ 300 300 200  pheight],'Renderer', 'painters')
% A=bsxfun(@minus, allC2(:,:,selC), pre_value2);B=bsxfun(@plus, allC2(:,:,selC), pre_value2);
% imagesc(-wind:wind,freq2.freq,smooth2a(nanmedian(A,3)./nanmedian(B,3),1,30));colormap(jet)
% axis xy
% set(gca,'CLim',[-0.08 0.08]) % -+10% modulation
% ylim([1 100])
% print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Sp_LFP_HIGHMOT' num2str(typC) '.pdf'])
%  end