%cd('\\engnas.bu.edu\research\eng_research_handata\Mohammed_Abumuaileq\movement_tracking\For_Eric\')
clear all
pheight=160
savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\git_scripts\Fig4_mot\'
clear allCOHs allCOHsS V1 V2 V3 V4 V5 V6
pathn='\\engnas.bu.edu\research\eng_research_handata\Mohammed_Abumuaileq\movement_tracking\For_Eric\raw_motion_signals2\'
nk=0;
 tt=[];tt2=[]; 
foldernames{1}='31571';
foldernames{2}='4244';
foldernames{3}='23087';
foldernames{4}='23144';

for mm=1:length(foldernames)
    cd([pathn foldernames{mm}])
    for m=1:5;
        try
            if m==1;
                sesR=dir('*VID_RIGHT_motion_signals.mat');
                sesL=dir('*VID_LEFT_motion_signals.mat');
            elseif m==2
                sesR=dir('*VID1_RIGHT_motion_signals.mat');
                sesL=dir('*VID1_LEFT_motion_signals.mat');
            elseif m==3
                sesR=dir('*VID2_RIGHT_motion_signals.mat');
                sesL=dir('*VID2_LEFT_motion_signals.mat');
            elseif m==4
                sesR=dir('*VID3_RIGHT_motion_signals.mat');
                sesL=dir('*VID3_LEFT_motion_signals.mat');
            elseif m==5
                sesR=dir('*VID4_RIGHT_motion_signals.mat');
                sesL=dir('*VID4_LEFT_motion_signals.mat');
            end
            
      for jh=1:length(sesL)
          
            FS=120
            load(sesL(jh).name)
            V1=zscore(back_paw_d(1:1:end-0));
            V2=zscore(front_paw_d(1:1:end-0));
            
            load(sesR(jh).name)
           V3=zscore(back_paw_d(1:1:end-0));
            % V3=zscore(posB(2,nn3:nn3+length(back_paw_d)));%back_paw_d(1:1:end-0));
            V4=zscore(front_paw_d(1:1:end-0));
            V5=(interp1(1:length(ball_speed),ball_speed,1:1./6:length(ball_speed)));%
           V6 =fastsmooth(abs(hilbert(diff(V1))),500,1,1);%fastsmooth((V5),50,1,1);
           
           %%%%%%%%%%%%%%%
           
           
           %%%%%%%%%%%%%%%
           
           
            wind=1.5*FS;
            
            numWin=floor(length(V1)/wind);
            %%%
            lfp=[];nj=0;
            for  id=1:numWin
                selT=1 +(wind*(id-1)) :  (wind*id) ;
                if mean(V6(selT))>0.1
                    nj=nj+1;
                lfp.trial{nj}(1,:)= V1(  selT);
                lfp.trial{nj}(2,:)= V2(  selT);
                lfp.trial{nj}(3,:)= V3(  selT);
                lfp.trial{nj}(4,:)= V4(  selT);
                lfp.trial{nj}(5,:)= V5(  selT);
                lfp.time{nj}= (1:size(selT,2))./FS;
                end;end
            lfp.label= {'A1','A2','A3','A4','A5'}
            V5=V5(end:-1:1);
            lfp2=[];nj=0;
            for  id=1:numWin
                selT=1 +(wind*(id-1)) :  (wind*id) ;
                  if mean(V6(selT))>0.1
                    nj=nj+1;
                lfp2.trial{nj}(1,:)= V1(  selT);
                lfp2.trial{nj}(2,:)= V2(  selT);
                lfp2.trial{nj}(3,:)= V3(  selT);
                lfp2.trial{nj}(4,:)= V4(  selT);
                lfp2.trial{nj}(5,:)= V5(  selT);
                lfp2.time{nj}= (1:size(selT,2))./FS;end
            end
            lfp2.label= {'A1','A2','A3','A4','A5'}
            
            %%%%%%%%%%%%%%%%%%%%%

       cfg=[];
 cfg.bpfilter      = 'yes'
 cfg.bpfreq    =[1.5 5]
  cfg.bpfiltord    =2;
  cfg.hilbert='angle'
lfpF=ft_preprocessing(cfg,lfp) 

          
  for id2=1:length(lfpF.trial)         
     ph= lfpF.trial{id2}(1,:);                 
  frameSize=31;polyOrder=2;diffOrder=1;
inPhase=unwrap(ph,pi/4);
inFreq=savitzkyGolayFilt(inPhase,polyOrder,diffOrder,frameSize);
inFreq(end)=[];
tt=[tt,-1*120.*inFreq/(2*pi)]; 
     ph= lfpF.trial{id2}(5,:);                 
inPhase=unwrap(ph,pi/4);
inFreq=savitzkyGolayFilt(inPhase,polyOrder,diffOrder,frameSize);
inFreq(end)=[];
tt2=[tt2,-1*120.*inFreq/(2*pi)]; 
  end      
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            warning off
            cfg = []; %block_type == cfg.blk
            cfg.method ='mtmfft'; %'mvar';
            cfg.output ='fourier';
            cfg.taper='hanning';
            cfg.keeptapers ='yes';
            cfg.keeptrials ='yes';
            cfg.trials='all';cfg.tapsmofrq =6;%
            cfg.channel= ['all']; %chans=cfg.channel;
            cfg.foi= [0.5:0.5:20];
            freq2 = ft_freqanalysis(cfg, lfp);
            freq2S = ft_freqanalysis(cfg, lfp2);
            cfg=[];
            cfg.method='coh'
            coh=ft_connectivityanalysis(cfg,freq2)
            cohS=ft_connectivityanalysis(cfg,freq2S)
            
            
            %figure,plot(squeeze(coh.cohspctrm(1,5,:)));hold on,plot(squeeze(coh.cohspctrm(2,5,:)))
            nk=nk+1;
            allCOHs(1,:,nk)=squeeze(coh.cohspctrm(1,5,:));
             allCOHs(2,:,nk)=squeeze(coh.cohspctrm(2,5,:));
             allCOHs(3,:,nk)=squeeze(coh.cohspctrm(3,5,:));
               allCOHs(4,:,nk)=squeeze(coh.cohspctrm(4,5,:));
           % allCOHsS{nk}=cohS.cohspctrm;
               allCOHsS(1,:,nk)=squeeze(cohS.cohspctrm(1,5,:));
             allCOHsS(2,:,nk)=squeeze(cohS.cohspctrm(2,5,:));
             allCOHsS(3,:,nk)=squeeze(cohS.cohspctrm(3,5,:));
               allCOHsS(4,:,nk)=squeeze(cohS.cohspctrm(4,5,:));
%                   allCOHs(1,:,nk)=squeeze(coh.wppcspctrm(1,5,:));
%              allCOHs(2,:,nk)=squeeze(coh.wppcspctrm(2,5,:));
%              allCOHs(3,:,nk)=squeeze(coh.wppcspctrm(3,5,:));
%                allCOHs(4,:,nk)=squeeze(coh.wppcspctrm(4,5,:));
%            % allCOHsS{nk}=cohS.cohspctrm;
%                allCOHsS(1,:,nk)=squeeze(cohS.wppcspctrm(1,5,:));
%              allCOHsS(2,:,nk)=squeeze(cohS.wppcspctrm(2,5,:));
%              allCOHsS(3,:,nk)=squeeze(cohS.wppcspctrm(3,5,:));
%                allCOHsS(4,:,nk)=squeeze(cohS.wppcspctrm(4,5,:));
            
      end    
        end
    end
end


clear allst
for t=1:1000
    for fx=1:size(allCOHs,2)
   plvs= [squeeze(allCOHs(1,fx,:)) ;squeeze(allCOHsS(1,fx,:))];
    order_p=[ ones(1,size(allCOHs,3)),ones(1,size(allCOHs,3)).*2];
z=randperm(length(plvs));
plvs=plvs(z);
[h,p,ci,stats] =ttest(plvs(order_p==1), plvs(order_p==2));
allst(fx,t)=abs(stats.tstat);
    end
end
clear allstT
for fx=1:size(allCOHs,2)
[h,p,ci,stats] =ttest(squeeze(allCOHs(1,fx,:)), squeeze(allCOHsS(1,fx,:)));
allstT(fx)=abs(stats.tstat);
end
coh_SH_99=prctile(allst',95);
%%%%
sign_line= zeros(1,size(allCOHs,2)).*NaN;
 sign_line(find((allstT-coh_SH_99)>0))=1;

ff=freq2.freq,
M=nanmean(allCOHs(1,:,:),3);
Ms=nanstd(allCOHs(1,:,:),[],3)./sqrt(nk);
 figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters') 
plot(ff,M,'COlor',[0.9 0.5 .9])
fill_error_area2(ff,M,Ms,[0.5 0.5 0.5]);hold on,
M2=nanmean(allCOHsS(1,:,:),3);
Ms2=nanstd(allCOHsS(1,:,:),[],3)./sqrt(nk);
hold on,,plot(ff,M2,'k')
fill_error_area2(ff,M2,Ms2,[0.5 0.5 0.5]);hold on,
hold on,plot(ff,(sign_line.*max(allCOHs(:)))*0.75,'y.-','Linewidth',3,'COlor',[ 1 0.85 0])
plot(ff,M,'COlor',[0.4 0.8 0])
ylim([0.0 0.7]);xlim([0.5 18]);%set(gca,'Xscale','log')

print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'SFC_BALLmotion_Backright' '.pdf'])

