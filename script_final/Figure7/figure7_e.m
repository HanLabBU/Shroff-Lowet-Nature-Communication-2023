%cd('\\engnas.bu.edu\research\eng_research_handata\Mohammed_Abumuaileq\movement_tracking\For_Eric\')

pheight=160
savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\git_scripts\Fig4_mot\'
clear allCOHs allCOHsS V1 V2 V3 V4 V5 V6
nk=0;
 tt=[];tt2=[]; 
foldernames{1}='31571';
foldernames{2}='4244';
foldernames{3}='23087';
foldernames{4}='23144';
foldernames{5}='31613';
 clear aM aM2   aM3
    mm=0;   mm2=0; mm3=0;
for mmv=[ 1  3  5]% sessions with paw tracking and LFP
    
    pathn='\\engnas.bu.edu\research\eng_research_handata\Mohammed_Abumuaileq\movement_tracking\For_Eric\raw_motion_signals3\'

    
    cd([pathn foldernames{mmv}])
    sesloc{1}=[2 3 4];
        sesloc{2}=[];
           sesloc{3}=[ 2 3 4 5];
                  sesloc{5}=[2 3 ];
          % 23087_20220520-2   23087_20220520-1
    for m= sesloc{mmv}
        try
            if mmv==1
            if m==1;
            elseif m==2
                sesR=dir('*504_VID1_RIGHT_motion_signals.mat'); sesL=dir('*504_VID1_LEFT_motion_signals.mat');
               fx= '31571_20220504-1.plx';  %%%%
            elseif m==3
                sesR=dir('*504_VID2_RIGHT_motion_signals.mat'); sesL=dir('*504_VID2_LEFT_motion_signals.mat');
                  fx= '31571_20220504-2.plx'; %%%
                elseif m==4
                sesR=dir('*608_VID1_RIGHT_motion_signals.mat'); sesL=dir('*608_VID1_LEFT_motion_signals.mat');
                  fx= '31571_20220608-1.plx'; %%%    
                  
            end
            elseif mmv==3
                 if m==1;
            elseif m==2
                sesR=dir('*520_VID1_RIGHT_motion_signals.mat'); sesL=dir('*520_VID1_LEFT_motion_signals.mat');
               fx= '23087_20220520-1.plx';  %%%%
            elseif m==3
                sesR=dir('*520_VID2_RIGHT_motion_signals.mat'); sesL=dir('*520_VID2_LEFT_motion_signals.mat');
                  fx= '23087_20220520-2.plx'; %%%
                    elseif m==4
                sesR=dir('*608_VID1_RIGHT_motion_signals.mat'); sesL=dir('*608_VID1_LEFT_motion_signals.mat');
                  fx= '23087_20220608-1.plx'; %%%
                         elseif m==5
                sesR=dir('*608_VID2_RIGHT_motion_signals.mat'); sesL=dir('*608_VID2_LEFT_motion_signals.mat');
                  fx= '23087_20220608-2.plx'; %%%
                  
                  
                 end 
            
                  elseif mmv==5
                 if m==1;
            elseif m==2
                sesR=dir('*520_VID1_RIGHT_motion_signals.mat'); sesL=dir('*520_VID1_LEFT_motion_signals.mat');
               fx= '31613_20220520-1.plx';  %%%%
               
            elseif m==3
           sesR=dir('*608_VID3_RIGHT_motion_signals.mat'); sesL=dir('*608_VID3_LEFT_motion_signals.mat');
               fx= '31613_20220608-3.plx';  %%%%
            end 
            
            end
            
      for jh=1:length(sesL)
          
            FS=1000
            load(sesL(jh).name)
            V1=zscore(back_paw_d(1:1:end-0));
            V2=zscore(front_paw_d(1:1:end-0));
            %nn3=nn2;
            load(sesR(jh).name)
            V3=zscore(back_paw_d(1:1:end-0));
            % V3=zscore(posB(2,nn3:nn3+length(back_paw_d)));%back_paw_d(1:1:end-0));
            V4=zscore(front_paw_d(1:1:end-0));
            V5=(interp1(1:length(ball_speed),ball_speed,1:1./6:length(ball_speed)));%
            V6 =fastsmooth(abs(hilbert(diff(V1))),500,1,1);%fastsmooth((V5),50,1,1);
            
           
             V1a=(interp1(1:length(V1),V1,1:1./8.3333:length(V1)));%
                V2a=(interp1(1:length(V2),V2,1:1./8.3333:length(V2)));%
            V3a=(interp1(1:length(V3),V3,1:1./8.3333:length(V3)));%
            V4a=(interp1(1:length(V4),V4,1:1./8.3333:length(V4)));%
           V5a=(interp1(1:length(V5),V5,1:1./8.3333:length(V5)));%
            V6a=(interp1(1:length(V6),V6,1:1./8.3333:length(V6)));%
         
           %%%%%%%%%%%%%%%
           
    %     clear all
    pp=['\\engnas.bu.edu\research\eng_research_handata\Mohammed_Abumuaileq\movement_tracking\For_Eric\LFP\' foldernames{mmv} '\'];
OpenedFileName=[ pp fx];

[adfreq, n, ts, fn, ad] = plx_ad(OpenedFileName, 1);
% get some counts
[tscounts, wfcounts, evcounts, slowcounts] = plx_info(OpenedFileName,1);
% and finally the events
[u,nevchannels] = size( evcounts );
if ( nevchannels > 0 )
    % need the event chanmap to make any sense of these
    [u,evchans] = plx_event_chanmap(OpenedFileName);
    for iev = 1:3%nevchannels
        if ( evcounts(iev) > 0 )
            evch = evchans(iev);
            [nevs{iev}, tsevs{iev}, svdummy] = plx_event_ts(OpenedFileName, evch);
        end;end;end
[nev,evnames] = plx_event_names(OpenedFileName);
ad1=ad(1:40:end);
%figure,plot(ad1(round(tsevs{3}(1).*1000):1:end))  
 ad1=     ad1(round(tsevs{3}(1).*1000):1:end) ; 
 ad1=ad1(1:length(V1a));
        Fn = 1000/2;FB=[ 1 3];
               [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
            LFPg= ((filtfilt(B,A,  ad1)));    
       %  figure('Color','w'),
       %  plot(zscore(V1a),'k'); hold on,plot(zscore(LFPg(1:1:end))+1,'Color',[0.7 0 0.3])
         
          if 0
         [s,w,t]=    spectrogram(ad1,500,480,[1:100],1000); 
              figure
              imagesc(abs(s))
              axis xy;colormap(jet)
              
          end
           
           
           %%%%%%%%%%%%%%%
           
           
            wind=1*FS;
            
            numWin=floor(length(V1a)/wind);
            %%%
            lfp=[];nj=0;
            for  id=1:numWin
                selT=1 +(wind*(id-1)) :  (wind*id) ;
                if mean(V6a(selT))>0.2
                    nj=nj+1;
                lfp.trial{nj}(1,:)= V1a(  selT);
                lfp.trial{nj}(2,:)= V2a(  selT);
                lfp.trial{nj}(3,:)= V3a(  selT);
                lfp.trial{nj}(4,:)= ad1(  selT);
                lfp.trial{nj}(5,:)= V5a(  selT);
                lfp.time{nj}= (1:size(selT,2))./FS;
                end;end
            lfp.label= {'A1','A2','A3','A4','A5'}
        ad1f=ad1(end:-1:1);
    %   V5a=V5a(end:-1:1); 
                    lfp2=[];nj=0;
            for  id=1:numWin
                selT=1 +(wind*(id-1)) :  (wind*id) ;
                  if  mean(V6a(selT))>0.2
                    nj=nj+1;
                lfp2.trial{nj}(1,:)= V1a(  selT);
                lfp2.trial{nj}(2,:)= V2a(  selT);
                lfp2.trial{nj}(3,:)= V3a(  selT);
                lfp2.trial{nj}(4,:)= ad1f(  selT);
                lfp2.trial{nj}(5,:)= V5a(  selT);
                lfp2.time{nj}= (1:size(selT,2))./FS;end
            end
            lfp2.label= {'A1','A2','A3','A4','A5'}
         

%             
            %%%%%%%%%%%%%%%%%%%%%
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
    allCOHs(1,:,nk)=squeeze(coh.cohspctrm(1,4,:));
       allCOHs(2,:,nk)=squeeze(coh.cohspctrm(5,4,:)); 
    allCOHs(3,:,nk)=squeeze(cohS.cohspctrm(1,4,:));
              
      end    
        end
    end
end

clear allstT2
for fx=1:size(allCOHs,2)
[h,p,ci,stats] =ttest(squeeze(allCOHs(2,fx,:)), squeeze(allCOHs(3,fx,:)));
allstT2(fx)=abs(stats.tstat);
end


savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\git_scripts\Fig4_mot\'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pheight=160;

tt=freq2.freq,
M=nanmean(allCOHs(1,:,:),3);
Ms=nanstd(allCOHs(1,:,:),[],3)./sqrt(nk);
 figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters') 
plot(tt,M,'b','COlor',[0.0 0.3 0])
fill_error_area2(tt,M,Ms,[0.5 0.5 0.5]);hold on,
plot(tt,M,'b','COlor',[0.8 0.2 0],'Linewidth',1.5)
M=nanmean(allCOHs(2,:,:),3);
Ms=nanstd(allCOHs(2,:,:),[],3)./sqrt(nk);
hold on,,plot(tt,M,'k')
fill_error_area2(tt,M,Ms,[0.5 0.5 0.5]);hold on,
plot(tt,M,'b','COlor',[0.2 0.8 0.3],'Linewidth',1.5)
M=nanmean(allCOHs(3,:,:),3);
Ms=nanstd(allCOHs(3,:,:),[],3)./sqrt(nk);
hold on,,plot(tt,M,'k','Linewidth',1)
fill_error_area2(tt,M,Ms,[0.5 0.5 0.5]);hold on,
;xlim([1 17]);
%set(gca,'Xscale','log')
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'SFC_LFP_BALLmotion_Backright' '.pdf'])
