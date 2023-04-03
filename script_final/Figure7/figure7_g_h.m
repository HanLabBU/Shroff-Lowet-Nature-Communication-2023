%% Spike trigerred power plots 

clear all
%close all
%% Load files 
cell_id=[];Is_ratio3=[];d_id=[];

pT=[];pT2=[];  nn=0;nn1=0;
CHANNEL=1; DD=[]; EE=[];
i_motion_threshold=300; % image motion threshold
 sw_t=3;
for gh=[ 1 2  ]
    
 if gh==1
        %cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Hua-an_Tseng/Data/SomArchon_Striatum/Linear_movement_PinnedBall/good_data_ALL/ChAT/')
        cd('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata2\Hua-an Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\ChAT\')
          %      cd('Z:\eng_research_handata2\Hua-an_Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\ChAT\')
      dyn_type=   [ 0,2,1,1,1,1,1,1,1,0,1,2,1,1,1,1,1,0,1,1,1,1,1,1,1,2,1,2,0,2,0,1]     

        if sw_t==1 % delta-identified
            chs=[find(dyn_type==1)];
        elseif sw_t==2 % regular-identified 
            chs=  [ find(dyn_type==2)];
        elseif sw_t==3   % all
            chs=[find(dyn_type>-1)];
            
        end
        %%%%%%%%%%%% identification vector, 1=delta, 2=regular/other, 0= low SNR
      dyn_type=   [ 0,2,1,1,1,1,1,1,1,0,1,2,1,1,1,1,1,0,1,1,1,1,1,1,1,2,1,2,0,2,0,1]     
        %%%%%%%%%%%%%%
    else
        %cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Hua-an_Tseng/Data/SomArchon_Striatum/Linear_movement_PinnedBall/good_data_ALL/MSN/')
       cd('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata2\Hua-an Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\MSN\')
            %   cd('Z:\eng_research_handata2\Hua-an_Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\MSN\')
dyn_type= [2,2,2,2,0,2,2,0,2,2,0,1,1,2,0,1,2,1,1,2,0,1,1,1,1,2,2,2,2,0,2]

        if sw_t==1      % delta-identified
            chs=[find(dyn_type==1)];
        elseif sw_t==2     % regular-identified or theta bursting
            chs=  [  find(dyn_type==2) ]
        elseif sw_t==3      % all
            chs=[find(dyn_type>-1)];
            
        end
        %%%%%%%%%%%% identification vector, 1=delta, 2=regular/other, 0= low SNR
dyn_type= [2,2,2,2,0,2,2,0,2,2,0,1,1,2,0,1,2,1,1,2,0,1,1,1,1,2,2,2,2,0,2]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
ses= dir('*.mat')

for sesd=[chs]  %%%1:length(ses) %[1:25 27:32]%length(ses);%[1:15 18:31] %:length(ses)  %22 exclude 16 17
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
if length(aligned.trace)>3% Actual - low movement

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
    cfg.foi= [2:4:100];
    cfg.toi=lfp.time{1}(1:1:end) ;
   
    cfg.width =5;
   
    cfg.t_ftimwin =ones(1,length(cfg.foi))*0.15;
    freq2 = ft_freqanalysis(cfg, lfp);
    
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 
 rang=[];
    wind=400;
    clear aM aM2  
    mm=0;   mm2=0; mm3=0;AA=[];BB=[];CC=[];
    for id=1:size(  freq2.fourierspctrm,1)
        % Confirm 3 here -  Vm, no difference between M and MG ?
        M=abs(squeeze(freq2.fourierspctrm(id,CHANNEL,:,:)))      ;
        imag_mot= lfp.trial{id}(6,:); % GET image motion data
        M(:,  imag_mot>i_motion_threshold)=NaN;  % Segements above threshold are NaN
        
        
        MG=abs(squeeze(freq2.fourierspctrm(id,CHANNEL,:,1:end)))      ;
        
        MG(:,  imag_mot>i_motion_threshold)=NaN;% Segements above threshold are NaN
       
        
        s4= lfp.trial{id}(2,:)>0;%s4=find(s); % spike times
       vv= lfp.trial{id}(3,:);
              Fn = FS/2;FB=[ 1 3];
                                         [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                                         LFPg= ((filtfilt(B,A,  vv)));
                             
                                   Va=  angle(hilbert( LFPg));
                                       vv= lfp.trial{id}(1,:);
              Fn = FS/2;FB=[ 70 140];
                                         [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                                         LFPg= ((filtfilt(B,A,  vv)));
                                   GP=  abs(hilbert( LFPg));
                                        Fn = FS/2;FB=[ 10 40];
                                         [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                                         LFPg= ((filtfilt(B,A,  vv)));
                                   BP=  abs(hilbert( LFPg));
                                %  VA= circ_dist(Va,pi);
                                
                                 rang=[rang,Va(s4)];
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
    AA=[AA,Va(HT)];BB=[BB,BP(HT)];CC=[CC,GP(HT)];

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
    
        DD=[DD,circ_corrcl(AA,BB)]; EE=[EE,circ_corrcl(AA,CC)];
    nn=nn+1;
    % Mean across all windows
    tphase=AA;clear delB delG
    phase_sel=linspace(-pi,pi,63);
    for ph=1:length(phase_sel)-1

         self= find(tphase >= phase_sel(ph)  &  tphase <= phase_sel(ph+1) );
         delB(ph)=nanmean(BB(self));     delG(ph)=nanmean(CC(self));
    end  
   % figure,plot( phase_sel(1:end-1),fastsmooth(delB,5,1,1))
   [z z1]=findpeaks(fastsmooth(BB,100,1,1));
     [z z2]=findpeaks(fastsmooth(CC,100,1,1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    allC(:,:,nn)=  nanmean(abs(aM),3);  
    allPLV(nn)= (nanmean(exp(1i.*rang)));
      allAngB(nn)=circ_mean(AA(z1)');%  angle(sum(BB./sum(BB).*exp(1i.*AA)));
         allAngBP(nn)=  circ_r(AA(z1)');%abs(sum(BB./sum(BB).*exp(1i.*AA)));
               allAngG(nn)=  circ_mean(AA(z2)');;%angle(sum(CC./sum(CC).*exp(1i.*AA)));
         allAngGP(nn)= circ_r(AA(z2)'); %abs(sum(CC./sum(CC).*exp(1i.*AA)));
end
if size(aM,3)>10 %& size(aM2,2)>10
    nn1=nn1+1;
    allC2(:,nn1)=  nanmean((aM2),2);
   
       allC3(:,nn1)= fastsmooth(nanmean((aM3),2),50,1,1).*1000;
    allC3(:,nn1)= allC3(:,nn1)./mean(allC3(:,nn1));
    % allC3(:,nn1)= allC3(:,nn1)./max(allC3(:,nn1));
  
      cell_id=[cell_id,gh];
       d_id=[d_id,  dyn_type(sesd) ];
 end
    
end
end
end
end
%cd('./Shuffle mat files/Figure 3/Pre_spike_window_norm2_sub_by_union_no_abs/Median')
savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\git_scripts\fig4_mot\'
pheight=160;

 selC=d_id>0&cell_id>0;%&cell_id==2;
tt=-wind:wind;
tsel=find(tt> -751 & tt <750);
pheight=160;
pre_window= allC(:,tsel, selC);
pre_value2=nanmean(pre_window,2);



figure('COlor','w','Position', [ 300 300 200  pheight],'Renderer', 'painters')
polarscatter((allAngB(cell_id>0)),abs(allAngBP(cell_id>0)),'k','filled','MarkerFacealpha',0.07);hold on,
%hold on,polarhistogram(allAngB(cell_id>0),12,'FaceColor',[0.6 0.6 0.2],'Normalization','probability','Facealpha',0.4)
polarscatter((allAngB(cell_id>0)),abs(allAngBP(cell_id>0)),'k','filled','MarkerFacealpha',0.5);hold on,
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'polarscatt_Beta' num2str(sw_t) '.pdf'])

[pval m] = circ_otest(allAngB(cell_id>0),1)

%Omnibus or Hodges-Ajne test for non-uniformity of circular data.
% n=59, p=0.0116

figure('COlor','w','Position', [ 300 300 200  pheight],'Renderer', 'painters')
polarscatter((allAngG(cell_id>0)),abs(allAngGP(cell_id>0)),'k','filled','MarkerFacealpha',0.07);hold on,
%hold on,polarhistogram(allAngG(cell_id>0),12,'FaceColor',[0.3 0.6 0.6],'Normalization','probability','Facealpha',0.4)
polarscatter((allAngG(cell_id>0)),abs(allAngGP(cell_id>0)),'k','filled','MarkerFacealpha',0.5);hold on,
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'polarscatt_HGAM' num2str(sw_t) '.pdf'])

[pval m] = circ_otest(allAngG(cell_id>-3), 1)
%polarscatter(angle(allPLV(cell_id==2)),abs(allPLV(cell_id==2)),'.b')
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Firing_motion_deltaPhase' num2str(sw_t) '.pdf'])

%Omnibus or Hodges-Ajne test for non-uniformity of circular data.
% n=59, p=0.0048