%addpath(genpath('\\engnas.bu.edu\Research\eng_research_handata\EricLowet\Scripts\'))


clear all
pT=[];pT2=[];  nn=0;
AC=[];AC2=[];
AC=0;AC2=0;
N1=0;N2=0;
i_motion_threshold=300; % image motion threshold
 sw_t=1;
for gh=[ 1 2]
   if gh==1
        %cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Hua-an_Tseng/Data/SomArchon_Striatum/Linear_movement_PinnedBall/good_data_ALL/ChAT/')
        cd('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata2\Hua-an Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\ChAT\')
       %         cd('Z:\eng_research_handata2\Hua-an_Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\ChAT\')
 dyn_type= [ 0 2 1 1 1 1 1 1 1 1    1 2 1 1 2 1 1 0 1 1   1 0 1 1 0 2 1 2 0 0  0 0];
 dyn_type= [ 0,2,1,1,1,1,1,1,1,0,1,2,1,1,1,1,1,0,1,1,1,1,1,1,1,2,1,2,0,2,0,1] 
 
 if sw_t==1 % delta-identified
            chs=[find(dyn_type==1)]
        elseif sw_t==2 % regular-identified 
            chs=  [ find(dyn_type==2) ];
        elseif sw_t==3   % all
            chs=[2 12 15  26 28,3,4,5,6,7,8,9,10,11,13,14,16,17,19,20,21,23,24,27,29];
            
        end
        %%%%%%%%%%%% identification vector, 1=delta, 2=regular/other, 0= low SNR
        dyn_type= [ 0 2 1 1 1 1 1 1 1 1    1 2 1 1 2 1 1 0 1 1   1 0 1 1 0 2 1 2 0 0  0 0];
       
        %%%%%%%%%%%%%%
    else
        %cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Hua-an_Tseng/Data/SomArchon_Striatum/Linear_movement_PinnedBall/good_data_ALL/MSN/')
        cd('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata2\Hua-an Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\MSN\')
       %        cd('Z:\eng_research_handata2\Hua-an_Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\MSN\')
dyn_type= [2,2,2,2,0,2,2,0,2,2,0,1,1,2,0,1,2,1,1,2,0,1,1,1,1,2,2,2,2,0,2, 0 0 0 0]

        if sw_t==1      % delta-identified
           chs=[find(dyn_type==1)];
        elseif sw_t==2     % regular-identified or theta bursting
            chs=  [ find(dyn_type==2) ]
        elseif sw_t==3      % all
            chs=[[ 2 4 6 7 14 26 27 28 29 9 10 12 13  15 16  18 19  22 23 24 25]];
            
        end
        %%%%%%%%%%%% identification vector, 1=delta, 2=regular/other, 0= low SNR
        dyn_type= [ 0 2 0 2 0 2 2 0 2 2    0 1 1 2 1 1 0 1 1 0    0 1 1 1 1 2 2 2 2 0     0 0 0 0];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end


ses= dir('*.mat')

for sesd=[chs] %[1:25 27:32]%length(ses);%[1:15 18:31] %:length(ses)  %22 exclude 16 17
    sesd
 
    load(ses(sesd).name)
try
%     
%    for id=1:length( aligned.trace_nosp)
%        
%        
%        
%        
%    end
   snr_thres= 0;

if length(aligned.trace)>1
   FS=1000;
   lfp=[];mk=0;
SN=[]
for in =1:length(aligned.trace)
    windF=length(aligned.trace{in});%ceil(FS*1);
    clear d A
   d(1,:)= aligned.lfp{in}(:,1:end); % LFP
  
    A=zeros( 1,length(aligned.mvmt{in}(1:end)));
    HM= aligned.f_mvmt{in}(1:end);HM(HM>=50)=1000;
    ;HM(HM<50)=-1000;
    % LM= aligned.mvmt{in}(1:end);;LM(LM<40)=1000;
    % LM(LM>=40)=-1000;
 
    d(4,:)=HM;   d(5,:)=HM.*-1;
  A(aligned.vectors{in,8})=1;
    B=zeros( 1,length(aligned.mvmt{in}));
    B(aligned.vectors{in,3})=1;
  
    A=aligned.spike{in}(:,1:end);A(isnan(A))=0;
    
    % A= spRES{in}.roaster;A(isnan(A))=0;
   d(2,:)=A;  %% SPIKES
   
    d(3,:)=  aligned.trace_nosp{in}(:,1:end) ;  %Vm
    
    
      d(6,:)=B;
      
        d(7,:)=  aligned.imgmotion{in};  % combined X-Y image motion
      d(7,[1 2 ])=d(7,3);
      d(7,end:-1:end)=d(7,end-2);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
   d(7,:)= [ 0 fastsmooth(abs(hilbert(diff(fastsmooth(d(7,:),10,1,1)))),100,1,1)];
  
      
        SN= [SN;spRES{in}.spike_snr{:}];;
   % d(1,:)= dataM.trial{in}(1,9:end);
   Nr= floor(length(d)./windF);
      if mean(spRES{in}.spike_snr{:}) >snr_thres
   for nh= 1:Nr
 mk=mk+1;
 lfp.trial{mk}(1:3,:) = zscore(d(1:3,(1  + ((windF)* (nh-1)):   windF * (nh))   ),[],2);
  lfp.trial{mk}(4:5,:) = (d(4:5,(1  + ((windF)* (nh-1)):   windF * (nh))   ));
  lfp.trial{mk}(6,:) = (d(6,(1  + ((windF)* (nh-1)):   windF * (nh))   ));
   lfp.trial{mk}(7,:) = (d(7,(1  + ((windF)* (nh-1)):   windF * (nh))   ));
  lfp.time{mk}= (1:size(lfp.trial{mk},2))./FS;
   end  ;end 

end


   n=0;
for ind=1:size(lfp.trial{1},1)
    n=n+1;
    lfp.label(n)= { [ 'CH' num2str(ind)]};
end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 
 
cfg=[];
 cfg.bpfilter      = 'yes'
 cfg.bpfreq    =[1 6]
  cfg.bpfiltord    =2;
  cfg.hilbert='angle'
lfp2=ft_preprocessing(cfg,lfp)

ifs=[];
for xc=1:5
    ph= lfp2.trial{xc}(3,:); 
frameSize=201;
polyOrder=2;
diffOrder=1;
inPhase=unwrap(ph,pi/4);
inFreq=savitzkyGolayFilt(inPhase,polyOrder,diffOrder,frameSize);
inFreq(end)=[];
tt=-1*1000.*inFreq/(2*pi);tt(tt<0.5 | tt>10)=NaN;
ifs=[ifs, tt(100:end-100)];
end
 [n1 n2]= hist(ifs,[0.8:0.2:5.2])
 
 %[ data,dataM, timed, tr]= delta_phase_TFR(freq2, lfp2,3,lfp,flipC);
 clear allc  allc2
for tr=1:4
F=aligned.spike{tr};F(isnan(F))=0;F=fastsmooth(F,10,1,1);
F1=zscore(aligned.trace_nosp{tr});%F(isnan(F))=0;

[c,lags]=xcorr(F,3000,'Coeff');[c1,lags]=xcorr(F1,3000,'Coeff');
allc(:,tr)=c(3003:end);allc2(:,tr)=c1(3003:end);
end

 %%%%%%

    %  figure,imagesc(abs(nanmean(exp(1i.*angle(aM)),3))); axis xy
    if 1%size(aM,3)>0
    nn=nn+1;
%allP(:,:,nn)=  nanmean(dataM(:,:,:),3);%nanmean(abs(aM),3);
%allPR(:,:,nn)=  nanmean(abs(aM2),3);
 allC(:,nn)=  n1(2:end-1)/sum(n1(2:end-1));%  (nanmean(exp(1i.*angle(aM)),3));%
 allC2(:,1,nn)=  nanmean(allc,2);%  (nanmean(exp(1i.*angle(aM)),3));%
 allC2(:,2,nn)=  nanmean(allc2,2);
%allS(:,nn)=fastsmooth(Sprob,4,1,1);

 %AC=[AC; permute(abs(aM),[3 1 2])];
%   AC=AC + nansum(abs(aM),3);
%     N1=N1+size(aM,3);
%  allCR2(:,:,nn)= nanmean(abs(aM2),3);;
%  % AC2=[AC2;  permute(abs(aM2),[3 1 2])];
%     AC2=AC2 + nansum(abs(aM2),3);
%     N2=N2+size(aM2,3);
%  allC2(:,:,nn)= nanmean(abs(aM),3);%abs(nanmean(abs(aM)./nanmean(abs(aM),3).*exp(1i.*angle(aM)),3)) ;
%   allC3(:,:,nn)=  nanmean(abs(aM3),3);%(nanmean(exp(1i.*angle(aM)),3));%
%  allCR4(:,:,nn)= nanmean(abs(aM4),3);;
% sigC(:,nn)= nanmean(sigM,2);
    end
    
    

    
    
end
end
end
end


%%%%
savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\git_scripts\fig4\'
pheight=160;
%%%%%%%%%%%%%%%%%%%%%%






% figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
% imagesc(timed,freq2.freq,smooth2a(bsxfun(@rdivide,nanmean(allC,3), nanmean(nanmean(allC,3),2)) ,1,3)  ); axis xy
% set(gca,'CLim',[ 0.9 1.1])
% colormap(jet)
timed=(1:0.2:5)*1;
Fcurve=nanmean(allC,2);Fcurve=Fcurve-min(Fcurve);
figure('COlor','w','Position', [ 300 400 180 160],'Renderer', 'painters')
plot(timed,Fcurve,'m','Linewidth',1)
fill_error_area2(timed,Fcurve,nanstd(allC,[],2)./sqrt(size(allC,2)),[0.2 0.2 0.2 ] )
axis tight
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'delta_fr_dist_' num2str(gh) '.pdf'])
%
clear FWHM FWHM2 peakF
for  fg=1:size(allC,2)
Fcurve=allC(:,fg);Fcurve=Fcurve-min(Fcurve);
Fcurve=Fcurve./max(Fcurve);
[n1 ]=minmax(timed(Fcurve>=0.5))
[t1 t2]=max(Fcurve);
FWHM(fg)= n1(2)-n1(1);
FWHM2(fg)= n1(2)./n1(1);
peakF(fg)= timed(t2);
end

mean(FWHM2)
std(FWHM2)
mean(FWHM)
std(FWHM)
mean(peakF)
std(peakF)
V2=FWHM
V3=FWHM2;V1=peakF;
  figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters') 
%b1=bar(1,nanmean(V1),'Facecolor',[ 0.4 0.4 0.4]);hold on,
h=boxplot( [V1 ;V2]'  ,[ ones(length(V1),1)+0 ;ones(length(V2),1).*2],  'notch','on', 'Widths',0.6 , 'colors',[ 0.2 0.2 0.2], 'symbol','.k')
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'bar_quant_delt_fres'  '.pdf'])

figure('COlor','w','Position', [ 300 400 120 pheight-10],'Renderer', 'painters')
violinplotSTR(V1',[1.3 ],'ViolinColor', [0.4 0.5 .9])
hold on,
violinplotSTR(V2',[ 1.9],'ViolinColor', [0.4 0.8 0.7]);ylim([ 1 3.2])
line([ 0.7 2.5], [ 0  0],'COlor', [0 0 0 ],'Linewidth',0.5)
ylim([ 1 3.2])
xlim([ 0.7 2.5])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'bar_quant_delt_fres'  '.pdf'])

figure('COlor','w','Position', [ 300 400 200 126],'Renderer', 'painters')
plot(lags(3003:end),nanmean(allC2(:,1,:),3),'k')
fill_error_area2(lags(3003:end)./1000,nanmean(allC2(:,1,:),3),nanstd(allC2(:,1,:),[],3)./sqrt(size(allC2,3)),[0.5 0.5 0.5])
axis tight
xlim([0 2500]);ylim([-0.05 0.2])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'autocrr_sp_' num2str(gh) '.pdf'])

figure('COlor','w','Position', [ 300 400 200 160],'Renderer', 'painters')
plot(lags(3003:end)./1000,nanmean(allC2(:,2,:),3),'r')
fill_error_area2(lags(3003:end)./1000,nanmean(allC2(:,2,:),3),nanstd(allC2(:,2,:),[],3)./sqrt(size(allC2,3)),[0.5 0.5 0.7])
plot(-lags(3003:end),nanmean(allC2(:,1,:),3),'k')
fill_error_area2(-lags(3003:end)./1000,nanmean(allC2(:,1,:),3),nanstd(allC2(:,1,:),[],3)./sqrt(size(allC2,3)),[0.5 0.5 0.5])
axis tight
xlim([-2.5 2.5]);ylim([-0.2 0.3])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'autocrr_VMpike_' num2str(gh) '.pdf'])


