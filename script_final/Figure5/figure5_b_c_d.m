
clear all
cell_id=[];Is_ratio3=[];d_id=[];

pT=[];pT2=[];  nn=0;sw_t=3;
snr_thres=5;
for gh=1:2
    
 if gh==1
        cd('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata2\Hua-an Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\ChAT\')
   dyn_type=   [ 0,2,1,1,1,1,1,1,1,0,1,2,1,1,1,1,1,0,1,1,1,1,1,1,1,2,1,2,0,2,0,1]  
            chs=[find(dyn_type>0)];
        %%%%%%%%%%%%%%
    else
        cd('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata2\Hua-an Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\MSN\')
        %%%%%%%%%%%% identification vector, 1=delta, 2=regular/other, 0= low SNR
 dyn_type= [2,2,2,2,0,2,2,0,2,2,0,1,1,2,0,1,2,1,1,2,0,1,1,1,1,2,2,2,2,0,2,0  0 0 0]
            chs=[[ find(dyn_type>0)]];
    end
ses= dir('*.mat')


for sesd= [chs]
    sesd
 
    load(ses(sesd).name)


if length(aligned.trace)>0
   FS=1000;
   lfp=[];mk=0;

for in =1:length(aligned.trace)
    wind=length(aligned.trace{in});%ceil(FS*1);
    clear d A
  d(1,:)= aligned.lfp{in}(:,1:end); % LFP
 
             %%%%%%%%%%%%%
                               
     A=aligned.spike{in}(:,1:end);A(isnan(A))=0;
 
   d(2,:)=A;
 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  d(4,:)=  aligned.imgmotion{in};  % combined X-Y image motion
      d(4,[1 2 ])=d(4,3);
      d(4,end:-1:end)=d(4,end-2);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
   d(4,:)= [ 0 fastsmooth(abs(hilbert(diff(fastsmooth(d(4,:),10,1,1)))),100,1,1)];
   
            A=zeros( 1,length(aligned.mvmt{in}(1:end)));
        HM= aligned.f_mvmt{in}(1:end);
        HM(HM>=50)=1000;
        HM(HM<50)=-1000;
        d(5,:)=HM;   
       d(6,:)=HM.*-1;
   
   Nr= floor(length(d)./wind);
    if mean(spRES{in}.spike_snr{:}) >=snr_thres% 
   for nh= 1:Nr
 mk=mk+1;
 lfp.trial{mk}(1:3,:) = zscore(d(1:3,(1  + ((wind)* (nh-1)):   wind * (nh))   ),[],2);
  lfp.trial{mk}(4,:) = (d(4,(1  + ((wind)* (nh-1)):   wind * (nh))   ));
 lfp.trial{mk}(5:6,:) = (d(5:6,(1  + ((wind)* (nh-1)):   wind * (nh))   ));
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
 
           warning off
  cfg = []; %
    cfg.method ='wavelet'; %
    cfg.output ='fourier';
     cfg.taper='hanning';
    cfg.keeptapers ='yes';
    cfg.keeptrials ='yes';
    cfg.trials='all';cfg.tapsmofrq =6;%
     cfg.channel= ['all']; %chans=cfg.channel;
    cfg.foi= [1:1:100];
     cfg.toi=lfp.time{1}(1:1:end) ;
     cfg.width =5;
freq2 = ft_freqanalysis(cfg, lfp);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 
 
 wind=10;
    clear aM aM2
    mm=0;mm2=0;
    for id=1:size(  freq2.fourierspctrm,1)
M=(squeeze(freq2.fourierspctrm(id,1,:,:)))      ;
M2=(squeeze(freq2.fourierspctrm(id,1,:,:)))      ;
s= lfp.trial{id}(2,:)>0;s=find(s);
high_mot= lfp.trial{id}(5,:)>0; low_mot= lfp.trial{id}(6,:)>0; 
for id2=1:length(s)

    if s(id2) >wind  & s(id2) +wind < length(M)  & high_mot(s(id2))>0
        mm=mm+1;
aM(:,:,mm)=   M(:,s(id2)-wind:s(id2)+wind);
    end
end

for id2=1:length(s)
    if s(id2) >wind  & s(id2) +wind < length(M)  & low_mot(s(id2))>0
        mm2=mm2+1;
        aM2(:,:,mm2)=  M2(:,s(id2)-wind:s(id2)+wind);
    end
end

    end
    

   
    if sum(~isnan(aM(1, wind+1,:)))>10  & sum(~isnan(aM2(1,wind+1,:)))>10
    nn=nn+1;
allP(:,:,nn)=  nanmean(abs(aM),3);
 Z= abs(nanmean(abs(aM)./nanmean(abs(aM),3).*exp(1i.*angle(aM)),3)) ;
 for b1=1:size(aM,1)
 NT= sum(~isnan(aM(b1,wind,:)));
  T=Z(b1,:).^2; 
  allC2(b1,:,nn)= (((1/(NT-1))*((T.*NT-1)))) ; %high motion
 end 
 
  Z= abs(nanmean(abs(aM2)./nanmean(abs(aM2),3).*exp(1i.*angle(aM2)),3)) ;
 for b1=1:size(aM2,1)
 NT= sum(~isnan(aM2(b1,wind,:)));
  T=Z(b1,:).^2; 
  allC(b1,:,nn)= (((1/(NT-1))*((T.*NT-1)))) ; % low motopm
 end 
     
  cell_id=[cell_id,gh];
       d_id=[d_id,  dyn_type(sesd) ];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%
   

  %  end 

end
end
end

% figure('COlor','w')
% imagesc([],freq2.freq,(nanmean((allP),3)).*repmat(freq2.freq.^0.1,size(allP,2),1)'  );colormap(jet)
% axis xy
savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\git_scripts\cohVm\'
pheight=160;

tt=-wind:wind;
tsel=tt<-100;




nr=0;clear allCOH allCOHs allCOH2 allCOH2s
for x=freq2.freq
nr=nr+1;
id=wind+1;
allCOH(nr)= nanmean(nanmean((allC(nr,id,:)),3),2);
allCOHs(nr)= nanmean(nanstd((allC(nr,id,:)),[],3),2)./sqrt(size(allC,3));
allCOH2(nr)= nanmean(nanmean((allC2(nr,id,:)),3),2);
allCOH2s(nr)= nanmean(nanstd((allC2(nr,id,:)),[],3),2)./sqrt(size(allC2,3));
end
%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id=wind+1;
clear allst
for t=1:1000
    t
for b1=1:size(allC,1)
plvs=[squeeze(allC(b1,id,:)); squeeze(allC2(b1,id,:))];
order_p=[ ones(1,size(allC,3)),ones(1,size(allC,3)).*2];
z=randperm(length(plvs));
plvs=plvs(z);
[h,p,ci,stats] =ttest(plvs(order_p==1), plvs(order_p==2));
allst(b1,t)=abs(stats.tstat);
end
end
clear allstT
for b1=1:size(allC,1)
[h,p,ci,stats] =ttest(squeeze(allC(b1,wind+1,:)), squeeze(allC2(b1,wind+1,:)));
allstT(b1)=abs(stats.tstat);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coh_SH_99=prctile(allst',95);
%%%%
sign_line= zeros(1,size(allCOH,2)).*NaN;
 sign_line(find((allstT-coh_SH_99)>0))=1;
savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\git_scripts\fig4\'
pheight=160;
figure('COlor','w','Position', [ 300 300 200 pheight])
plot(freq2.freq,allCOH,'b','COlor',[0.9 0.5 0]);fill_error_area2(freq2.freq, allCOH, allCOHs,[0.5 0.5 0.5])
plot(freq2.freq,allCOH2,'r','COlor',[0.4 0.8 0]);fill_error_area2(freq2.freq, allCOH2, allCOH2s,[0.5 0.5 0.5])
hold on,plot(freq2.freq,(sign_line.*max(allCOH2(:)))*1.5,'y.-','Linewidth',3,'COlor',[ 1 0.85 0])
ylim([min(allCOH2(:)) max(allCOH2(:))*1.4]);
plot(freq2.freq,allCOH,'b','COlor',[0.9 0.5 .9],'Linewidth',1);fill_error_area2(freq2.freq, allCOH, allCOHs,[0.5 0.5 0.5])
plot(freq2.freq,allCOH2,'r','COlor',[0.4 0.8 0],'Linewidth',1);fill_error_area2(freq2.freq, allCOH2, allCOH2s,[0.5 0.5 0.5])
set(gca,'Xscale','log')
xlim([1 100]);ylim([ -0.005 0.05])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'SFC_LFP_low_high_' num2str(sw_t) '.pdf'])



%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 fsel=find(round(freq2.freq)>=1 & round(freq2.freq)<=4);% 2Hz peak selection
 selC=cell_id==1&d_id==1;
V1= squeeze(nanmean(nanmean((allC2( fsel,id, selC)),2),1))- squeeze(nanmean(nanmean((allC( fsel,:, selC)),2),1));
  selC=cell_id==1&d_id==2;
 V2= squeeze(nanmean(nanmean((allC2( fsel,id, selC)),2),1))- squeeze(nanmean(nanmean((allC( fsel,id, selC)),2),1));
 selC=cell_id==2&d_id==1;
V3= squeeze(nanmean(nanmean((allC2( fsel,id, selC)),2),1))- squeeze(nanmean(nanmean((allC( fsel,id, selC)),2),1));
 selC=cell_id==2&d_id==2;
V4= squeeze(nanmean(nanmean((allC2( fsel,id, selC)),2),1))- squeeze(nanmean(nanmean((allC( fsel,id, selC)),2),1));
 
figure('COlor','w','Position', [ 300 400 180 pheight],'Renderer', 'painters') 
b1=bar(1,nanmean(V1),'Facecolor',[ 0.8 0 0]);hold on,
set(b1,'FaceAlpha',0.4)
 % b1=bar(2,nanmean(V2),'Facecolor',[ 0 0 0.8]);hold on,
set(b1,'FaceAlpha',0.4)
  b1=bar(4,nanmean(V3),'Facecolor',[ 0.8 0 0]);hold on,
set(b1,'FaceAlpha',0.4)
  b1=bar(5,nanmean(V4),'Facecolor',[ 0 0 0.8]);hold on,
set(b1,'FaceAlpha',0.4)
errorbar([1   ],nanmean(V1,1), nanstd(V1)./sqrt(size(V1,1)),'.k');
errorbar([4   ],nanmean(V3,1), nanstd(V3)./sqrt(size(V3,1)),'.k');
errorbar([5   ],nanmean(V4,1), nanstd(V4)./sqrt(size(V4,1)),'.k');



figure('COlor','w','Position', [ 300 400 180 150],'Renderer', 'painters')
violinplotSTR(V1,[1.3 ],'ViolinColor', [ 0.9 0 0.0])
hold on,
%violinplot2(V2,[1.7 ],'ViolinColor', [ 0. 0. 0.9])
violinplotSTR(V3,[ 2.1],'ViolinColor', [ 0.9 0. 0.])
violinplotSTR(V4,[ 2.6],'ViolinColor', [ 0. 0. 0.9])
line([ 0.8 2.9], [ 0  0],'COlor', [0 0 0 ],'Linewidth',0.5)
 print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'LowHIgh_COH_BARQUANT.pdf'])


 [h,p,ci,stats] =ttest(V1)
% df=20,  0.0106   0.0146   0.0082
 [h,p,ci,stats] =ttest(V3)
% df=8, 0.0756    0.0396  0.0677
[h,p,ci,stats] =ttest(V4)
% df=15, 0.19   0.9585    0.9114

     [h,p,ci,stats] =ttest2(V1,V4)
%df=35, tst=2.7,         0.0478     0.0524     0.0267
   [h,p,ci,stats] =ttest2(V3,V4)
%df=23, tst=3.8,         0.1967     0.1240    0.0886
   [h,p,ci,stats] =ttest2(V1,V3)
%df=28, tst=3.8,         0.4753   0.6071   0.5859
 