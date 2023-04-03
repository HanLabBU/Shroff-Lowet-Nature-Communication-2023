%addpath(genpath('\\engnas.bu.edu\Research\eng_research_handata\EricLowet\Scripts\'))
%addpath(genpath('Z:\EricLowet\Scripts\'))

clear all
pT=[];pT2=[];  nn=0;
 INT1=[];INT2=[];
sw_t=2;
for gh=[1 2 ]
    
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
for sesd=chs;%
    sesd
  HIGHPOW=[];LOWPOW=[];clear spRES
    load(ses(sesd).name)
try
   
%    end
   snr_thres=3;

if length(aligned.trace)>0
   FS=830;
   lfp=[];mk=0;
SN=[]
for in =1:length(aligned.trace)
    windF=length(aligned.trace{in});%ceil(FS*1);
    clear d A
    A= spRES{in}.roaster;A(isnan(A))=0;
   d(1,:)= aligned.lfp{in}(:,1:end); % LFP
  
    A=zeros( 1,length(aligned.mvmt{in}(1:end)));
    HM= aligned.f_mvmt{in}(1:end);HM(HM>=50)=1000;
    ;HM(HM<50)=-1000;
    % LM= aligned.mvmt{in}(1:end);;LM(LM<40)=1000;
    % LM(LM>=40)=-1000;
 
    d(4,:)=HM;   d(5,:)=HM.*-1;
  %A(aligned.vectors{in,4})=1;  %
  
  
  %  A=aligned.spike{in}(:,1:end);A(isnan(A))=0;
    
     A= spRES{in}.roaster;A(isnan(A))=0;
       B=zeros( 1,length(aligned.mvmt{in}));    %MOTION
 %  B(aligned.vectors{in,3})=1;
   d(2,:)=A;%B;  %% SPIKES
    d(3,:)= zscore(aligned.trace{in}(:,1:end));  %Vm
        SN= [SN;spRES{in}.spike_snr{:}];;
   % d(1,:)= dataM.trial{in}(1,9:end);
   Nr= floor(length(d)./windF);
      if mean(spRES{in}.spike_snr{:}) >snr_thres
   for nh= 1:Nr
 mk=mk+1;
 lfp.trial{mk}(1:3,:) = zscore(d(1:3,(1  + ((windF)* (nh-1)):   windF * (nh))   ),[],2);
  lfp.trial{mk}(4:5,:) = (d(4:5,(1  + ((windF)* (nh-1)):   windF * (nh))   ));
  lfp.time{mk}= (1:size(lfp.trial{mk},2))./FS;
   end  ;end 

end


   n=0;
for ind=1:size(lfp.trial{1},1)
    n=n+1;
    lfp.label(n)= { [ 'CH' num2str(ind)]};
end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 
           warning off
  cfg = []; %block_type == cfg.blk
    cfg.method ='wavelet'; %'mvar';
    cfg.output ='fourier';
     cfg.taper='hanning';
    cfg.keeptapers ='yes';
    cfg.keeptrials ='yes';
    cfg.trials='all';cfg.tapsmofrq =5;%
     cfg.channel= [3 ]; %chans=cfg.channel;
    cfg.foi= [0.5:30.5:100];
     cfg.toi=lfp.time{1}(1:1:end) ;
     cfg.width =5;
    cfg.t_ftimwin =[ones(1,length(cfg.foi))*0.5];
freq2 = ft_freqanalysis(cfg, lfp);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 REF=1;

for in=1:length(spRES)

A=spRES{in}.roaster;A(isnan(A))=0;

v=find(A);vx= v(3:end)-v(2:end-1);vx2=v(2:end)-v(1:end-1)  ;
INT1=[INT1, vx(1:end) ];
INT2=[INT2, vx2(1:end-1)  ];
end



clear spikTyp
wind=51;
for tria = 1:length(spRES)
   s=  spRES{tria}. spike_idx{1}  ;
    for in = 1:length(s)
        if length(s) >1
       if in==1
        d=(s(in+1)-  s(in));
        if d <wind;   spikTyp{tria}(in)= 1 ;
        else; spikTyp{tria}(in)= 0 ;end
      
       elseif in == length(s)
                 d=( s(in)-s(in-1) );
        if d <wind;   spikTyp{tria}(in)= 1 ;
        else; spikTyp{tria}(in)= 0 ;end  
       else
                d1=(s(in+1)-  s(in));         d=( s(in)-s(in-1) );
        if d <wind |d1 <wind ;   spikTyp{tria}(in)= 1 ;
        else; spikTyp{tria}(in)= 0 ;end 
               
       end
        else
          spikTyp{tria}=[];  
    end
    end
end



 wind=600;alls=[];
    clear aM aM2  aM3 aM4 
    mm=0;   mm2=0; mm3=0; mm4=0;
    for id=1:size(  freq2.fourierspctrm,1)
M=(squeeze(freq2.fourierspctrm(id,REF,:,:)))      ;
s= lfp.trial{id}(2,:)>0;%s=find(s);
%alls=[alls ,s];
HT=find(lfp.trial{id}(4,:)>0);LT=find(lfp.trial{id}(5,:)>0);
%vv= lfp.trial{id}(3,:);
  NT=length(s); %TL= size(MG,2);
   

HIGHPOW=[HIGHPOW,s(HT(5:end-5))];
LOWPOW=[LOWPOW,s(LT(5:end-5))];

    end
    
    

    

%     %  figure,imagesc(abs(nanmean(exp(1i.*angle(aM)),3))); axis xy
%     if size(aM,3)>1
if length(HIGHPOW)> 10000  & length(LOWPOW) >10000
     nn=nn+1;
     allP(:,1,nn)= nanmean(abs(HIGHPOW),2);
     allP(:,2,nn)= nanmean(abs(LOWPOW),2);
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     allFIR(nn,1)=(sum(HIGHPOW)./length(HIGHPOW)).*828;
     allFIR(nn,2)=(sum(LOWPOW)./length(LOWPOW)).*828;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
     [n1 n2]=hist(diff(find(HIGHPOW)),[2:2:180]);
     n2=n2*1.2;
    allFIR_IB(nn,1) = 1000/(sum(n2(1:end-1).*n1(1:end-1))./sum(n1(1:end-1)))
    
      [n1 n2]=hist(diff(find(LOWPOW)),[2:2:180]);
     n2=n2*1.2;
    allFIR_IB(nn,2) = 1000/(sum(n2(1:end-1).*n1(1:end-1))./sum(n1(1:end-1)))
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
        [n1 n2]=hist(diff(find(HIGHPOW)),[400:4:1800]);
     n2=n2*1.2;
    allFIR_BB(nn,1) = 1000/(sum(n2(1:end-1).*n1(1:end-1))./sum(n1(1:end-1)))
    
        [n1 n2]=hist(diff(find(LOWPOW)),[400:4:1800]);
     n2=n2*1.2;
    allFIR_BB(nn,2) = 1000/(sum(n2(1:end-1).*n1(1:end-1))./sum(n1(1:end-1)))
    
    
    
end
% allP(:,:,nn)=  nansum(abs(aM),3);
% trwin(nn)= size(aM,3);trwin2(nn)= size(aM2,3);
% %allPR(:,:,nn)=  nanmean(abs(aM2),3);
%  allC(:,:,nn)=  nanmean(abs(aM),3);%       (nanmean(exp(1i.*angle(aM)),3));%
%  %allC(:,:,nn)=      (nanmean(exp(1i.*angle(aM)),3));%
%  allCR2(:,:,nn)= nansum(abs(aM2),3);;
%  allC2(:,:,nn)= nanmean(abs(aM),3);%abs(nanmean(abs(aM)./nanmean(abs(aM),3).*exp(1i.*angle(aM)),3)) ;
%   allC3(:,nn)=  fastsmooth(nanmean((aM3),2),1,1,1);%(nanmean(exp(1i.*angle(aM)),3));%
%  allCR4(:,:,nn)= nanmean(abs(aM4),3);;
%     end
    
    
 
            
        end
%
end
    


end
end


% save('nw_data_shuffle_motion_CHAT_off')
%%%%
savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\git_scripts\fig2\'
pheight=160;
%%%%%%%%%%%%%%%%%%%%%%
  figure('COlor','w','Position', [ 300 200 400 pheight],'Renderer', 'painters') 
  subplot(1,2,1)
hist(INT2,[2:3:120]),%set(gca,'Xscale','log');
xlim([ 1 110])
h = findobj(gca,'Type','patch');
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = [0 0.2 0.2];
  subplot(1,2,2)
hist(INT2,[10:10:3000]),%set(gca,'Xscale','log');
xlim([ 40 2000])
h = findobj(gca,'Type','patch');
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = [0 0.2 0.2];

%   figure('COlor','w','Position', [ 300 200 250 pheight],'Renderer', 'painters') 
%   subplot(1,1,1),plot(INT1(1:3:end)+randn(1,length(INT1(1:3:end))).*0.3,INT2(1:3:end)+randn(1,length(INT1(1:3:end))).*0.3,'.k','Markeralpha',0.5)
% set(gca,'Xscale','log','Yscale','log')
% axis tight;axis([ 5 2000 5 2000])
% line([ 5 2000],[5 2000],'COlor',[ 0.8 0 0],'Linewidth',1)
% 
  figure('COlor','w','Position', [ 300 200 200 pheight],'Renderer', 'painters') 
hist(INT2,[2:10:1012]),set(gca,'Xscale','log');axis tight
xlim([ 7 1000])
h = findobj(gca,'Type','patch');
h.FaceColor = [0.5 0.5 0.5];
h.EdgeColor = [0.2 0.2 0.2];
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'ISI_2log_' num2str(sw_t) '.pdf'])


  figure('COlor','w','Position', [ 300 200 200 pheight],'Renderer', 'painters') 
hist(INT1,[2:10:1012]),set(gca,'Xscale','log');axis tight
xlim([ 7 1000])
h = findobj(gca,'Type','patch');
h.FaceColor = [0.5 0.5 0.5];
h.EdgeColor = [0.2 0.2 0.2];
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'ISI_1log_' num2str(sw_t) '.pdf'])
xlim([100 1000])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'ISI_1log_' num2str(sw_t) 'zoom.pdf'])


 figure('COlor','w','Position', [ 300 200 250 pheight],'Renderer', 'painters') ;
scatter1 = scatter(INT1(1:4:end)+randn(1,length(INT1(1:4:end))).*0.3,INT2(1:4:end)+randn(1,length(INT1(1:4:end))).*0.3,5,'MarkerFaceColor','k','MarkerEdgeColor','k'); 
scatter1.MarkerFaceAlpha = .7;
scatter1.MarkerEdgeAlpha = .1;
set(gca,'Xscale','log','Yscale','log')
axis tight;axis([ 5 1500 5 1500])
line([ 5 1500],[5 1500],'COlor',[ 0.8 0 0],'Linewidth',1)
box on
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'ISI_return_map_' num2str(sw_t) '.pdf'])

% clear allC
% n1=0;figure
% for ix=2:wind:1000
%     n1=n1+1;n2=0;
%     for ix2=2:wind:1000
%         n2=n2+1;zz=find((INT1>ix  & INT1 <ix+wind) &  (INT2>ix2  & INT2 <ix2+wind) );
%    nn1= length(  find((INT1>ix  & INT1 <ix+wind) &  (INT2>ix2  & INT2 <ix2+wind) ) );
%   allC(n1,n2)=nn1;
%   hold on, plot(INT1(zz)+randn(1,length(INT1(zz))).*0.3,INT2(zz)+randn(1,length(INT1(zz))).*0.3,'.k')
%     end
% end
% set(gca,'Xscale','log','Yscale','log')
% axis([ 5 1000 5 1000])


% figure,imagesc(allC)
% axis xy
% 
% figure
% subplot(1,1,1);
% plot(aligned.trace{1},'k');axis tight
% % subplot(2,1,2);
% plot(aligned.trace{2});axis tight



Vx=aligned.trace{1};Vx=Vx-fastsmooth(Vx,5000,1,1);

%figure,hist(Vx,30)
A=[];B=[];
for id=1:length(aligned.trace_nosp)
   vv= aligned.spike{id};
   vv(isnan(vv))=0;
A=[A, zscore((aligned.trace_nosp{id}))];
B=[B, zscore(vv)];
end

% 
% [c,lags]=xcorr(zscore(A),zscore(B),1000,'Coeff');
% 
%   figure('COlor','w','Position', [ 300 400 290 pheight],'Renderer', 'painters') 
%  % c(1000:1001)=0;
% plot(lags.*1.2,fastsmooth(c,2,1,1),'k')


 [h,p,ci,stats] = ttest(allFIR(:,1),allFIR(:,2))
 
 [h,p,ci,stats] = ttest(allFIR_IB(:,1),allFIR_IB(:,2))
  
   %[h,p,ci,stats] = ttest(allFIR_BB(:,1),allFIR_BB(:,2))
 
  V1=allFIR(:,1);
 V2=allFIR(:,2);
  V1b=allFIR_IB(:,1);
  V2b=allFIR_IB(:,2);
 V1c=allFIR_BB(:,1);
  V2c=allFIR_BB(:,2);
   M=[V1,V2, V1b, V2b];
     figure('COlor','w','Position', [ 300 400 180 160],'Renderer', 'painters') 
b1=bar(1,nanmean(V1),'Facecolor',[ 0.8 0 0]);hold on,
set(b1,'FaceAlpha',0.4)
  b1=bar(2,nanmean(V2),'Facecolor',[ 0 0 0.8]);hold on,
set(b1,'FaceAlpha',0.4)
  b1=bar(4,nanmean(V1b),'Facecolor',[ 0.8 0 0]);hold on,
set(b1,'FaceAlpha',0.4)
  b1=bar(5,nanmean(V2b),'Facecolor',[ 0. 0 0.8]);hold on,
set(b1,'FaceAlpha',0.4)
%   b1=bar(7,nanmean(V1c),'Facecolor',[ 0.8 0 0]);hold on,
% set(b1,'FaceAlpha',0.4)
%   b1=bar(8,nanmean(V2c),'Facecolor',[ 0. 0 0.8]);hold on,
% set(b1,'FaceAlpha',0.4)
errorbar([1 2 4 5 ],nanmean(M,1), nanstd(M)./sqrt(size(M,1)),'.k');
  
 [h,p,ci,stats] = ttest(V1./V2,1)
 [h,p,ci,stats] = ttest(V1b./V2b,1)