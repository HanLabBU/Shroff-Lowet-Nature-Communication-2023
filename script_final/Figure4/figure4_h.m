
%addpath(genpath('Z:\EricLowet\'))

clear all
pT=[];pT2=[];  nn=0;
cell_id=[];Is_ratio3=[];d_id=[];
allM=[];allM2=[]; sw_t=3;
snr_thres=5;
for gh=[ 1 2]
   if gh==1
        %cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Hua-an_Tseng/Data/SomArchon_Striatum/Linear_movement_PinnedBall/good_data_ALL/ChAT/')
        cd('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata2\Hua-an Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\ChAT\')
       %         cd('Z:\eng_research_handata2\Hua-an_Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\ChAT\')

dyn_type= [ 0,2,1,1,1,1,1,1,1,0,1,2,1,1,1,1,1,0,1,1,1,1,1,1,1,2,1,2,0,2,0,1] 
 
 if sw_t==1 % delta-identified
            chs=[find(dyn_type==1)]
        elseif sw_t==2 % regular-identified 
            chs=  [ find(dyn_type==2) ];
        elseif sw_t==3   % all
            chs=[find(dyn_type>0)];
            
        end
        %%%%%%%%%%%%%%
    else
        %cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Hua-an_Tseng/Data/SomArchon_Striatum/Linear_movement_PinnedBall/good_data_ALL/MSN/')
        cd('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata2\Hua-an Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\MSN\')
       %        cd('Z:\eng_research_handata2\Hua-an_Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\MSN\')

dyn_type= [2,2,2,2,0,2,2,0,2,2,0,1,1,2,0,1,2,1,1,2,0,1,1,1,1,2,2,2,2,0,2,0  0 0 0]

        if sw_t==1      % delta-identified
           chs=[find(dyn_type==1)];
        elseif sw_t==2     % regular-identified or theta bursting
            chs=  [ find(dyn_type==2) ]
        elseif sw_t==3      % all
            chs=[[find(dyn_type>0)]];
            
        end
        
    end


ses= dir('*.mat')

for sesd=chs   
    sesd
 
    load(ses(sesd).name)


if size(aligned.vectors,1)>0
   FS=1000;
   lfp=[];mk=0;

for in =1:size(aligned.vectors,1)%length(aligned.trace)
    wind=length(aligned.trace{in});%ceil(FS*1);
    clear d A
  d(1,:)= aligned.lfp{in}(:,1:end); % LFP
  
  % d(1,:)= spRES{in}.denoise_trace;%
   vv= aligned.trace_nosp{in}(:,1:end);  %Vm
             %%%%%%%%%%%%%
                                        %% denoise
                                         Fn = FS/2;FB=[ 75 88];
                                         [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                                         LFPg= ((filtfilt(B,A,  vv)));
                                     vv= vv-LFPg;
    H=aligned.spike{in}(:,1:end);H(isnan(H))=0;
     A=aligned.spike{in}(:,1:end);A(isnan(A))=0;
  %     A= spRES{in}.roaster;A(isnan(A))=0;
   d(2,:)=A;
   % d(1,:)= dataM.trial{in}(1,9:end);
         B=zeros( 1,length(aligned.mvmt{in}));    %MOTION
    B(aligned.vectors{in,3})=1;
   d(7,:)=B;
            B=zeros( 1,length(aligned.mvmt{in}));    %MOTION
    B(aligned.vectors{in,4})=1;
   d(8,:)=B;
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
                       if mean(spRES{in}.spike_snr{:}) >snr_thres% 

   for nh= 1:Nr
 mk=mk+1;
 lfp.trial{mk}(1:3,:) = zscore(d(1:3,(1  + ((wind)* (nh-1)):   wind * (nh))   ),[],2);
  lfp.trial{mk}(4,:) = (d(4,(1  + ((wind)* (nh-1)):   wind * (nh))   ));
 lfp.trial{mk}(5:6,:) = (d(5:6,(1  + ((wind)* (nh-1)):   wind * (nh))   ));
 lfp.trial{mk}(7,:) = (d(7,(1  + ((wind)* (nh-1)):   wind * (nh))   ));
 lfp.trial{mk}(8,:) = (d(8,(1  + ((wind)* (nh-1)):   wind * (nh))   ));

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
 
 
 
 wind=400;
    clear aM aM2
    mm=0;mm2=0;
    for id=1:size( lfp.trial,2)
s= lfp.trial{id}(7,:)>0;s=find(s);
ss= lfp.trial{id}(8,:)>0;ss=find(ss);
s2= lfp.trial{id}(2,:)>0;%s=find(s);
 imag_mot= lfp.trial{id}(4,:); 
high_mot= lfp.trial{id}(5,:)>0; low_mot= lfp.trial{id}(6,:)>0; 

  NT=length(s); TL= length(s2);
    
                VV=randperm(TL);VV(VV<wind | VV+wind > TL)=[];
                sh=VV(1:NT);
                
for id2=1:length(s)
    if s(id2) >wind  & s(id2) +wind < length(s2)  %& high_mot(s(id2))>0
        mm=mm+1;
aM(:,mm)=  fastsmooth(s2(s(id2)-wind:s(id2)+wind),100,1,1);
    end
end
% for id2=1:length(ss)
%     if ss(id2) >wind  & ss(id2) +wind < length(s2)  %& low_mot(s(id2))>0
%         mm2=mm2+1;
% aM2(:,mm2)=  fastsmooth(s2( sh(id2)-wind: sh(id2)+wind),50,1,1);
%     end
% end

    end
    try
      allM=[allM, aM];
    allM2=[allM2, aM];end


 %  figure,imagesc(nanmean(abs(aM),3)); axis xy
    
    %  figure,imagesc(abs(nanmean(exp(1i.*angle(aM)),3))); axis xy
    if 1 %sum(~isnan(aM(2,101,:)))>20  & sum(~isnan(aM2(2,101,:)))>20
    nn=nn+1;
%allP(:,:,nn)=  nanmean(abs(aM),3);
try
allC(:,nn)=nansum(aM,2);
%allC2(:,nn)=nansum(aM2,2);
alltrig(nn)=size(aM,2);
%alltrig2(nn)=size(aM2,2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iter=1:500
    iter
 wind=400;
    clear aM aM2
    mm=0;mm2=0;
    for id=1:size( lfp.trial,2)
s= lfp.trial{id}(7,:)>0;s=find(s);
%ss= lfp.trial{id}(8,:)>0;ss=find(ss);
s2= lfp.trial{id}(2,:)>0;%s=find(s);
  NT=length(s); TL= length(s2);
                VV=randperm(TL);VV(VV<wind | VV+wind > TL)=[];
                sh=VV(1:NT);               
for id2=1:length(sh)
    if sh(id2) >wind +1 & sh(id2) +wind < length(s2)  %& low_mot(s(id2))>0
        mm2=mm2+1;
aM2(:,mm2)=  fastsmooth(s2( sh(id2)-wind: sh(id2)+wind),50,1,1);
    end
end
    end
  allC3(:,nn,iter)=nansum(aM2,2);
alltrig3(nn,iter)=size(aM2,2);
    end
 
    
    

    end
       cell_id=[cell_id,gh];
       d_id=[d_id,  dyn_type(sesd) ];
    %%%%%%%%%%%%%%%%%%%%%%
   

  %  end 
    
 
end
end
end


%cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\git_scripts\fig4\')
savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\git_scripts\fig4\'
%savepath='Z:\EricLowet\git_scripts\fig4\'

pheight=160;

tt=-wind:wind;
tsel=tt<-100;

figure,plot(tt,fastsmooth(nanmean(allM,2),100,1,1))
hold on,plot(tt,fastsmooth(nanmean(allM2,2),100,1,1))

 selC=d_id==2;
clear MM
for fg=1:size(allC3,3)
MM(:,fg)= nanmean(bsxfun(@rdivide, squeeze((allC3(:,selC,fg))),(alltrig3(selC,fg))' ),2);
end
M_97=prctile(MM',97.5);
M_3=prctile(MM',2.5);
  figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
 M= fastsmooth(nanmean(allC(:,selC)./alltrig(selC),2),50,1,1).*1000;
 Ms= fastsmooth(nanstd(allC(:,selC)./alltrig(selC),[],2)./sqrt(length(find(selC))),50,1,1).*1000;
plot(tt,M,'k','Linewidth',1.4,'COlor',[ 0.2 0.2 0.2]); hold on,
%fill_error_area2(tt,M,Ms,[0.5 0.5 0.5]);hold on,
plot(tt, nanmean(MM.*1000,2),'k','Linewidth',1,'COlor',[1 1 1 ])
fill_error_area2(tt, nanmean(MM.*1000,2), nanstd(MM.*1000,[],2).*2,[1 0.96 0.9]);hold on,
fill_error_area2(tt,M,Ms,[0.6 0.6 0.6]);hold on,
plot(tt,M,'k','Linewidth',1.5,'COlor',[ 0.2 0.2 0.2]); hold on,
xlim([-300 300]);ylim([ 2 16])
x2 = [tt, fliplr(tt)];
 print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'REG_motion_onset.pdf'])
 

 
 selC=d_id==1;%&cell_id==2;%& cell_id==1;
clear MM
for fg=1:size(allC3,3)
MM(:,fg)= nanmean(bsxfun(@rdivide, squeeze((allC3(:,selC,fg))),(alltrig3(selC,fg))' ),2);
end
M_97=prctile(MM',97.5);
M_3=prctile(MM',2.5);
  figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
 M= fastsmooth(nanmean(allC(:,selC)./alltrig(selC),2),50,1,1).*1000;
 Ms= fastsmooth(nanstd(allC(:,selC)./alltrig(selC),[],2)./sqrt(length(find(selC))),50,1,1).*1000;
plot(tt,M,'k','Linewidth',1.4,'COlor',[ 0.2 0.2 0.2]); hold on,
%fill_error_area2(tt,M,Ms,[0.5 0.5 0.5]);hold on,
plot(tt, nanmean(MM.*1000,2),'k','Linewidth',1,'COlor',[1 1 1 ])
fill_error_area2(tt, nanmean(MM.*1000,2), nanstd(MM.*1000,[],2).*2,[1 0.96 0.9]);hold on,
fill_error_area2(tt,M,Ms,[0.6 0.6 0.6]);hold on,
plot(tt,M,'k','Linewidth',1.5,'COlor',[ 0.2 0.2 0.2]); hold on,
xlim([-300 300]);ylim([ 2 16])
x2 = [tt, fliplr(tt)];
 print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Delta_motion_onset_CHAT.pdf'])
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
selwinR= round(100*1.2);selwinL= round(0*1.2);
selC=d_id==1& cell_id==1;
V1=(nansum(allC(wind+selwinL:wind+selwinR,selC),1)./alltrig(selC))-(nansum(allC(wind-selwinR:wind-selwinL,selC),1)./alltrig(selC));
selC=d_id==2& cell_id==1;
V2=(nansum(allC(wind+selwinL:wind+selwinR,selC),1)./alltrig(selC))-(nansum(allC(wind-selwinR:wind-selwinL,selC),1)./alltrig(selC));
selC=d_id==1& cell_id==2;
V3=(nansum(allC(wind+selwinL:wind+selwinR,selC),1)./alltrig(selC))-(nansum(allC(wind-selwinR:wind-selwinL,selC),1)./alltrig(selC));
selC=d_id==2& cell_id==2;
V4=(nansum(allC(wind+selwinL:wind+selwinR,selC),1)./alltrig(selC))-(nansum(allC(wind-selwinR:wind-selwinL,selC),1)./alltrig(selC));
% 
%   figure('COlor','w','Position', [ 300 400 180 pheight],'Renderer', 'painters') 
% b1=bar(1,nanmean(V1),'Facecolor',[ 0.8 0 0]);hold on,
% set(b1,'FaceAlpha',0.4)
%   b1=bar(2,nanmean(V2),'Facecolor',[ 0 0 0.8]);hold on,
% set(b1,'FaceAlpha',0.4)
%   b1=bar(4,nanmean(V3),'Facecolor',[ 0.8 0 0]);hold on,
% set(b1,'FaceAlpha',0.4)
%   b1=bar(5,nanmean(V4),'Facecolor',[ 0 0 0.8]);hold on,
% set(b1,'FaceAlpha',0.4)
% errorbar([1   ],nanmean(V1), nanstd(V1)./sqrt(size(V1,2)),'.k');
% errorbar([2   ],nanmean(V2), nanstd(V2)./sqrt(size(V2,2)),'.k');
% errorbar([4   ],nanmean(V3), nanstd(V3)./sqrt(size(V3,2)),'.k');
% errorbar([5   ],nanmean(V4), nanstd(V4)./sqrt(size(V4,2)),'.k');

figure('COlor','w','Position', [ 300 400 180 150],'Renderer', 'painters')
violinplotSTR(V1',[1.3 ],'ViolinColor', [ 0.9 0 0.0])
hold on,
%violinplot2(V2,[1.7 ],'ViolinColor', [ 0. 0. 0.9])
violinplotSTR(V3',[ 2.1],'ViolinColor', [ 0.9 0. 0.])
violinplotSTR(V4',[ 2.6],'ViolinColor', [ 0. 0. 0.9])
line([ 0.8 2.9], [ 0  0],'COlor', [0 0 0 ],'Linewidth',0.5)
ylim([-2 4])

 print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'F_motion_onset_BARQUNAT.pdf'])

%%%%%%%%%%%

 [h,p,ci,stats] = ttest(V1)
% df=20,      0.0189
 [h,p,ci,stats] = ttest(V3)
 % df=8,     0.0986
  [h,p,ci,stats] = ttest(V4)
 % df=15,   0.5290
 
 
     [h,p,ci,stats] =ttest2(V1,V4)
%df=35, tst=2.7,      0.0250
   [h,p,ci,stats] =ttest2(V3,V4)
%df=23, tst=3.8,        0.0412
   [h,p,ci,stats] =ttest2(V1,V3)
%df=28, tst=3.8,          0.6417
 
 
