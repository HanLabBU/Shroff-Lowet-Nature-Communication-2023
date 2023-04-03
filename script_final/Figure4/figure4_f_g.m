
%addpath(genpath('Z:\EricLowet\'))

clear all
pT=[];pT2=[];  nn=0;
cell_id=[];Is_ratio3=[];d_id=[];
allM=[];allM2=[]; 
sw_t=2;
nn2=0;
for gh=[1 2 ]
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
        %%%%%%%%%%%% identification vector, 1=delta, 2=regular/other, 0= low SNR
        dyn_type= [ 0 2 1 1 1 1 1 1 1 1    1 2 1 1 2 1 1 0 1 1   1 0 1 1 0 2 1 2 0 0  0 0];
       
        %%%%%%%%%%%%%%
    else
        %cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Hua-an_Tseng/Data/SomArchon_Striatum/Linear_movement_PinnedBall/good_data_ALL/MSN/')
        cd('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata2\Hua-an Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\MSN\')
       %        cd('Z:\eng_research_handata2\Hua-an_Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\MSN\')
        dyn_type= [ 0 2 0 2 0 2 2 0 2 2    0 1 1 2 1 1 0 1 1 0    0 1 1 1 1 2 2 2 2 0     0 0 0 0];
dyn_type= [2,2,2,2,0,2,2,0,2,2,0,1,1,2,0,1,2,1,1,2,0,1,1,1,1,2,2,2,2,0,2,0  0 0 0]

        if sw_t==1      % delta-identified
           chs=[find(dyn_type==1)];
        elseif sw_t==2     % regular-identified or theta bursting
            chs=  [ find(dyn_type==2) ]
        elseif sw_t==3      % all
            chs=[[find(dyn_type>0)]];
            
        end
        %%%%%%%%%%%% identification vector, 1=delta, 2=regular/other, 0= low SNR
        dyn_type= [ 0 2 0 2 0 2 2 0 2 2    0 1 1 2 1 1 0 1 1 0    0 1 1 1 1 2 2 2 2 0     0 0 0 0];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end


ses= dir('*.mat')

for sesd=[chs]  %[3,4,5,6,7,8,9,10,11,13,14,16,17,19,20,21,23,24,27,29 12 15  26 28  ]] %1:length(ses);%[1:15 18:31] %:length(ses)  %22 exclude 16 17
    sesd
 
    load(ses(sesd).name)

%     

%    
try 
if size(aligned.vectors,1)>0
   FS=1000;
   lfp=[];mk=0;

for in =1:size(aligned.vectors,1)%length(aligned.trace)
    wind=length(aligned.trace{in});%ceil(FS*1);
    clear d A
vv=aligned.lfp{in};
vv=aligned.trace_nosp{in};
            Fn = FS/2;FB=[ 2 4];
               [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                            LFPg= ((filtfilt(B,A,  vv)));
                        vv= LFPg;
                        
  d(1,:)= vv; % LFP
   d(1,:)= locmotion_raw{in}.speed;%
  %d(1,:)= aligned.mvmt{in};
  % d(1,:)= spRES{in}.denoise_trace;%
   vv=  locmotion_raw{in}.speed;  %Vm
             %%%%%%%%%%%%%
                                        %% denoise
                                         Fn = FS/2;FB=[ 2 6];
                                         [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                                         LFPg= ((filtfilt(B,A,  vv)));
                                     vv= LFPg;
    d(3,:) =vv;    
    %   d(1,:)= vv;
   %  A=zeros( 1,length(aligned.mvmt{in}(1:end)));
  %  A(aligned.vectors{in,6})=1;
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
   for nh= 1:Nr
 mk=mk+1;
 lfp.trial{mk}(1:3,:) = zscore(d(1:3,(1  + ((wind)* (nh-1)):   wind * (nh))   ),[],2);
   lfp.trial{mk}(1,:) = (d(1,(1  + ((wind)* (nh-1)):   wind * (nh))   ));
  lfp.trial{mk}(4,:) = (d(4,(1  + ((wind)* (nh-1)):   wind * (nh))   ));
 lfp.trial{mk}(5:6,:) = (d(5:6,(1  + ((wind)* (nh-1)):   wind * (nh))   ));
 lfp.trial{mk}(7,:) = (d(7,(1  + ((wind)* (nh-1)):   wind * (nh))   ));
 lfp.trial{mk}(8,:) = (d(8,(1  + ((wind)* (nh-1)):   wind * (nh))   ));

  lfp.time{mk}= (1:size(lfp.trial{mk},2))./FS;
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
s= lfp.trial{id}(7,:)>0;s=find(s)+0;
ss= lfp.trial{id}(8,:)>0;ss=find(ss);
s2= lfp.trial{id}(2,:)>0;%s=find(s);
motD= lfp.trial{id}(1,:);%s=find(s);
motDA= angle(hilbert((lfp.trial{id}(3,:))));%s=find(s);
 imag_mot= lfp.trial{id}(4,:); 
high_mot= lfp.trial{id}(5,:)>0; low_mot= lfp.trial{id}(6,:)>0; 

  NT=length(s); TL= length(s2);
    
                VV=randperm(TL);VV(VV<wind | VV+wind > TL)=[];
                sh=VV(1:NT);
          swin=100;      rwin=400;  
for id2=1:length(s)
    if s(id2) >wind+  swin  & s(id2) +wind < length(s2)-rwin  %& high_mot(s(id2))>0
        mm=mm+1;
     rr=   find(motDA(s(id2)-  swin:s(id2)+rwin)>-0.15 &motDA(s(id2)-  swin:s(id2)+rwin)<0.15);
     if ~isempty(rr)
     rr=rr(1)-  swin;
     else, rr=0;end;rr=0;
     aM(:,mm)=  fastsmooth(s2([s(id2)-wind:s(id2)+wind]+rr),100,1,1);
     aM2(:,mm)= (motD([s(id2)-wind:s(id2)+wind]+rr));
    end
end


if  0
   figure,plot(motD); hold on,
   plot(motDA)
    hold on,plot(s(id2)+rr,1,'om')
    
end
% for id2=1:length(ss)
    end
    try
      allM=[allM, aM];
    allM2=[allM2, aM2];end


 %  figure,imagesc(nanmean(abs(aM),3)); axis xy
    
    %  figure,imagesc(abs(nanmean(exp(1i.*angle(aM)),3))); axis xy
    if 1 %sum(~isnan(aM(2,101,:)))>20  & sum(~isnan(aM2(2,101,:)))>20
    nn=nn+1;
%allP(:,:,nn)=  nanmean(abs(aM),3);
try
allC(:,nn)=nansum(aM,2);
allC2(:,nn)=nansum(aM2,2);
alltrig(nn)=size(aM,2);
%alltrig2(nn)=size(aM2,2);
if size(aM,2)>=5
    nn2=nn2+1;
    for x1=1:5
DATA(x1).A(:,nn2)=zscore(fastsmooth(aM(80:end-250,x1),60,1,1));
DATA(x1).times= (1:size(aM(80:end-250,x1),1));
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    

    end
       cell_id=[cell_id,gh];
       d_id=[d_id,  dyn_type(sesd) ];
    %%%%%%%%%%%%%%%%%%%%%%
   

   end 
    
 
end
end
end


%cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\git_scripts\fig4\')
savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\git_scripts\fig4\'
%savepath='Z:\EricLowet\git_scripts\fig4\'

pheight=160;

tt=-wind:wind;
tsel=tt<-100;

clear allX
for ind=1:size(allC,2)
   allX(:,ind)= zscore(fastsmooth(allC(50:end-50,ind),60,1,1)); 
end
allX=allX(100:500,:);



 selC= d_id>0;
clear MM
  figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
 M= fastsmooth(nanmean(allC(:,selC)./alltrig(selC),2),50,1,1).*1000;
 Ms= fastsmooth(nanstd(allC(:,selC)./alltrig(selC),[],2)./sqrt(length(find(selC))),50,1,1).*1000;
 M2= fastsmooth(nanmean(allC2(:,selC)./alltrig(selC),2),50,1,1).*.1+0;
 M2s= fastsmooth(nanstd(allC2(:,selC)./alltrig(selC),[],2)./sqrt(length(find(selC))),50,1,1).*.1;
 plot(tt,M,'k','Linewidth',1.4,'COlor',[ 0.2 0.2 0.2]); hold on,
fill_error_area2(tt,M,Ms,[0.5 0.5 0.5]);hold on,
plot(tt,M,'k','Linewidth',1.5,'COlor',[ 0.2 0.2 0.2]); hold on,
fill_error_area2(tt,M2,M2s,[0.5 0.5 0.5]);hold on,
plot(tt,M2,'r','Linewidth',1.5,'COlor',[ 0.3 0.7 0.7]); hold on,
 xlim([-300 300])
 ylim([1 15])
  print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'MOT_motion_onset_ ' num2str(sw_t) '.pdf'])

%  selC=cell_id==2;
% clear MM
%   figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
%  M= fastsmooth(nanmean(allC(:,selC)./alltrig(selC),2),50,1,1).*1000;
%  Ms= fastsmooth(nanstd(allC(:,selC)./alltrig(selC),[],2)./sqrt(length(find(selC))),50,1,1).*1000;
%  M2= fastsmooth(nanmean(allC2(:,selC)./alltrig(selC),2),50,1,1).*5+3;
%  M2s= fastsmooth(nanstd(allC2(:,selC)./alltrig(selC),[],2)./sqrt(length(find(selC))),50,1,1).*5;
%  plot(tt,M,'k','Linewidth',1.4,'COlor',[ 0.2 0.2 0.2]); hold on,
% fill_error_area2(tt,M,Ms,[0.5 0.5 0.5]);hold on,
% plot(tt,M,'k','Linewidth',1.5,'COlor',[ 0.2 0.2 0.2]); hold on,
% fill_error_area2(tt,M2,M2s,[0.5 0.5 0.5]);hold on,
% plot(tt,M2,'r','Linewidth',1.5,'COlor',[ 0.5 0.2 0.2]); hold on,
%  %print(gc