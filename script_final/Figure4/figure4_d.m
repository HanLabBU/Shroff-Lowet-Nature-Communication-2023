%% Spike trigerred power plots

clear all
%close all
%% Load files =[];
Is_ratio=[];Is_ratio2=[];P_ratio=[];cell_id=[];Is_ratio3=[];d_id=[];
  nn=0;nn1=0;
%%%%%%%%%%%%%%%
sw_t=1; % switch cell seection, 1==delta, 2==regular, 3=both
%%%%%%%%%%%%%%%%%
for gh=[1 2 ]  % 1= CHAT, 2=SYN, 
 
   
    
    
    if gh==1
        cd('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata2\Hua-an Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\ChAT\')
                      dyn_type=   [ 0,2,1,1,1,1,1,1,1,0,1,2,1,1,1,1,1,0,1,1,1,1,1,1,1,2,1,2,0,2,0,1] 
                        

   chs=find(dyn_type==1)    ;    
   
        %%%%%%%%%%%%%%
    else
        %cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Hua-an_Tseng/Data/SomArchon_Striatum/Linear_movement_PinnedBall/good_data_ALL/MSN/')
        cd('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata2\Hua-an Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\MSN\')
        dyn_type= [2,2,2,2,0,2,2,0,2,2,0,1,1,2,0,1,2,1,1,2,0,1,1,1,1,2,2,2,2,0,2,0  0 0 0]
      

          chs=find(dyn_type==1);
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    
    
    ses= dir('*.mat')
    
    for sesd=[ chs]
        sesd
        load(ses(sesd).name) %LOAD
        
        try
            snr_thres= 0;%
            if length(aligned.trace)>0% Actual - low movement
                
                FS=1000;
                lfp=[];mk=0;
                SN=[];
                
                %% Put in a matrix d - LFP, spikes, Vm, High motion and Low motion
                
                for in =1:length(aligned.trace)
                    windF= length(aligned.trace{in});ceil(FS*1);
                    clear d A
                    d(1,:)= aligned.lfp{in}(:,1:end); % LFP
                    A=aligned.spike{in}(:,1:end);A(isnan(A))=0;
                    A= spRES{in}.roaster;
                    A=aligned.spike{in}(:,1:end);
                    A(isnan(A))=0;
                    d(2,:)=A;  %% SPIKES
                   % d(3,:)= aligned.trace_nosp{in}(:,1:end);  %Vm
                        LFPtr=   aligned.trace_nosp{in}(:,1:end);
                                  Fn = 1000/2;FB=[ 83 86];
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
                  %  A(aligned.vectors{in,8})=1;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    d(6,:)=  aligned.imgmotion{in};  % combined X-Y image motion
                    d(6,[1 2 ])=d(6,3); % remove edge effect
                    d(6,end:-1:end)=d(6,end-2); % remove edge effect
                    
                    d(6,:)= [ 0 fastsmooth(abs(hilbert(diff(fastsmooth(d(6,:),10,1,1)))),100,1,1)];
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    %% Make it suitable for fieldtrip by zscoring
                    if mean(spRES{in}.spike_snr{:}) >0  & length([spRES{in}.spike_snr{:}])>0
                        for nh= 1:Nr
                           % if  sum((d(4,(1  + ((windF)* (nh-1)):   windF * (nh))   ))>0)<100
                            mk=mk+1;
                            lfp.trial{mk}(1:3,:) = zscore(d(1:3,(1  + ((windF)* (nh-1)):   windF * (nh))   ),[],2);
                            lfp.trial{mk}(4:5,:) = (d(4:5,(1  + ((windF)* (nh-1)):   windF * (nh))   )); % locomotion
                            lfp.trial{mk}(6,:) = (d(6,(1  + ((windF)* (nh-1)):   windF * (nh))   ));% image motion
                            lfp.time{mk}= (1:size(lfp.trial{mk},2))./FS;end
                       % end
                    end
                end
                n=0;
                for ind=1:size(lfp.trial{1},1)
                    n=n+1;
                    lfp.label(n)= { [ 'CH' num2str(ind)]};
                end
                
                if length(lfp.trial) >1 %  & mean(SN) >  5 & length(SN)>5
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    ff=[];
                    for i=1:length(lfp.trial)
                     ff=[ff, sum(lfp.trial{i}(2,:)>0)];
                    end
                    Is_ratio3=[Is_ratio3; mean(ff)];
                    
 INT1=[];INT2=[];                    
for in=1:length(spRES)
    high_mot= lfp.trial{in}(4,1:end-1)>0; low_mot= lfp.trial{in}(5,1:end-1)>0; 
A=spRES{in}.roaster;A(isnan(A))=0;
v=find(A(high_mot>0));vx=v(2:end)-v(1:end-1)  ;
v=find(A(low_mot>0));vx2=v(2:end)-v(1:end-1)  ;
INT1=[INT1, vx ];
INT2=[INT2, vx2 ];
end

if 1 % length(INT1)>10 & length(INT1)>10
nn=nn+1;

                [n1 n2]=hist(INT1,[2:2:180]);
    allFIR_IB(nn,1) = 1000/(sum(n2(1:end-1).*n1(1:end-1))./sum(n1(1:end-1)))
                   [n1 n2]=hist(INT2,[2:2:180]);
    allFIR_IB(nn,2) = 1000/(sum(n2(1:end-1).*n1(1:end-1))./sum(n1(1:end-1)))
  
                  [n1 n2]=hist(INT1,[400:4:1800]);
    allFIR_BB(nn,1) = 1000/(sum(n2(1:end-1).*n1(1:end-1))./sum(n1(1:end-1)))
                  [n1 n2]=hist(INT2,[400:4:1800]);
    allFIR_BB(nn,2) = 1000/(sum(n2(1:end-1).*n1(1:end-1))./sum(n1(1:end-1)))
                    cell_id=[cell_id,gh];
                    d_id=[d_id,  dyn_type(sesd) ];
end
                end
            end
        end
    end
end


cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\git_scripts\fig2\')
savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\git_scripts\fig2\'

pheight=160;


colors = [ [0.4 0.8 0];[0.9 0.5 .9] ];;
V1=   allFIR_BB(:,2);V2=   allFIR_BB(:,1);
%   figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters') 
% %b1=bar(1,nanmean(V1),'Facecolor',[ 0.4 0.4 0.4]);hold on,
% h=boxplot( [V1 ;V2]  ,[ ones(length(V1),1)+0 ;ones(length(V2),1).*2],  'notch','on', 'Widths',0.6 , 'colors',[ 0.2 0.2 0.2], 'symbol','.k')
% h = findobj(gca,'Tag','Box');
% for j=1:length(h)
%     patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
% end
% xlim([.4 2.6])
% ylim([1.8 2.7])
% print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'FR_INTRA_barquantLOWHIG.pdf'])

figure('COlor','w','Position', [ 300 400 120 pheight-10],'Renderer', 'painters')
 violins1 =violinplotSTR(V1,[1.3 ],'ViolinColor', [0.9 0.5 .9])
hold on,
 violins2 =violinplotSTR(V2,[ 1.9],'ViolinColor', [0.4 0.8 0])
line([ 0.8 2.4], [ 0  0],'COlor', [0 0 0 ],'Linewidth',0.5)
XX=violins1.ScatterPlot.XData;XX2=violins2.ScatterPlot.XData;
line([ XX'  XX2']',[V1 ,V2]','color',[ 0.7 0.7 0.7 ])
%ylim([-0.1 0.2])
xlim([.8 2.4])
ylim([1.6 2.8])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'FR_INTRA_barquantLOWHIG.pdf'])

%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'FR_INTRA_barquant.pdf'])

[h,p,ci,stats] = ttest(V1,V2)
%df=30, 0.005

colors = [ [0.4 0.8 0];[0.9 0.5 .9] ];;
V1=   allFIR_IB(:,2);V2=   allFIR_IB(:,1);
%   figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters') 
% %b1=bar(1,nanmean(V1),'Facecolor',[ 0.4 0.4 0.4]);hold on,
% h=boxplot( [V1 ;V2]  ,[ ones(length(V1),1)+0 ;ones(length(V2),1).*2],  'notch','on', 'Widths',0.6 , 'colors',[ 0.2 0.2 0.2], 'symbol','.k')
% h = findobj(gca,'Tag','Box');
% for j=1:length(h)
%     patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
% end
% xlim([.4 2.6])
% ylim([10 60])
% print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'FR_INTER_barquantLOWHIG.pdf'])


figure('COlor','w','Position', [ 300 400 120 pheight-10],'Renderer', 'painters')
 violins1 =violinplotSTR(V1,[1.3 ],'ViolinColor', [0.9 0.5 .9])
hold on,
 violins2 =violinplotSTR(V2,[ 1.9],'ViolinColor', [0.4 0.8 0])
line([ 0.8 2.4], [ 0  0],'COlor', [0 0 0 ],'Linewidth',0.5)
XX=violins1.ScatterPlot.XData;XX2=violins2.ScatterPlot.XData;
line([ XX'  XX2']',[V1 ,V2]','color',[ 0.7 0.7 0.7 ])
%ylim([-0.1 0.2])
ylim([5 60]);%
xlim([ 0.8 2.4])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'FR_INTER_barquantLOWHIG.pdf'])




%set(b1,'FaceAlpha',0.4)
%errorbar([1   ],nanmean(V1,1), nanstd(V1)./sqrt(size(V1,1)),'.k');
%ylim([0.15 0.3])
%   6.5548e-04
%, df=42
 %print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'FR_rate_delta_reg_barquant.pdf'])

[h,p,ci,stats] = ttest(V1,V2)
%df=30,  1.7127e-05

