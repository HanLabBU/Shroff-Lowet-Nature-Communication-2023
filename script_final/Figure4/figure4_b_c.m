clear all
%close all
%% Load data
inc=[]; dec=[];
% [files,paths]= uigetfile('*.mat','MultiSelect','on');
low_speed_spike_rate_all=[];
hi_speed_spike_rate_all=[];
low_speed_rate=0;
hi_speed_rate=0;
    select_delta_neurons=1;
for gh=1:2
    
    if gh==1
%cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Hua-an_Tseng/Data/SomArchon_Striatum/Linear_movement_PinnedBall/good_data_ALL/ChAT/')
    cd('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata2\Hua-an Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\ChAT\')
 
     %%%%%%%%%%%% identification vector, 1=delta, 2=regular/other, 0= low SNR
        dyn_type= [ 0,2,1,1,1,1,1,1,1,0,1,2,1,1,1,1,1,0,1,1,1,1,1,1,1,2,1,2,0,2,0,1] ;
   
        %%%%%%%%%%%%%%
    else
%cd('/home/hanlabadmins/eng_handata/eng_research_handata2/Hua-an_Tseng/Data/SomArchon_Striatum/Linear_movement_PinnedBall/good_data_ALL/MSN/')
cd('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata2\Hua-an Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\good_data_ALL\MSN\')
  
%%%%%%%%%%%% identification vector, 1=delta, 2=regular/other, 0= low SNR
        dyn_type= [2,2,2,2,0,2,2,0,2,2,0,1,1,2,0,1,2,1,1,2,0,1,1,1,1,2,2,2,2,0,2,0  0 0 0];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

ses= dir('*.mat');
if  select_delta_neurons
n_files=size(find(dyn_type==1 ),2);
select_idx=find(dyn_type==1 );
else
  n_files=size(find(dyn_type==2 ),2);
select_idx=find(dyn_type==2 );  
end


inc_responsive = zeros(numel(n_files),1);
dec_responsive = zeros(numel(n_files),1);




for sesd=1:n_files
    sesd
load(ses(select_idx(sesd)).name)



a= aligned.vectors;
low_spike_rate=NaN(1,size(aligned.vectors,1));
hi_spike_rate=NaN(1,size(aligned.vectors,1));
%% Observed value 

for trial=1:size(aligned.vectors,1)

spikes=aligned.spike{trial};
hi_start_index=a{trial,3};
hi_stop_index= a{trial,4};
low_start_index=a{trial,5};
low_stop_index= a{trial,6};  


% Low speed spike rate
low_spike_count=0;
low_speed_frame_count=0;
for j=1:numel(low_start_index)
    low_speed_frames=spikes(low_start_index(j):low_stop_index(j));
    low_spike_count=low_spike_count+nansum(low_speed_frames);
    low_speed_frame_count=low_speed_frame_count+numel(low_speed_frames);
end
low_spike_rate(trial)= low_spike_count/low_speed_frame_count;

% High speed spike rate
hi_spike_count=0;
hi_speed_frame_count=0;
for i=1:numel(hi_start_index)
    hi_speed_frames=spikes(hi_start_index(i):hi_stop_index(i));
    hi_spike_count=hi_spike_count+nansum(hi_speed_frames);
    hi_speed_frame_count=hi_speed_frame_count+numel(hi_speed_frames);
end
hi_spike_rate(trial)= hi_spike_count/hi_speed_frame_count;

end
% Average across trials 
avg_low_spike_rate=nanmean(low_spike_rate);
avg_hi_spike_rate=nanmean(hi_spike_rate);
% Difference 
actual_diff= avg_hi_spike_rate-avg_low_spike_rate;


%% Shuffled value
shuffle_diff=zeros(1000,1);
for n=1:1000
% shuffle_low_spike_rate=NaN(1,size(aligned.vectors,1));
% shuffle_hi_spike_rate=NaN(1,size(aligned.vectors,1));shuffle_
shuffle_low_spike_count_all=0;
shuffle_low_speed_frame_tot=0;
shuffle_hi_spike_count_all=0;
shuffle_hi_speed_frame_tot=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shuff_trial=randperm(size(aligned.vectors,1)); %%%% ADDED line !!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for trial=1:size(aligned.vectors,1)
trial=shuff_trial(trial);  %%%% ADDED LINE !!!!!!!!!!!!!!!!!!!!!!
spikes=aligned.spike{trial};
hi_start_index=a{trial,3};
hi_stop_index= a{trial,4};
low_start_index=a{trial,5};
low_stop_index= a{trial,6};  

% if numel(hi_start_index)<=0 || numel(low_start_index)<=0
%     continue
% end
% Low speed spike rate

shuffle_low_spike_count=0;
low_speed_frame_count=0;
for j=1:numel(low_start_index)
    low_speed_frames=randperm(length(spikes),low_stop_index(j)-low_start_index(j));
    shuffle_low_spike_count=shuffle_low_spike_count+nansum(spikes(low_speed_frames));
    low_speed_frame_count=low_speed_frame_count+numel(low_speed_frames);
end
shuffle_low_spike_count_all=shuffle_low_spike_count_all+shuffle_low_spike_count;
shuffle_low_speed_frame_tot=shuffle_low_speed_frame_tot+low_speed_frame_count;

% High speed spike rate
shuffle_hi_spike_count=0;
hi_speed_frame_count=0;
for i=1:numel(hi_start_index)
    hi_speed_frames=randperm(length(spikes),hi_stop_index(i)-hi_start_index(i));
    shuffle_hi_spike_count=shuffle_hi_spike_count+nansum(spikes(hi_speed_frames));
    hi_speed_frame_count=hi_speed_frame_count+numel(hi_speed_frames);
end
shuffle_hi_spike_count_all=shuffle_hi_spike_count_all+shuffle_hi_spike_count;
shuffle_hi_speed_frame_tot=shuffle_hi_speed_frame_tot+hi_speed_frame_count;
end
% Average across trials 
shuffle_avg_low_spike_rate=shuffle_low_spike_count_all/shuffle_low_speed_frame_tot;
shuffle_avg_hi_spike_rate=shuffle_hi_spike_count_all/shuffle_hi_speed_frame_tot;
% Difference 
shuffle_diff(n)= shuffle_avg_hi_spike_rate-shuffle_avg_low_spike_rate;
end
threshold=prctile(shuffle_diff',97.5);
inc_responsive(sesd)= actual_diff>threshold;

if 0%inc_responsive(sesd) ==1
figure(sesd)
hist(shuffle_diff)
line([actual_diff actual_diff], [0 300],'Color','g')
line([prctile(shuffle_diff',97.5) prctile(shuffle_diff',97.5)], [0 300],'Color','r')
end

threshold2=prctile(shuffle_diff',2.5);
dec_responsive(sesd)= actual_diff<threshold2;

low_speed_spike_rate_all =[low_speed_spike_rate_all,avg_low_spike_rate];
hi_speed_spike_rate_all =[hi_speed_spike_rate_all,avg_hi_spike_rate];
%    
end
inc=[inc ,inc_responsive];
dec=[dec ,dec_responsive];
end
 
cd('./Shuffle mat files/Figure 1')

inc_responsive=inc;
dec_responsive=dec;
percent=[low_speed_spike_rate_all;hi_speed_spike_rate_all];
% % Convert to spikes/s
percent=percent'*1000;
figure
boxplot(percent,'notch','on')
hold on
idx=find(inc_responsive==1);
idx2=find(dec_responsive==1);
idx3 = find(inc_responsive==0 & dec_responsive==0);
hold on
plot([percent(idx,1),percent(idx,2)]','.-','color','r','Markersize',10,'MarkerEdgeColor','k')
plot([percent(idx2,1),percent(idx2,2)]','.-','color','b','Markersize',10,'MarkerEdgeColor','k')
plot([percent(idx3,1),percent(idx3,2)]','.-','color',[0.5 0.5 0.5],'Markersize',10,'MarkerEdgeColor','k')
xticks([1 2])
xticklabels({'Low movement','High movement'})
ylabel('Spike rate (Hz)')
%hold off
[h,p,ci,stats] = ttest(percent(:,1),percent(:,2));
 if p<0.05 
     plot([1,2],[max(max(percent(:,1)),max(percent(:,2)))+0.5,max(max(percent(:,1)),max(percent(:,2)))+0.5],'k')
     text(1.5,max(max(percent(:,1)),max(percent(:,2)))+0.7,'*')
 end
ylim([-1 max(max(percent(:,1)),max(percent(:,2)))+1.5])
hold off
% saveas(gcf, 'CHAT_motion_resp_cells_box_with_space.fig')
% saveas(gcf, 'CHAT_motion_resp_cells_box_with_space.jpg')
% saveas(gcf, 'CHAT_motion_resp_cells_box_with_space.pdf')
% %saveas(gcf, 'MSN_motion_resp_cells_box_with_space.eps')


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\git_scripts\fig4\'

colors=[ [0.4 0.8 0];[0.9 0.5 .9]]
pheight=160
  figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters') 

percent=[low_speed_spike_rate_all;hi_speed_spike_rate_all];
% % Convert to spikes/s
percent=percent'*1000;
boxplot(percent,'notch','on','COlor',[ 0.2 0.2 0.2],'Widths',0.6, 'symbol','.k' )
hold on
idx=find(inc_responsive==1);
idx2=find(dec_responsive==1);
idx3 = find(inc_responsive==0 & dec_responsive==0);
hold on
%plot([percent(idx,1),percent(idx,2)]','-','color',[0.4 0.4 0.4],'Markersize',10,'MarkerEdgeColor','k')
%plot([percent(idx2,1),percent(idx2,2)]','-','color',[0.4 0.4 0.4],'Markersize',10,'MarkerEdgeColor','k')
%plot([percent(idx3,1),percent(idx3,2)]','-','color',[0.4 0.4 0.4],'Markersize',10,'MarkerEdgeColor','k')
xticks([1 2])
xticklabels({'L','H'})
%ylabel('Spike rate (Hz)')
%hold off
[h,p,ci,stats] = ttest(percent(:,1),percent(:,2))    ;
%df=30, 0.0026,   , 12,3,16
%df=20,  0.5561,   11,4,6
 if p<0.05 
  %   plot([1,2],[max(max(percent(:,1)),max(percent(:,2)))+0.5,max(max(percent(:,1)),max(percent(:,2)))+0.5],'k')
     text(1.5,max(max(percent(:,1)),max(percent(:,2)))+0.7,'*')
 end
ylim([-1 max(max(percent(:,1)),max(percent(:,2)))+1.5])
boxplot(percent,'notch','on','COlor',[ 0.2 0.2 0.2],'Widths',0.6 , 'symbol','.k')
hold off
h = findobj(gca,'Tag','Box');
for j=1:2
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end



xlim([.4 2.6])
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Fire_effect_low_high_' num2str(select_delta_neurons) '.pdf'])



figure('COlor','w','Position', [ 300 400 120 pheight-10],'Renderer', 'painters')
 violins1 =violinplotSTR(percent(:,1),[1.3 ],'ViolinColor', [0.9 0.5 .9])
hold on,
 violins2 =violinplotSTR(percent(:,2),[ 1.9],'ViolinColor', [0.4 0.8 0])
XX=violins1.ScatterPlot.XData;XX2=violins2.ScatterPlot.XData;
line([ XX'  XX2']',[percent(:,1),percent(:,2)]','color',[ 0.7 0.7 0.7 ])
line([ 0.8 2.4], [ 0  0],'COlor', [0 0 0 ],'Linewidth',0.5)
xlim([.8 2.4])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Fire_effect_low_high_' num2str(select_delta_neurons) '.pdf'])
