clear all;
ses= dir('*.mat');
%% Example data: Striatal cholinergic interneuron with delta rhythmic Vm and spiking and synapsin-label striatal neuron with regular spiking
%% Plots example data set, compute ISI histogram


%% cd to folder with the example data

exampleNeuron=2; % 1= ChI delta, 2= Syn, Non-delta


load(ses(exampleNeuron).name) %LOAD

FS=1000; % Sampling rate
data=[];mk=0;


for in =1:length(aligned.trace)
    clear d 
    d(1,:)= aligned.trace{in};
    A=aligned.spike{in}(:,1:end);A(isnan(A))=0;
    d(2,:)=A;  %% SPIKES
    d(3,:)=  aligned.trace_nosp{in}(:,1:end);  %Vm
    
    %% GET LOCMOMOTION INFO
    A=zeros( 1,length(aligned.mvmt{in}(1:end)));
    HM= aligned.f_mvmt{in}(1:end);
    HM(HM>=50)=1000;
    HM(HM<50)=-1000;
    d(4,:)=HM;
    d(5,:)=HM.*-1;
    traceNoise(in)= spRES{in}.trace_noise;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    d(6,:)=  aligned.imgmotion{in};  % combined X-Y image motion
    d(6,[1 2 ])=d(6,3); % remove edge effect
    d(6,end:-1:end)=d(6,end-2); % remove edge effect
    d(6,:)= [ 0 fastsmooth(abs(hilbert(diff(fastsmooth(d(6,:),10,1,1)))),100,1,1)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Data structure
            mk=mk+1;
            data.trial{mk}(1:3,:) = d(1:3,: );
            data.trial{mk}(4:5,:) = d(4:5,:  ); % locomotion
            data.trial{mk}(6,:) = d(6,: );% image motion
            data.time{mk}= (1:size(data.trial{mk},2))./FS;
  

end
n=0;
for ind=1:size(data.trial{1},1)
    n=n+1;
    data.label(n)= { [ 'CH' num2str(ind)]};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INT1=[];INT2=[]; % Compute ISI intervals n , n-1
for in=1:size(data.trial,2)   
    A=data.trial{in}(2,:);    
    v=find(A);vx= v(3:end)-v(2:end-1);vx2=v(2:end)-v(1:end-1)  ;
    INT1=[INT1, vx(1:end) ];
    INT2=[INT2, vx2(1:end-1)  ];
    
end


%% Plot first trial
trial_sel=1;
figure('COlor','w','Position', [ 300 200 500 160],'Renderer', 'painters') ;
plot(data.time{trial_sel},data.trial{trial_sel}(1,:)./traceNoise(trial_sel),'k')
axis tight
xlabel('Time ms')
ylabel('SBR')
print(gcf, '-dpdf' , '-r300' ,'-painters', [ 'PLOT_TRACE_exampe_neuron_' num2str(exampleNeuron) '_.pdf'])

%% Plot ISI return map
figure('COlor','w','Position', [ 300 200 240 190],'Renderer', 'painters') ;
scatter1 = scatter(INT1(1:1:end),INT2(1:1:end),8,'MarkerFaceColor','k','MarkerEdgeColor','k');
scatter1.MarkerFaceAlpha = .7;
scatter1.MarkerEdgeAlpha = .1;
set(gca,'Xscale','log','Yscale','log')
axis tight;axis([ 5 1500 5 1500])
line([ 5 1500],[5 1500],'COlor',[ 0.8 0 0],'Linewidth',1)
box on;xlabel('ISI n') ;ylabel(' ISI n-1')
title('ISI return map')
print(gcf, '-dpdf' , '-r300' ,'-painters', [ 'PLOT_ISI_RETURN_exampe_neuron_' num2str(exampleNeuron) '_.pdf'])
