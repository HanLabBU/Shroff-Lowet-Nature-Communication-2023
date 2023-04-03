savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\git_scripts\fig2\'

pheight=160;


      dyn_type= [ 0,2,1,1,1,1,1,1,1,0,1,2,1,1,1,1,1,0,1,1,1,1,1,1,1,2,1,2,0,2,0,1]
      

      D1=find(dyn_type==1);
            D2=find(dyn_type==2);
              Excluded=find(dyn_type==0); % 5 neurons
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 p1=   length(D1);%length(find(allProbs(2,:)./mean(allProbs(1,:),1)>0.15));
p2=   length(D2);%length(find(allProbs(2,:)./mean(allProbs(1,:),1)<=0.15));
 figure('COlor','w','Position',[ 300 400 250 pheight],'Renderer', 'painters'),
h=pie(  [ p1 p2],[ 0 1] )
 patchHand = findobj(h, 'Type', 'Patch'); 
set(patchHand, {'FaceColor'}, mat2cell([0.9 0  0;0  0 0.9], ones(2,1), 3))
patchHand(2).FaceAlpha = 0.7;
patchHand(1).FaceAlpha = 0.7;
 print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'pie_chart_CHAT.pdf'])


 %%%%%%%%
  dyn_type= [2,2,2,2,0,2,2,0,2,2,0,1,1,2,0,1,2,1,1,2,0,1,1,1,1,2,2,2,2,0,2,0  0 0 0]
      D1=find(dyn_type==1);
            D2=find(dyn_type==2);
Excluded=find(dyn_type==0); % 5 neurons
 
  p1=   length(D1);%length(find(allProbs(2,:)./mean(allProbs(1,:),1)>0.15));
p2=   length(D2);%length(find(allProbs(2,:)./mean(allProbs(1,:),1)<=0.15));
 figure('COlor','w','Position',[ 300 400 250 pheight],'Renderer', 'painters'),
h=pie(  [ p1 p2],[ 0 1] )
 patchHand = findobj(h, 'Type', 'Patch'); 
set(patchHand, {'FaceColor'}, mat2cell([0.9 0  0;0  0 0.9], ones(2,1), 3))
patchHand(2).FaceAlpha = 0.7;
patchHand(1).FaceAlpha = 0.7;
  print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'pie_chart_SYN.pdf'])

 
