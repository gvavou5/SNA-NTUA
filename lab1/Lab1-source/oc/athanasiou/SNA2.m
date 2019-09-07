clear all;close all;
%% erwthma A
%%%%%sort kai sygkrish apotelesmatwn
[x,y]= getNodeCoordinates(140);
XY = [x y];
REG =smallw(140,2,0);
RGER = erdrey(140,750);
C1=1000*rand([140,2]);
[RGG , rRGG] = FindTopology(C1 , 140 , 250);
SFBA = pref(140, 4);
SWWS =smallw(140,2,0.3);
%  REG
REGcent =brandesBetwCentr(REG); %%% vriskei to betweeness opws sthn prwth askhsh
[sor1 in1] =  sort(REGcent,'descend');
Regc=[sor1 in1] ;
for i = 1:140
  V =find(REG(i,:) ==1);
  S = [i V];
  Adj_temp = subgraph(REG,S);
 Cent =brandesBetwCentr(Adj_temp);
  Ego_centREG(i) = Cent(1);
 end
  [sorted1 is1] = sort(Ego_centREG,'descend');
egoRegc=[sorted1; is1] ;
egoRegc=egoRegc';
 % RGER
 RGERcent =brandesBetwCentr(RGER);
[sor2 in2] =  sort(RGERcent,'descend');
rgerc=[sor2 in2]; 
for i = 1:140
  V =find(RGER(i,:) ==1);
  S = [i V];
  Adj_temp = subgraph(RGER,S);
Cent =brandesBetwCentr(Adj_temp);
  Ego_centRGER(i) = Cent(1);
 end
   [sorted2 is2] = sort(Ego_centRGER,'descend');
egorger=[sorted2 ;is2];
egorger=egorger';
 % RGG
 RGGcent =brandesBetwCentr(RGG);
[sor3 in3] =  sort(RGGcent,'descend');
rggc=[sor3 in3];
 for i = 1:140
  V =find(RGG(i,:) ==1);
  S = [i V];
  Adj_temp = subgraph(RGG,S);
Cent =brandesBetwCentr(Adj_temp);
Ego_centRGG(i) = Cent(1);
 end
[sorted3 is3] = sort(Ego_centRGG,'descend');
egorgg=[sorted3; is3]; 
egorgg=egorgg'; 
% SFBA
 SFBAcent =brandesBetwCentr(SFBA);
[sor4 in4] =  sort(SFBAcent,'descend');
sfbace=[sor4 in4]
 for i = 1:140
  V =find(SFBA(i,:) ==1);
  S = [i V];
  Adj_temp = subgraph(SFBA,S);
Cent =brandesBetwCentr(Adj_temp);
  Ego_centSFBA(i) = Cent(1);
 end
    [sorted4 is4] = sort(Ego_centSFBA,'descend');
 egosf=[sorted4 ;is4];
egosf=egosf';
 % SWWS
 SWWScent =brandesBetwCentr(SWWS);
[sor5 in5] =  sort(SWWScent,'descend');
swcent=[sor5 in5]; 

for i = 1:140
  V =find(SWWS(i,:) ==1);
  S = [i V];
  Adj_temp = subgraph(SWWS,S);
Cent =brandesBetwCentr(Adj_temp);
  Ego_centSWWS(i) = Cent(1);
 end
    [sorted5 is5] = sort(Ego_centSWWS,'descend');
 egosw=[sorted5 ;is5];
 egosw=egosw';
 %% erwthma B
 G = cell(3,1);
 G(1) =mat2cell( importgml('football.gml'));
 G(2) = mat2cell(importgml('lesmis.gml'));
 G(3) = mat2cell(importgml('dolphins.gml'));
 DEG = cell(3,1);
 Clust = cell(3,1);
 Ego_cent = cell(3,1);
 for i = 1:3 
     A = cell2mat(G(i));
     %%%%%%% metatroph se mh kateuthynomeno
     if isdirected(A)
     temp1 = (A' == 1);
     temp2 = (A == 0);
     temp3 = temp1 + temp2;
     temp4 = (temp3 == 2);
     A = A + temp4;
     G(i) =mat2cell(A);
     end
     %%%%%%% vathmos
     TEMP =degrees(A);
     max(size(TEMP))
     DEG(i) = mat2cell(TEMP);
     figure;
     hist(TEMP)
     meanDEG(i) = mean(TEMP);
     %%%%%%% syntelesths omadopoihshs
     [CC1,CC2,CC]= clust_coeff(A); 
     Clust(i) =mat2cell( CC);
     meanClust(i) = mean(CC);
     [supremum,cumdist]=cumulativecentrality(CC,max(size(TEMP)));
     figure;
     plot(cumdist)
     %%%%%%% egokentrikothta
     L = length(TEMP);
     Ego = zeros(L,1);
     for j = 1: L
  v =find(A(j,:) ==1);
  s = [j v];
  adj_temp = subgraph(A,s);
  Atemp = (adj_temp^2).*(ones(size(adj_temp)) - adj_temp);
  U = triu(Atemp)-diag(diag(Atemp));    
  Ego(j) = sum(1./U(U~=0));
     end
    [supremum,cumdist]=cumulativecentrality(Ego,max(size(TEMP)));
     figure;
     plot(cumdist)
 Ego_cent(i) = mat2cell(Ego);        
  MeanEgo(i) = mean(Ego);
 end

%% erwthma C


 % FOOTBALL
Vfootball1=GCSpectralClust2(G{1},6,10);
h = PlotGraph(G{1},Vfootball1(:,6));
Q1(1) = QFModul(Vfootball1(:,6),G{1});
Vfootball2= GCModulMax1(G{1});
h = PlotGraph(G{1},Vfootball2);
Q2(1) = QFModul(Vfootball2,G{1});

% LESMIS
Vlesmis1=GCSpectralClust2(G{2},6,10);
h = PlotGraph(G{2},Vlesmis1(:,6));
Q1(2) = QFModul(Vlesmis1(:,6),G{2});
Vlesmis2= GCModulMax1(G{2});
h = PlotGraph(G{2},Vlesmis2);
Q2(2) = QFModul(Vlesmis2,G{2});

% DOLPHINS
Vdol1=GCSpectralClust2(G{3},6,10);
h = PlotGraph(G{3},Vdol1(:,6));
Q1(3) = QFModul(Vdol1(:,6),G{3});
Vdol2= GCModulMax1(G{3});
h = PlotGraph(G{3},Vdol2);
Q2(3) = QFModul(Vdol2,G{3});

% REG
VREG1=GCSpectralClust2(REG,6,10);
h = PlotGraph(REG,VREG1(:,6));
Q1(4) = QFModul(VREG1(:,6),REG);
VREG2= GCModulMax1(REG);
h = PlotGraph(REG,VREG2);
Q2(4) = QFModul(VREG2,REG);


% RGER
VRGER1=GCSpectralClust2(RGER,6,10);
h = PlotGraph(RGER,VRGER1(:,6));
Q1(5) = QFModul(VRGER1(:,6),RGER);
VRGER2= GCModulMax1(RGER);
h = PlotGraph(RGER,VRGER2);
Q2(5) = QFModul(VRGER2,RGER);

% RGG
VRGG1=GCSpectralClust2(RGG,6,10);
h = PlotGraph(RGG,VRGG1(:,6));
Q1(6) = QFModul(VRGG1(:,6),RGG);
VRGG2= GCModulMax1(RGG);
h = PlotGraph(RGG,VRGG2);
Q2(6) = QFModul(VRGG2,RGG);

% SFBA
VSFBA1=GCSpectralClust2(SFBA,6,10);
h = PlotGraph(SFBA,VSFBA1(:,6));
Q1(7) = QFModul(VSFBA1(:,6),SFBA);
VSFBA2= GCModulMax1(SFBA);
h = PlotGraph(SFBA,VSFBA2);
Q2(7) = QFModul(VSFBA2,SFBA);

% SWWS
VSWWS1=GCSpectralClust2(SWWS,6,10);
h = PlotGraph(SWWS,VSWWS1(:,6));
Q1(8) = QFModul(VSWWS1(:,6),SWWS);
VSWWS2= GCModulMax1(SWWS);
h = PlotGraph(SWWS,VSWWS2);
Q2(8) = QFModul(VSWWS2,SWWS);

%% NEWMANGIRVAN (argei paaaara poly)
%%%% den kserw pws na plotarw tis koinothtes edw!!!!!!!!!

[module_hist1,q1]=newman_comm_fast(G{1});

h=PlotGraph(G{1},q1);
[module_hist2,q2]=newman_comm_fast(G{2});
 h=PlotGraph(G{2},VSWWS1(:,6));
[module_hist3,q3]=newman_comm_fast(G{3});
 h=PlotGraph(G{3},VSWWS1(:,6));
[module_hist4,q4]=newman_comm_fast(REG);
 h=PlotGraph(REG,VSWWS1(:,6));
[module_hist5,q5]=newman_comm_fast(RGER);
 h=PlotGraph(RGER,VSWWS1(:,6));
[module_hist6,q6]=newman_comm_fast(RGG);
 h=PlotGraph(RGG,VSWWS1(:,6));
[module_hist7,q7]=newman_comm_fast(SFBA);
 h=PlotGraph(SFBA,VSWWS1(:,6));
[module_hist8,q8]=newman_comm_fast(SWWS);
 h=PlotGraph(SWWS,VSWWS1(:,6));


%{
[modules1,module_hist1,q1] = newmangirvan(G{1} , 6);
[modules2,module_hist2,q2] = newmangirvan(G{2} , 6);
[modules3,module_hist3,q3] = newmangirvan(G{3} , 6);
[modules4,module_hist4,q4] = newmangirvan(REG , 6);
[modules5,module_hist5,q5] = newmangirvan(RGER , 6);
[modules6,module_hist6,q6] = newmangirvan(RGG , 6);
[modules7,module_hist7,q7] = newmangirvan(SFBA , 6);
[modules8,module_hist8,q8] = newmangirvan(SWWS , 6);

%}
save('sna.mat');