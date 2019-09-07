clc;clear all;close all;

addpath(genpath('./lab2code'));
addpath(genpath('./epidemics_code'));

lesmis =importgml('lesmis.gml');
football =importgml('football.gml');
dolphins =importgml('dolphins.gml');

% check if graphs are directed or not and if it is make it undirected like
% Lab2

if isdirected(lesmis)
    temp=undir(lesmis,size(lesmis,1));
    lesmis=temp;
    clearvars temp;
end

if isdirected(football)
    temp=undir(football,size(football,1));
    football=temp;
    clearvars temp;
end

if isdirected(dolphins)
    temp=undir(dolphins,size(dolphins,1));
    dolphins=temp;
    clearvars B;
end

%% LESMIS

lesmis_max=0;
for i=0.7:0.1:0.9
    for j=0.1:0.1:0.2
        for k=1:3
            clust=cellcluster(MyGeneticAlgo(lesmis,i,j,k,30,300,5),77);
            if QFModul(clust,lesmis)>lesmis_max
                lesmis_max=QFModul(clust,lesmis);
                lesmis_maxcluster=clust;
            end
        end
    end
end        

PlotGraph(lesmis,lesmis_maxcluster);
javaframe = get(handle(gcf),'JavaFrame');
set(javaframe,'Maximized',1); 
label = ['Best Results for Lesmis with QFModule and Genetic Algorithm :' num2str(lesmis_max)];
title(label);
saveas(gcf,'Lesmis.png');

%% FOOTBALL

football_max=0;
for i=0.7:0.1:0.9
    for j=0.1:0.1:0.2
        for k=1:3
            clust=cellcluster(MyGeneticAlgo(football,i,j,k,30,300,5),115);
            if QFModul(clust,football)>football_max
                football_max=QFModul(clust,football);
                football_maxcluster=clust;
            end
        end
    end
end           

PlotGraph(football,football_maxcluster);
javaframe = get(handle(gcf),'JavaFrame');
set(javaframe,'Maximized',1); 
label = ['Best Results for Football with QFModule and Genetic Algorithm :' num2str(football_max)];
title(label);
saveas(gcf,'Football.png');

%% DOLPHINS

dolphins_max=0;
for i=0.7:0.1:0.9
    for j=0.1:0.1:0.2
        for k=1:3
            clust=cellcluster(MyGeneticAlgo(dolphins,i,j,k,30,300,5),62);
            if QFModul(clust,dolphins)>dolphins_max
                dolphins_max=QFModul(clust,dolphins);
                dolphins_maxcluster=clust;
            end
        end
    end
end

PlotGraph(dolphins,dolphins_maxcluster);
javaframe = get(handle(gcf),'JavaFrame'); 
set(javaframe,'Maximized',1); 
label = ['Best Results for Dolphins with QFModule and Genetic Algorithm :' num2str(dolphins_max)];
title(label);
saveas(gcf,'Dolphins.png');

