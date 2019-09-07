b=0;
c=0;
d=0;
for i=10:10:100
    c=c+1;
    for j=0.3:0.1:0.9
        d=d+1;
        for y=0.01:0.01:0.2
            b=b+1;
                optimproblem.options.PopulationSize = i;
                optimproblem.options.CrossoverFraction = j;
                optimproblem.options.MutationFcn{2} = y;
                t = ga(optimproblem);
            A(b,d,c)=sum(t);
           % print(a);
        end
        b=0;
    end
    d=0;
end
%%
for i=1:180
    x(i)=rand()*1000; % a random real number between 0-1 * norma . like a pointer
    y(i)=rand()*1000;
end

rgNodes = transpose([x ; y]);

%[football,rggNodeDegree]=rgg(rgNodes,115,250);

%[lesmis,rggNodeDegree]=rgg(rgNodes,77,250);
%[dolphins,rggNodeDegree]=rgg(rgNodes,62,250);
%Lets make the graphs matrix
lesmis =importgml('lesmis.gml');
football =importgml('football.gml');
dolphins =importgml('dolphins.gml');

if isdirected(football)
    B=undir(football,size(football,1));
    football=B;
    clearvars B;
end
if isdirected(lesmis)
    B=undir(football,size(lesmis,1));
    lesmis=B;
    clearvars B;
end
if isdirected(dolphins)
    B=undir(football,size(dolphins,1));
    dolphins=B;
    clearvars B;
end
%%
maxi=0;
for i=0.7:0.1:0.9
    for j=0.1:0.1:0.2
        for k=1:3
            res=cellcluster(genetic(lesmis,77,i,j,k),77);
            if QFModul(res,lesmis)>maxi
                maxi=QFModul(res,lesmis);
                maxcluster=res;
            end
        end
    end
end
    
%%            
%res=genetic(lesmis,77,0.7,0.2,2);
figure;
PlotGraph(lesmis,maxcluster);
title('Genetic algorithm for lesmis');
saveas(gcf,'Genetic_lesmis.png');
%% football
clearvars res;
maxif=0;
for i=0.7:0.1:0.9
    for j=0.1:0.1:0.2
        for k=1:3
            res=cellcluster(genetic(football,115,i,j,k),115);
            if QFModul(res,football)>maxif
                maxif=QFModul(res,football);
                maxclusterf=res;
                p_c_f=i;
                p_m_f=j;
                elitism_f=k;
            end
        end
    end
end
    
            
%res=genetic(lesmis,77,0.7,0.2,2);
figure;
PlotGraph(football,maxclusterf);
title('Genetic algorithm for football');
saveas(gcf,'Genetic_foot.png');
%% dolphins
clearvars res;
%maxid=0;
parfor i=1:3
    maxid(i)=0;   
    for j=0.1:0.1:0.2
        for k=1:3
            res=cellcluster(genetic(dolphins,62,0.6 +i/10 ,j,k),62);
            if QFModul(res,dolphins)>maxid(i)
                maxid(i)=QFModul(res,dolphins);
                maxclusterd{i}=res;
                p_m_d(i)=j;
                elitism_d(i)=k;
            end
        end
    end
end
    
[~,ind]=max(maxid);
%res=genetic(lesmis,77,0.7,0.2,2);
figure;
PlotGraph(dolphins,maxclusterd{ind});
title('Genetic algorithm for dolphins');
saveas(gcf,'Genetic_dol.png');
%football paraller
% football
clearvars res;
%maxid=0;
%%
parfor i=1:3
    maxif1(i)=0;   
    for j=0.1:0.1:0.2
        for k=1:3
            res=cellcluster(genetic(football,115,0.6 +i/10 ,j,k),115);
            if QFModul(res,football)>maxif1(i)
                maxif1(i)=QFModul(res,football);
                maxclusterf1{i}=res;
                p_m_f1(i)=j;
                elitism_f1(i)=k;
            end
        end
    end
end
[~,ind]=max(maxif1);
%res=genetic(lesmis,77,0.7,0.2,2);
figure;
PlotGraph(football,maxclusterf1{ind});
title('Genetic algorithm for football');
saveas(gcf,'Genetic_football.png');
%%
clearvars res;
%maxid=0;
parfor i=1:3
    maxil(i)=0;   
    for j=0.1:0.1:0.2
        for k=1:3
            res=cellcluster(genetic(lesmis,77,0.6 +i/10 ,j,k),77);
            if QFModul(res,lesmis)>maxil(i)
                maxil(i)=QFModul(res,lesmis);
                maxclusterl{i}=res;
                p_m_l(i)=j;
                elitism_l(i)=k;
            end
        end
    end
end
    
[~,ind]=max(maxil);
%res=genetic(lesmis,77,0.7,0.2,2);
figure;
PlotGraph(lesmis,maxclusterl{ind});
title('Genetic algorithm for lesmis');
saveas(gcf,'Genetic_lesmis2.png');
