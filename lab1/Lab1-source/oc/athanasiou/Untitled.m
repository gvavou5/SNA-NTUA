%% 
%%Erwthma A)
%density // peirazw
Lattice=smallw(140,2,0);
RandomErdos=erdrey(140,750);
for i=1:140
    x(i)=rand()*1000;
    y(i)=rand()*1000;
end

rgNodes140 = transpose([x ; y]);
clearvars x y;
RandomGeometrical=rgg(rgNodes140,140,i);
ScaleFree=pref(140,4);
Smallworld=smallw(140,2,0.3);
%get node coordinates kalo gia smallw lattice 
%gia ta alla
edges=0;
for i=1:140
    for j=1:140
        if Lattice(i,j)==1
            edges=edges+1;
        end
    end
end
[X,Y]=getNodeCoordinates(140);
X=X';
Y=Y';
field1 = 'Adj';  value1 = Lattice;
field2 = 'x';  value2 = X;
field3 = 'y';  value3 =Y ;
field4 = 'nv';  value4 = 140;
field5 = 'ne';  value5 = edges;
s = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5);
plotGraphBasic(s,4,1);
edges=0;
for i=1:140
    for j=1:140
        if RandomErdos(i,j)==1
            edges=edges+1;
        end
    end
end
for i=1:140
     X(i)=1000*rand;
     Y(i)=1000*rand;
end;
field1 = 'Adj';  value1 = RandomErdos;
field2 = 'x';  value2 = X;
field3 = 'y';  value3 =Y ;
field4 = 'nv';  value4 = 140;
field5 = 'ne';  value5 = edges;
s = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5);
plotGraphBasic(s,4,1);edges=0;
for i=1:140
    for j=1:140
        if RandomGeometrical(i,j)==1
            edges=edges+1;
        end
    end
end
[X,Y]=getNodeCoordinates(140);
X=X';
Y=Y';
field1 = 'Adj';  value1 = RandomGeometrical;
field2 = 'x';  value2 = X;
field3 = 'y';  value3 =Y ;
field4 = 'nv';  value4 = 140;
field5 = 'ne';  value5 = edges;
s = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5);
plotGraphBasic(s,4,1);
edges=0;
for i=1:140
    for j=1:140
        if ScaleFree(i,j)==1
            edges=edges+1;
        end
    end
end
[X,Y]=getNodeCoordinates(140);
X=X';
Y=Y';
field1 = 'Adj';  value1 =ScaleFree;
field2 = 'x';  value2 = X;
field3 = 'y';  value3 =Y ;
field4 = 'nv';  value4 = 140;
field5 = 'ne';  value5 = edges;
s = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5);
plotGraphBasic(s,4,1);edges=0;
for i=1:140
    for j=1:140
        if Smallworld(i,j)==1
            edges=edges+1;
        end
    end
end
[X,Y]=getNodeCoordinates(140);
X=X';
Y=Y';
field1 = 'Adj';  value1 = Smallworld;
field2 = 'x';  value2 = X;
field3 = 'y';  value3 =Y ;
field4 = 'nv';  value4 = 140;
field5 = 'ne';  value5 = edges;
s = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5);
plotGraphBasic(s,4,1);
%% 
%%Erwthma B 
Dlat=degrees(Lattice);
%[x,y,dista]=cumulativedist(Dlat,140);
figure(6);
hist(Dlat)
title('Histogram of degree distribution for Lattice');    


Dranderdos=degrees(RandomErdos);
figure(7);
hist(Dranderdos)
title('Histogram of degree distribution for Random Erdos');    


Dscalefree=degrees(ScaleFree);
figure(8);
hist(Dscalefree)
title('Histogram of degree distribution for ScaleFree');    


Drandgeo=degrees(RandomGeometrical);
figure(9);
hist(Drandgeo)
title('Histogram of degree distribution for Random Geometric');    


Dsworld=degrees(Smallworld);
figure(10);
hist(Dsworld)
title('Histogram of degree distribution for Smallworld');    


[supremumlat,cumdistlat,distlat]=cumulativedist(Dlat,140);
figure(1);
plot(cumdistlat)
title('Cummulative of degree distribution for Lattice');    

[supremumrerdos,cumdistrerdos,distrerdos]=cumulativedist(Dranderdos,140);
figure(2);
plot(cumdistrerdos)
title('Cummulative of degree distribution for Random Erdos');    

[supremumscfree,cumdistscfree,distscfree]=cumulativedist(Dscalefree,140);
figure(3);
plot(cumdistscfree)
title('Cummulative of degree distribution for Scale Free');    

[supremumrandgeo,cumdistrandgeo,distrandgeo]=cumulativedist(Drandgeo,140);
figure(4);
plot(cumdistrandgeo)
title('Cummulative of degree distribution for Random Geometric');    

[supremumsw,cumdistsw,distsw]=cumulativedist(Dsworld,140);
figure(5);
plot(cumdistsw)
title('Cummulative degree distribution for SmallWorld');    

Vlat= var(Dlat)
Mlat = mean(Dlat)

Vrerdos= var(Dranderdos)
Mrerdos=mean(Dranderdos)

Vscfree= var(Dscalefree) % histogram
Mscfree= mean(Dscalefree)

Vrgeo= var(Drandgeo)
Mrgeo= mean(Drandgeo)

Vsw= var(Dsworld)
Msw= mean(Dsworld)
%% erwthma G

for i=1:140
    for j=1:140
        Lattice(i,j)=10*rand;
    end
end
[NSLAT]=nodestrength(Lattice);
[kk,cumdistrandlat,ll]=cumulativedist(NSLAT,140);
figure(1);
plot (cumdistrandlat);
title('Cummulative Strength distribution for Lattice');

avgstrengthLAT=sum(NSLAT)/140
%h = cdfplot(NSLAT);
%[h,stats] = cdfplot(NSLAT);
for i=1:140
    for j=1:140
        RandomErdos(i,j)=10*rand;
    end
end
[NSRE]=nodestrength(RandomErdos);

[kk,cumdistre,ll]=cumulativedist(NSRE,140);
figure(2);
plot (cumdistre);
title('Cummulative Strength distribution for Random Erdos');

avgstrengthRerdos=sum(NSRE)/140
%h = cdfplot(NSRE);
%[h,stats] = cdfplot(NSRE);
for i=1:140
    for j=1:140
        ScaleFree(i,j)=10*rand;
    end
end

[NSSF]=nodestrength(ScaleFree);
[kk,cumdistsf,ll]=cumulativedist(NSSF,140);
figure(3);
plot (cumdistsf);
title('Cummulative Strength distribution for Scale Free');
avgstrengthSCALEfree=sum(NSSF)/140

%h = cdfplot(NSSF);
%[h,stats] = cdfplot(NSSF);
for i=1:140
    for j=1:140
        RandomGeometrical(i,j)=10*rand;
    end
end
[NSRG]=nodestrength(RandomGeometrical);
[kk,cumdistrg,ll]=cumulativedist(NSRG,140);
figure(4);
plot (cumdistrg);
title('Cummulative Strength distribution for Random Geometric');
avgstrengthRG=sum(NSRG)/140

%h = cdfplot(NSRG);
%[h,stats] = cdfplot(NSRG);
for i=1:140
    for j=1:140
        Smallworld(i,j)=10*rand;
    end
end
[NSSw]=nodestrength(Smallworld);
figure(5);
[kk,cumdistsw,ll]=cumulativedist(NSSw,140);
plot (cumdistsw);
title('Cummulative Strength distribution for Small world');
avgstrengthRG=sum(NSRG)/140

avgstrengthSmallworld=sum(NSSw)/140

%h = cdfplot(NSSw);
%[h,stats] = cdfplot(NSSw);
%% erwthma D
for i=1:5
    for j=1:2
        ComparisonM(i,j)=0;
    end
end

[Latapl]=ave_path_length(Lattice);
[SHlat] = graphallshortestpaths(Lattice);
ComparisonM(1,1)=Latapl;
ComparisonM(1,2)=var(SHlat(:));

[rerdapl]=ave_path_length(RandomErdos);
[SHrerdos]= graphallshortestpaths(RandomErdos);

ComparisonM(2,1)=rerdapl;
ComparisonM(2,2)=var(SHrerdos(:));

[sfreeapl]=ave_path_length(ScaleFree);
[SHscfree]=graphallshortestpaths(ScaleFree);
ComparisonM(3,1)=sfreeapl;
ComparisonM(3,2)=var(SHscfree(:));

[rgeoapl]=ave_path_length(RandomGeometrical);
[SHrgeo]=graphallshortestpaths(sparse(RandomGeometrical));
ComparisonM(4,1)=rgeoapl;
ComparisonM(4,2)=var(SHrgeo(:));

[swapl]=ave_path_length(Smallworld);
[SHsw]=graphallshortestpaths(Smallworld);
ComparisonM(5,1)=swapl;
ComparisonM(5,2)=var(SHsw(:));

%mesh timh variance
%% erwthma E
    [K,avgCLat,clLat]=clust_coeff(Lattice);
    [supremum,cumdist]=cumulativecentrality(clLat,140);
    figure(1);
    plot(cumdist)
    avgCLat
    
    [K,avgCRerdos,clRerdos]=clust_coeff(RandomErdos);
    [supremum,cumdist]=cumulativecentrality(clRerdos,140);   
    figure(2);
    plot(cumdist)
    avgCRerdos
    
    [K,avgCScalefree,clScalefree]=clust_coeff(ScaleFree);
    [supremum,cumdist]=cumulativecentrality(clScalefree,140);   
    figure(3);   
    plot(cumdist)
    avgCScalefree
    
    [K,avgCRandgeo,clRandgeo]=clust_coeff(sparse(RandomGeometrical));
    [supremum,cumdist]=cumulativecentrality(clRandgeo,140);   
    figure(4);
    plot(cumdist)
    avgCRandgeo
    
    [K,avgCSworld,clSworld]=clust_coeff(Smallworld);
    [supremum,cumdist]=cumulativecentrality(clSworld,140);   
    figure(5);
    plot(cumdist)
    avgCSworld
    

%% erwthma Z


figure;
regA=smallw(140,2,0);
degreg= degrees(full(regA));
[maximum,cumdistr] = cumulativecentrality(degreg,140);
stem([1:101],cumdistr,'filled','b');
figure;
closreg= closeness(full(regA));
[maximum,cumdistr] = cumulativecentrality(closreg,140);
stem([1:101],cumdistr,'filled','b');
figure;
betwreg=node_betweenness_faster(full(regA));
[maximum,cumdistr] = cumulativecentrality(betwreg,140);
stem([1:101],cumdistr,'filled','b');

figure;
eigenreg=eigencentrality(full(regA));
[maximum,cumdistr] = cumulativecentrality(eigenreg,140);
stem([1:101],cumdistr,'filled','b');

figure;


regA=erdrey(140,750);
degreg= degrees(full(regA));
[maximum,cumdistr] = cumulativecentrality(degreg,140);
stem([1:101],cumdistr,'filled','b');
figure;

closreg= closeness(full(regA));
[maximum,cumdistr] = cumulativecentrality(closreg,140);
stem([1:101],cumdistr,'filled','b');
figure;
betwreg=node_betweenness_faster(full(regA));
[maximum,cumdistr] = cumulativecentrality(betwreg,140);
stem([1:101],cumdistr,'filled','b');

figure;
eigenreg=eigencentrality(full(regA));
[maximum,cumdistr] = cumulativecentrality(eigenreg,140);
stem([1:101],cumdistr,'filled','b');


for i=1:140
     B(i,1)=1000*rand;
     B(i,2)=1000*rand;
end

figure;
regA=rgg(B,140,250);
degreg= degrees(full(regA));
[maximum,cumdistr] = cumulativecentrality(degreg,140);
stem([1:101],cumdistr,'filled','b');
figure;

closreg= closeness(full(regA));
[maximum,cumdistr] = cumulativecentrality(closreg,140);
stem([1:101],cumdistr,'filled','b');
figure;
betwreg=node_betweenness_faster(full(regA));
[maximum,cumdistr] = cumulativecentrality(betwreg,140);
stem([1:101],cumdistr,'filled','b');

figure;
eigenreg=eigencentrality(full(regA));
[maximum,cumdistr] = cumulativecentrality(eigenreg,140);
stem([1:101],cumdistr,'filled','b');


figure;
regA=pref(140,4);
degreg= degrees(full(regA));
[maximum,cumdistr] = cumulativecentrality(degreg,140);
stem([1:101],cumdistr,'filled','b');
figure;

closreg= closeness(full(regA));
[maximum,cumdistr] = cumulativecentrality(closreg,140);
stem([1:101],cumdistr,'filled','b');
figure;
betwreg=node_betweenness_faster(full(regA));
[maximum,cumdistr] = cumulativecentrality(betwreg,140);
stem([1:101],cumdistr,'filled','b');

figure;
eigenreg=eigencentrality(full(regA));
[maximum,cumdistr] = cumulativecentrality(eigenreg,140);
stem([1:101],cumdistr,'filled','b');



figure;
regA=smallw(140,2,0.3);

degreg= degrees(full(regA));
[maximum,cumdistr] = cumulativecentrality(degreg,140);
stem([1:101],cumdistr,'filled','b');
figure;

closreg= closeness(full(regA));
[maximum,cumdistr] = cumulativecentrality(closreg,140);
stem([1:101],cumdistr,'filled','b');
figure;
betwreg=node_betweenness_faster(full(regA));
[maximum,cumdistr] = cumulativecentrality(betwreg,140);
stem([1:101],cumdistr,'filled','b');

figure;
eigenreg=eigencentrality(full(regA));
[maximum,cumdistr] = cumulativecentrality(eigenreg,140);
stem([1:101],cumdistr,'filled','b');








%% erwthma H
clearvars A ConnPer100 ConnPer200;

c2=0;
c1=0;
Connectivity=[];
index=1;
for i=1:5
    for j=1:100
        Lattice1=smallw(100,i,0);
        if isconnected(Lattice1)
            c1=c1+1;
        end
        Lattice2=smallw(200,i,0);
        if isconnected(Lattice2)
            c2=c2+1;
        end
    end
    ConnPer100(index) = c1/100;
    ConnPer200(index) = c2/100;
     A(index)=i;
    index=index+1;
end
figure(1);
hold on;

plot(A,ConnPer100,'b');

plot(A,ConnPer200,'r');
title('connectivity percentage of Lattice to d');
xlabel('various values of d');
ylabel('Connectivity percentage');
% Prosthiki upomnimatos
legend('for 100 nodes', 'for 200 nodes','Location', 'SouthEast');
hold off;

clearvars A ConnPer100 ConnPer200;
c2=0;
c1=0;
index=1;
for i=0.1:0.1:0.9
   
    for j=1:100
        erR1=random_graph(100,i);
        if isconnected(erR1)
            c1=c1+1;
        end
        erR2=random_graph(200,i);
        if isconnected(erR2)
            c2=c2+1;
        end
          
    end
    ConnPer100(index)=c1/100;
    ConnPer200(index)=c2/100;
     A(index)=i;
    index=index+1;
end
figure(2);
hold on;

plot(A,ConnPer100,'b');

plot(A,ConnPer200,'r');
title('connectivity percentage for Erdos-Renyi to d');
xlabel('various values of d');
ylabel('Connectivity percentage');
% Prosthiki upomnimatos
legend('for 100 nodes', 'for 200 nodes','Location', 'SouthEast');
hold off;
clearvars A ConnPer100 ConnPer200;

c2=0;
c1=0;
 for i=100:100:800
    for j=1:100
        RandomErdos1=erdrey(100,i);
        if isconnected(RandomErdos1)
            c1=c1+1;
        end
        RandomErdos2=erdrey(200,i);
        if isconnected(RandomErdos2)
            c2=c2+1;
        end
    end
      ConnPer100(index)=c1/100;
    ConnPer200(index)=c2/100;
     A(index)=i;
    index=index+1;
end
figure(3);
hold on;

plot(A,ConnPer100,'b');

plot(A,ConnPer200,'r');
title('connectivity percentage for Gilbert to d');
xlabel('various values of d');
ylabel('Connectivity percentage');
% Prosthiki upomnimatos
legend('for 100 nodes', 'for 200 nodes','Location', 'SouthEast');
hold off;
clearvars A ConnPer100 ConnPer200;

c2=0;
c1=0;
index=1;
   for i=25:25:250
        for k=1:100
            for i=1:200
    x(i)=rand()*1000;
    y(i)=rand()*1000;
end

rgNodes200 = transpose([x ; y]);
clearvars x y;
         RandomGeometrical2=rgg(rgNodes200,200,i);
            if isconnected(RandomGeometrical2)
            c2=c2+1;
            end
        end
    ConnPer200(index)=c2/100;
    index=index+1;
    end
     %%c2=0;
     c1=0;
index=1;
for i=25:25:250
    for k=1:100
        for i=1:100
            x(i)=rand()*1000;
            y(i)=rand()*1000;
        end

         rgNodes100 = transpose([x ; y]);
         clearvars x y;
         RandomGeometrical1=rgg(rgNodes100,100,i);
        if isconnected(RandomGeometrical1)
            c1=c1+1;
        end
         
    end
         A(index)=i;

    ConnPer100(index)=c1/100;
    index=index+1;
end
figure(4);
hold on;

plot(A,ConnPer100,'b');

plot(A,ConnPer200,'r');
title('connectivity percentage to d');
xlabel('various values of d');
ylabel('Connectivity percentage');
% Prosthiki upomnimatos
legend('for 100 nodes', 'for 200 nodes','Location', 'SouthEast');
hold off;
clearvars A ConnPer100 ConnPer200;

c2=0;
c1=0;
index=1;    
for i=2:2:10
    for j=1:100
        ScaleFree1=pref(100,i);
        if isconnected( ScaleFree1)
            c1=c1+1;
        end
        ScaleFree2=pref(200,i);
        if isconnected(ScaleFree2)
            c2=c2+1;
        end
    end
      ConnPer100(index)=c1/100;
    ConnPer200(index)=c2/100;
     A(index)=i;
    index=index+1;
end
figure(4);
hold on;

plot(A,ConnPer100,'b');

plot(A,ConnPer200,'r');
title('connectivity percentage for Scale Free to d');
xlabel('various values of d');
ylabel('Connectivity percentage');
% Prosthiki upomnimatos
legend('for 100 nodes', 'for 200 nodes','Location', 'SouthEast');
hold off;
clearvars A ConnPer100 ConnPer200;

c2=0;
c1=0;
index=1;
for i=1:1:5
    for j=1:100
            Smallworld1=smallw(100,i,0.3);
        if isconnected(Smallworld1)
            c1=c1+1;
        end
            Smallworld2=smallw(200,i,0.3);
       if isconnected(Smallworld2)
            c2=c2+1;
        end
    end
      ConnPer100(index)=c1/100;
    ConnPer200(index)=c2/100;
     A(index)=i;
    index=index+1;
end
figure(5);
hold on;

plot(A,ConnPer100,'b');

plot(A,ConnPer200,'r');
title('connectivity percentage to d');
xlabel('various values of d');
ylabel('Connectivity percentage');
% Prosthiki upomnimatos
legend('for 100 nodes', 'for 200 nodes','Location', 'SouthEast');
hold off;
clearvars A ConnPer100 ConnPer200;

c2=0;
c1=0;
index=1;
    for j=0.1:0.1:0.7
        for l=1:100 
        Smallworld1=smallw(100,4,j);
        if isconnected(Smallworld1)
            c1=c1+1;
        end
            Smallworld2=smallw(200,4,j);
        if isconnected(Smallworld2)
            c2=c2+1;
        end
        end
       ConnPer100(index)=c1/100;
    ConnPer200(index)=c2/100;
     A(index)=i;
    index=index+1;
    end

figure(6);
hold on;

plot(A,ConnPer100,'b');

plot(A,ConnPer200,'r');
title('connectivity percentage to d');
xlabel('various values of d');
ylabel('Connectivity percentage');
% Prosthiki upomnimatos
legend('for 100 nodes', 'for 200 nodes','Location', 'SouthEast');
hold off;



%% erwthma I
l=1;
for i=0:0.1:1
    SW=smallw(140,2,i);
    apl(l)=ave_path_length(SW);
    [K,avgC(l),L]=clust_coeff(SW);
    l=l+1;
end
%%cl regula avpl tyxaiou
