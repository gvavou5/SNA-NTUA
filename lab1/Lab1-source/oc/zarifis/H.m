%% RG - ER
figure;
connect100=zeros(8,1);
for j=1:8
    for i=1:100
        regA=erdrey(100,j*100);
        if isconnected(regA)
            connect100(j) = connect100(j)+1;
        end
    end
 
end

connect200=zeros(8,1);
for j=1:8
    for i=1:100
        regA=erdrey(200,j*100);
        if isconnected(regA)
            connect200(j) = connect200(j)+1;
        end
    end
end
connect100=connect100/100;
connect200=connect200/100;
connect=[connect100 connect200];
x=[100:100:800];
%b=bar(x,connect);
plot(x,connect);
title('Connected graph percent RG-ER');
legend('100 nodes','200 nodes');


%% RG (Gilbert)
figure;
connect100=zeros(9,1);
for j=1:9
    for i=1:100
        regA=random_graph(100,0.1*j);
        if isconnected(regA)
            connect100(j) = connect100(j)+1;
        end
    end
    connect100(j)=connect100(j)/100;
end

connect200=zeros(9,1);
for j=1:9
    for i=1:100
        regA=random_graph(200,0.1*j);
        if isconnected(regA)
            connect200(j) = connect200(j)+1;
        end
    end
    connect200(j)=connect200(j)/100;
end

connect=[connect100 connect200];
x=[0.1:0.1:0.9];
%b=bar(x,connect);
plot(x,connect);
title('Connect percent RG-G');
legend('100 nodes','200 nodes');


%% RGG
clearvars connect100 connect200;
figure;
clearvars x y;
for i=1:100
    x(i)=rand()*1000;
    y(i)=rand()*1000;
end

rgNodes100 = transpose([x ; y]);
clearvars x y;
for i=1:200
    x(i)=rand()*1000;
    y(i)=rand()*1000;
end

rgNodes200 = transpose([x ; y]);
clearvars x y;

connect100=zeros(10,1);
for j=1:10
    for i=1:100
        regA=rgg(rgNodes100,100,j*25);
        if isconnected(regA)
            connect100(j) = connect100(j)+1;
        end
    end
    connect100(j)=connect100(j)/100;
end

connect200=zeros(10,1);
for j=1:10
    for i=1:100
        regA=rgg(rgNodes200,200,j*25);
        if isconnected(regA)
            connect200(j) = connect200(j)+1;
        end
    end
    connect200(j)=connect200(j)/100;
end

connect=[connect100 connect200];
x=[25:25:250];
%b=bar(x,connect);
%legend('100 nodes','200 nodes');

plot(x,connect);
title('Connect percent RGG')
legend('100 nodes','200 nodes');

%% SF
figure;
connect100=zeros(5,1);
for j=1:5
    for i=1:100
        regA=pref(100,j*2);
        if isconnected(regA)
            connect100(j) = connect100(j)+1;
        end
    end
    connect100(j)=connect100(j)/100;
end

connect200=zeros(5,1);
for j=1:5
    for i=1:100
        regA=pref(200,j*2);
        if isconnected(regA)
            connect200(j) = connect200(j)+1;
        end
    end
    connect200(j)=connect200(j)/100;
end

connect=[connect100 connect200];
x=[2:2:10];
%b=bar(x,connect);
plot(x,connect);
title('Connect percent SF');
legend('100 nodes','200 nodes');


%% SW
figure;
connect100=zeros(5,1);
for j=1:5
    for i=1:100
        regA=smallw(100,j*2,0.1);
        if isconnected(regA)
            connect100(j) = connect100(j)+1;
        end
    end
    connect100(j)=connect100(j)/100;
end

connect200=zeros(5,1);
for j=1:5
    for i=1:100
        regA=smallw(200,j*2,0.1);
        if isconnected(regA)
            connect200(j) = connect200(j)+1;
        end
    end
    connect200(j)=connect200(j)/100;
end

connect=[connect100 ];
x=[2:2:10];
%b=bar(x,connect);
plot(x,connect);
title('Connect percent SW with g=0.1');
legend('100 nodes','200 nodes');



figure;
connect100=zeros(7,1);
for j=1:7
    for i=1:100
        regA=smallw(100,2,j*0.1);
        if isconnected(regA)
            connect100(j) = connect100(j)+1;
        end
    end
    connect100(j)=connect100(j)/100;
end

connect200=zeros(7,1);
for j=1:7
    for i=1:100
        regA=smallw(200,2,j*0.1);
        if isconnected(regA)
            connect200(j) = connect200(j)+1;
        end
    end
    connect200(j)=connect200(j)/100;
end

connect=[connect100 ];
x=[0.1:0.1:0.7];
plot(x,connect);
title('Connect percent for SW  d=2');
legend('100 nodes','200 nodes');


%% REG

figure;
connect100=zeros(5,1);
for j=1:5
    for i=1:100
        regA=smallw(100,j*2,0);
        if isconnected(regA)
            connect100(j) = connect100(j)+1;
        end
    end
    connect100(j)=connect100(j)/100;
end

connect200=zeros(5,1);
for j=1:5
    for i=1:100
        regA=smallw(200,j*2,0);
        if isconnected(regA)
            connect200(j) = connect200(j)+1;
        end
    end
    connect200(j)=connect200(j)/100;
end

connect=[connect100 connect200];
x=[2:2:10];
plot(x,connect)
title('Connect Percent REG');
legend('100 nodes','200 nodes');
