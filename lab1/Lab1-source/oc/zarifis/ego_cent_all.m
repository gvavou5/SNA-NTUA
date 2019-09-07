function ego_rg = ego_cent_all(rgA,n)

for i=1:n
    T=[i,kneighbors(rgA,i,1)];
    temp= subgraph(rgA,T);
    len=size(full(T));
    ego_rg(i)= ego_net(full(temp),len(2));
    clearvars len temp T;
end

end

