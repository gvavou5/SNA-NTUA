function [ communities ] = MyGeneticAlgo(Adjacency,pc,pm,elitism,generations,chromosomes,FitMaxUnchanged)

	%Initialization
	nodes = size(Adjacency,1);
	B0 = cell(chromosomes,generations+1);
    current_generation = 1;
    FitnessPower = 1; % fitness power is 1 as said in the lab 
    
	neighbours = cell(nodes); %Neighbours of each node in a cell array
	for i=1:nodes
		neighbours{i} = find(Adjacency(i,:)==1);
    end
    
	for i=1:chromosomes
		for j=1:nodes
			B0{i,current_generation}(j) = neighbours{j}(randi(size(neighbours{j},2)));
		end
    end
	
	Adjacency1    = zeros(nodes,nodes); %New Adjacency Matrix
    fitnessMax    = -inf ;              % maximum of fitness function
	fitnessMaxCnt = 0;                  % if fitnessMax stays unchanged for 5 loops then break 
	
    
	while (current_generation <= generations+1 && fitnessMaxCnt < FitMaxUnchanged)  
		% fit_max_unchanged is 5 for this exercise
		for i = 1:chromosomes
			for j = 1:nodes
                % for every generation create the new Adjacency Matrix
                % according to B0
				Adjacency1(j,B0{i,current_generation}(j)) = 1;
				Adjacency1(B0{i,current_generation}(j),j) = 1;
            end
          
			connected_components = find_conn_comp(Adjacency1); %Find the connected components
            
			%calculation of fitness function and scores
			for j = 1:size(connected_components,2)
				Subgr{j}  = subgraph(Adjacency1,connected_components{j});
				CStemp(j) = fitness_function(Subgr{j},FitnessPower);
			end
			CS(i) = sum(CStemp); %CS for every chromosome stored here
			Adjacency1 = zeros(nodes,nodes);
		end
		
		%Elitism Selection
		[CS_sort,gg] = sort(CS,'descend');
        
		for i=1:elitism
			B0{i,current_generation+1} = B0{gg(i),current_generation};
        end
        
		%Roulette Method - Remaing Chromosomes
		CS_sum = sum(CS);
		for i = elitism+1:chromosomes 
			x = rand(1);
			k = 1;
			while (k < chromosomes && x < sum(CS(1:i))/CS_sum ) k=k+1; end
			B0{i,current_generation+1} = B0{k,current_generation};
		end
		
		%Crossover
		for i=1:2:chromosomes-1
			if (rand(1) <= pc)
				pos = randi(nodes-1);
				for k=pos+1:nodes
					aux = B0{i,current_generation+1}(k);
					B0{i,current_generation+1}(k) = B0{i+1,current_generation+1}(k);
					B0{i+1,current_generation+1}(k) = aux;				
				end
			end
		end
		  
		%Mutation
		for i=1:chromosomes
			for k=1:nodes
				if (rand(1) <= pm)
					B0{i,current_generation+1}(k) = neighbours{k}(randi(size(neighbours{k},2)));
				end
			end
		end
		
		if (CS_sort(1) == fitnessMax)
			fitnessMaxCnt = fitnessMaxCnt + 1;
		else
			fitnessMax = CS_sort(1);
			fitnessMaxCnt = 1;
		end
		
		current_generation = current_generation+1;
    end
    
	% Get the best result for the previous generation
	current_generation = current_generation - 1;
	Adjacency1 = zeros(nodes,nodes);
	for j = 1:nodes
		Adjacency1(j,B0{gg(1),current_generation}(j)) = 1;
		Adjacency1(B0{gg(1),current_generation}(j),j) = 1;
    end
    
	communities = find_conn_comp(Adjacency1); % cell array of the communities
    
end
