%%
%TODO***  make a function that creates V_adj, node_name,assigns G_comp.Nodes.Cluster{1}
%** always check if the cluster one should only contain 1 node as per the
%** V1 can only belong to cluster t1
%research paper
clear;
clc;
close all;
%%
V_adj =    [0     10    0     0;
            10     0    15     0;
            0     15     0     20;
            0      0     20     0];
        
        
 node_name =  {'V1','V2','V3','V4'};  
 
  
 G_init = digraph(V_adj, node_name);
 
 figure;
 plot(G_init,'EdgeLabel', G_init.Edges.Weight);
 title('Given Graph');
  
  
 [sz_r_G , sz_c_G] = size(V_adj);
 
 %%
 
 V_comp = V_adj; % completing the graph
 
 for i= 1:sz_r_G
     for j = 1:sz_c_G
        if ((j~=i) && (V_comp(i,j)==0)) 
          % if(i<j)
               V_comp(i,j) =  distances(G_init, i,j);%V_comp(i,j-1) + V_comp(j-1,j);
           %elseif (i>j)
          %     V_comp(i,j) =  ;V_comp(i-1,j) + V_comp(i,i-1);
          % end
            
        end
        
     end
 end
 
 G_comp = digraph(V_comp, node_name);
 %%
 figure;
 plot(G_comp, 'EdgeLabel', G_comp.Edges.Weight);
 title('Given Completed Graph');

%%
%node cluster relationship


Cluster_cell = cell(length(V_adj),1);

%G_comp.Nodes(:,'Cluster') = [];
G_comp.Nodes.Cluster = Cluster_cell;

%v1 => t1 ; v2 => t2&t5;  v3 => t3; v4 => t4&t5; 
%putting cluster data using only indices to do that
G_comp.Nodes.Cluster{1} = 1;
G_comp.Nodes.Cluster{2} = [5,2];
G_comp.Nodes.Cluster{3} = 3;
G_comp.Nodes.Cluster{4} = [5, 4];

%add replica nodes
%size( G_comp.Nodes.Cluster{1})
%str=sprintf('Similar Image - %d', top_5+1);

%asd = G_comp.Edges.EndNodes
%cx = ismember(asd,'V21')
%fasd = asd(cx(:,1),2)
%asd_w = G_comp.Edges.Weight
%fasd_w = asd_w(cx(:,1))
[sz_r_clst, sz_c_clust]= size(G_comp.Nodes.Cluster);
%%

%add new nodes which belong to more than one cluster
node_count = sz_r_clst;
for i = 1:sz_r_clst
    [r_rep, c_rep]= size(G_comp.Nodes.Cluster{i}); % checking if repeating nodes are there
    pre_cluster = G_comp.Nodes.Cluster{i};
    
    for j= 1:c_rep
        
        if(j == 1)
            str1 = G_comp.Nodes.Name{i};
            str  = sprintf('%s-%d', str1, G_comp.Nodes.Cluster{i}(j)); 
            G_comp.Nodes.Name{i} = str;
            G_comp.Nodes.Cluster{i} = G_comp.Nodes.Cluster{i}(j);
        else
            node_count = node_count + 1;
            str  = sprintf('%s-%d', str1, pre_cluster(j));
            G_comp = addnode(G_comp, str);            
            G_comp.Nodes.Cluster{node_count} = pre_cluster(j);
        end
    
    end
    %G_comp.Nodes.Name{2} = 'V20';
end

%%
given_nodes = G_comp.numnodes;
figure;
plot(G_comp);
title('Given Completed Graph with nodes added');


% add edges to the new nodes which are made in previous step

[sz_r_clst, sz_c_clust]= size(G_comp.Nodes.Cluster);

pre_edges = G_comp.Edges;%edges before new edges for new nodes are added
pre_endnode_tab = pre_edges.EndNodes;
pre_edgeweight_tab = pre_edges.Weight;

for i = 1:sz_r_clst
    
   if(ismember(G_comp.Nodes.Name{i}, pre_edges.EndNodes) == 0)
      % edge_match = ismember(pre_endnode_tab, G_comp.Nodes.Name{i}); 
       edge_match = ~cellfun(@isempty,(strfind(pre_endnode_tab,strtok(G_comp.Nodes.Name{i}, '-')))); %% strtok(- delimeter) on 'V1-5' -> 'V1' then find this 'V1' in the table output is logical array contaning atleast V1
       fin_node = pre_endnode_tab(edge_match(:,1), 2);
       start_node = cell(size(fin_node));   % ***step might take lot of time to run
       start_node(:) =  {G_comp.Nodes.Name{i}};
       edge_weight  = pre_edgeweight_tab(edge_match(:,1), 1);
       G_comp = addedge(G_comp, start_node, fin_node, edge_weight);
       G_comp = addedge(G_comp, fin_node, start_node, edge_weight); % directed edge back
   end
    
end

plot(G_comp, 'EdgeLabel', G_comp.Edges.Weight);
node_table = G_comp.Nodes;
%%
% adding zero cost edges between the new nodes

for i = 1:sz_r_clst
    
    if(node_table.Cluster{i}~=0)
        edge_match_copy = ~cellfun(@isempty,(strfind(node_table.Name,strtok(node_table.Name{i}, '-')))); % strtok(- delimeter) on 'V1-5' -> 'V1' then find this 'V1' in the table output is logical array contaning atleast V1
        if((sum(edge_match_copy(:)) > 1)) % check for repeated node
            copy_node = node_table.Name(edge_match_copy(:)) ;        
            node_table.Cluster(edge_match_copy) = {0};
            [r_copy c_copy] = size(copy_node);
            for j = 1:r_copy
                for k = 1:r_copy
                    if (k~=j)
                       G_comp = addedge(G_comp, copy_node(j), copy_node(k), 0); % adding zero cost edges 
                       % ***step will take a lot of time  can be altered               
                    end
                end
            end
            
        end
    end
    
end

plot(G_comp, 'EdgeLabel', G_comp.Edges.Weight);




%%
%adding the gadget nodes and edges


G_gadget = digraph([], []);
G_gadget = addnode(G_gadget, 'A,1-1');
G_gadget = addnode(G_gadget, 'C,1-1');

cluster_init = cell(2,1);
G_gadget.Nodes.Cluster = cluster_init;


G_gadget = addedge(G_gadget, {'C,1-1'}, {'A,1-1'}, 0);

G_gadget.Nodes.Cluster{1} = 1;
G_gadget.Nodes.Cluster{2} = 1;

[r_nodesize c_nodesize] = size(G_comp.Nodes);

exitvar =1;
prenode_tab = G_comp.Nodes;

for k = 1:r_nodesize
    % ensure that not run for cluster '1'
    if((k~=1) && (prenode_tab.Cluster{k}~=0)) % for k ==1 already created above and second condition along with  prenode_tab.Cluster(match_clus) = {0}; below ensures code is not repeated for already created nodes
        match_clus = ismember(cell2mat(prenode_tab.Cluster), [prenode_tab.Cluster{k}]); % logical array of same cluster members
        match_clusind = find(match_clus); % index where same cluster exist G_comp.Nodes
        num_nodecurclus = size(match_clusind,1); % number of nodes in current cluster
        for j =1:num_nodecurclus
            str_num = strtok(strtok(prenode_tab.Name{match_clusind(j)}, '-'), 'V'); %  'V' as delimeter strtok, check strtok usage above; output node number 
            str_clus = sprintf('%d', prenode_tab.Cluster{match_clusind(j)}); % output cluster number

            str_a  = sprintf('A,%s-%s', str_num, str_clus); 
            str_b  = sprintf('B,%s-%s', str_num, str_clus); 
            str_c  = sprintf('C,%s-%s', str_num, str_clus); 
            cur_clus = prenode_tab.Cluster{match_clusind(j)}; % cur_clus is the current cluster
            addabc_table = table({str_a str_b str_c}', {cur_clus cur_clus cur_clus}', 'VariableNames', {'Name' 'Cluster'});

            G_gadget = addnode(G_gadget, addabc_table);
            if(j==1)
                str_e  = sprintf('e-%s', str_clus); % e node added only once for a cluster
                adde_table = table({str_e}', {cur_clus}', 'VariableNames', {'Name' 'Cluster'});
                G_gadget = addnode(G_gadget, adde_table);
                str_aFirst = str_a;
                str_Oldc = str_c;           % the new c becomes old c for next cycle of abc...is used to connect cycles
            end
            add_edgetable = table({str_a str_b; str_b str_c; str_b str_a; str_c str_e; str_e str_b}, [0 0 0 0 0]', 'VariableNames',{'EndNodes','Weight'}); % adding a->b->c and b->a and c-> e and e-> b
            G_gadget = addedge(G_gadget, add_edgetable);
            if(j>1)
                con_prevGadgetTab =  table({str_Oldc str_a}, 0, 'VariableNames',{'EndNodes','Weight'});% connect current abc to previous abc via an edge from previous c to current c
                G_gadget = addedge(G_gadget, con_prevGadgetTab); % connecting abc_old-> abd_new
                str_Oldc = str_c;           % the new c becomes old c for next cycle
            end
            
            if(j==num_nodecurclus) % its the last node in current cluster then connect it like abc_last -> a_first
                con_lastToFirst = table({str_c str_aFirst}, 0, 'VariableNames',{'EndNodes','Weight'});
                G_gadget = addedge(G_gadget, con_lastToFirst); % connect it like abc_last -> a_first
            end
        
        end
        
        %declare the nodes and add edges abc->abc
        prenode_tab.Cluster(match_clus) = {0}; %this will ensure that matchclus is not run for same cluster again
    end
    
    
    
    
end
figure;

plot(G_gadget, 'EdgeLabel', G_gadget.Edges.Weight);

%**************************** check the edge between V2-5 and V4-5 these
%********should not exist after G-S(abc) transformation


for i = 1:size(G_comp.Edges.EndNodes,1)
    init_NodeClus = strsplit(cell2mat(strtok(G_comp.Edges.EndNodes(i,1), 'V')), '-');
    end_NodeClus  = strsplit(cell2mat(strtok(G_comp.Edges.EndNodes(i,2), 'V')), '-');
    
    if(~isequal(init_NodeClus{2}, end_NodeClus{2}))
        str_a_inter = sprintf('A,%s-%s', init_NodeClus{1}, init_NodeClus{2}); % a node for intercluster edge  - init_NodeClus{1} is node name and init_Nodeclus{2} is cluster
        str_c_inter = sprintf('C,%s-%s', end_NodeClus{1}, end_NodeClus{2}); % c node for intercluster edge
        interEdge_table =  table({str_a_inter str_c_inter}, G_comp.Edges.Weight(i), 'VariableNames',{'EndNodes','Weight'});
        G_gadget = addedge(G_gadget, interEdge_table); % 
    end
end



figure;

LWidths = 5*(G_gadget.Edges.Weight+1)/max(G_gadget.Edges.Weight);
p_gadget = plot(G_gadget, 'Layout','circle','EdgeLabel', G_gadget.Edges.Weight, 'LineWidth',LWidths);
p_gadget.NodeColor = 'r';


%%
%*** directed to undirected graph

 % G_gadget2 will be converted to undirected its a copy of G_gadget but
 % with undirected edges
 temp_G = G_gadget;
 
 G_gadget2 = graph([], []);


 [r_Gnodes c_Gnodes] = size(temp_G.Nodes);
 [r_Gedge c_Gedge] = size(temp_G.Edges);
 
 undir_table = table(cell(r_Gnodes*3,1), 'VariableNames',{'Name'});
 undir_edgetab = table(cell(r_Gnodes*2,2), zeros(r_Gnodes*2,1), 'VariableNames', {'EndNodes','Weight'}); % has only internal edges  1;A,1-1 -> 2;A,1-1 and 2;A,1-1 -> 3;A,1-1 of this type
 % for every node we have 2 edges(to make it undirected) and all already
 % existing edges in G_gadget we have a edge in G_gadget2
 tab_counter = 0; % its actually tab node counter
 tab_edgecounter = 1 ; % its tab edge
 
% adding new nodes for undirected graph A,1-1 becomes 1;A,1-1   2;A,1-1 3;A,1-1
% also adding edges 1;A,1-1 -> 2;A,1-1 and 2;A,1-1 -> 3;A,1-1

for i = 1:r_Gnodes
    
    
    for j=1:3
        tab_counter = tab_counter+1;
        str_cur = sprintf('%d;%s', j, G_gadget.Nodes.Name{i}); % A,1-1 becomes 1;A,1-1   2;A,1-1 3;A,1-1
        undir_table(tab_counter, 'Name') =  table({str_cur}, 'VariableNames',{'Name'});% making the undirected table and will add to the graph together as this saves time
        
    end
    
    
    undir_edgetab(tab_edgecounter:(tab_edgecounter+1), :)  = table([undir_table{(tab_counter-2),'Name'} undir_table{(tab_counter-1),'Name'}...
                                         ;undir_table{(tab_counter-1),'Name'} undir_table{tab_counter,'Name'}] , [0 0]',...
                                         'VariableNames',{'EndNodes','Weight'});   % edges 1;A,1-1 -> 2;A,1-1 and 2;A,1-1 -> 3;A,1-1
    
    tab_edgecounter = tab_edgecounter+2; % as we are adding 2 rows in the table at once
    
    
end

G_gadget2 = addnode(G_gadget2, undir_table); % we had made the tables above just putting those in the graph
G_gadget2 = addedge(G_gadget2, undir_edgetab); %we had made the tables above just putting those in the graph
%%

% ** we will add edges for an edge A,1-1 to C,1-1 the new edge 

undir_edgetab_ex =  table(cell(r_Gedge,2), zeros(r_Gedge,1), 'VariableNames', {'EndNodes','Weight'}); % external edge tab


for i = 1:r_Gedge
    
    
    
    str_1node = sprintf('%d;%s', 1, temp_G.Edges.EndNodes{i,1}); 
    
    str_3node = sprintf('%d;%s', 3, temp_G.Edges.EndNodes{i,2}); 
    
    
    
    undir_edgetab_ex(i, :) =   table({str_1node str_3node}, temp_G.Edges.Weight(i), 'VariableNames', {'EndNodes','Weight'});
    
    
    
    
end

figure;
G_gadget2 = addedge(G_gadget2, undir_edgetab_ex); %we had made the tables above just putting those in the graph

LWidths2 = 5*(G_gadget2.Edges.Weight+1)/max(G_gadget2.Edges.Weight);
p_gadget2 = plot(G_gadget2, 'Layout','circle','EdgeLabel', G_gadget2.Edges.Weight, 'LineWidth',LWidths2);
p_gadget2.NodeColor = 'r';


%% 
%completed graph just so that we can feed it into tsp solver

[s_gadget2 t_gadget2] = findedge(G_gadget2); % s - start -> t- target this is a list of edges in terms of index


tot_edgeweight = sum(G_gadget2.Edges.Weight);
gadget2_compMat = (tot_edgeweight*50+100)*ones(G_gadget2.numnodes, G_gadget2.numnodes) - diag((tot_edgeweight*50+99)*ones(1,G_gadget2.numnodes));



for i = 1:G_gadget2.numedges
    
     gadget2_compMat(s_gadget2(i),t_gadget2(i))  =  G_gadget2.Edges.Weight(i);
     gadget2_compMat(t_gadget2(i), s_gadget2(i))  =  G_gadget2.Edges.Weight(i);
        
     
    
end
    



[Out_sol] = TSP_tour_EXPLICIT(gadget2_compMat,'/home/ashishkb/softwares/concorde/concorde_build/TSP/concorde');
Out_solName = G_gadget2.Nodes.Name(Out_sol);

%%
% following we are doing to highlight the solution
s_sol = zeros(G_gadget2.numnodes,1);
t_sol = zeros(G_gadget2.numnodes,1);

for i = 2:G_gadget2.numnodes
    s_sol(i-1) = Out_sol(i-1);
    t_sol(i-1) = Out_sol(i);
end

s_sol(G_gadget2.numnodes) = t_sol(end-1);
t_sol(G_gadget2.numnodes) = s_sol(1);
%%
%plot the tour in the graph
figure;
LWidths2 = 5*(G_gadget2.Edges.Weight+1)/max(G_gadget2.Edges.Weight);
p_gadget2 = plot(G_gadget2, 'Layout','circle','EdgeLabel', G_gadget2.Edges.Weight, 'LineWidth',LWidths2);
p_gadget2.NodeColor = 'r';
highlight(p_gadget2,s_sol,t_sol, 'NodeColor','g','EdgeColor','r'); % if given edges solution is wrong this will not highlight the solution

%1;B,2-5 means ->'1' in the start is part of 3 nodes made for undirected  'B'-> means gadget applied for abc as in the paper next '2' means he node number in given graph next '5' means the cluster the node belongs to 
% {'1;A,1-1';'2;A,1-1';'3;A,1-1';'1;C,1-1';'2;C,1-1';'3;C,1-1';'1;A,2-2';'2;A,2-2';'3;A,2-2';'1;B,2-2';'2;B,2-2';'3;B,2-2';'1;e-2';'2;e-2';'3;e-2';'1;C,2-2';'2;C,2-2';'3;C,2-2';'1;A,2-5';'2;A,2-5';'3;A,2-5';'1;B,2-5';'2;B,2-5';'3;B,2-5';'1;e-5';'2;e-5';'3;e-5';'1;C,4-5';'2;C,4-5';'3;C,4-5';'1;B,4-5';'2;B,4-5';'3;B,4-5';'1;A,4-5';'2;A,4-5';'3;A,4-5';'1;C,2-5';'2;C,2-5';'3;C,2-5';'1;A,3-3';'2;A,3-3';'3;A,3-3';'1;B,3-3';'2;B,3-3';'3;B,3-3';'1;e-3';'2;e-3';'3;e-3';'1;C,3-3';'2;C,3-3';'3;C,3-3';'1;A,4-4';'2;A,4-4';'3;A,4-4';'1;B,4-4';'2;B,4-4';'3;B,4-4';'1;e-4';'2;e-4';'3;e-4';'1;C,4-4';'2;C,4-4';'3;C,4-4'}
% 
% 
%             V1-1                          ->                                         V2-2                                                                                   ->         V2-5                                                       ->              V4-5                                                                         ->      V2-5                             ->              V3-3                                                                                         ->                       V4-4                                                                                                    ->   V1-1  
%            V1-1  -> V2-2  ->  V2-5    ->   V4-5 -> V2-5   ->    V3-3   ->    V4-4  ->   V1-1  



%%
% get the final solution:

    

 fin_sol = cell(ceil(length(Out_solName)/9),1);  % final solution in a cell array divide by 9 because we are removing abc=3(nodes) and 123 = 3(used for undirected) therefore 3*3 = 9
 fin_rm_redunt = cell(length(given_nodes),1); % remove redundant and repeated nodes from fin_sol 
 sol_count = 1;
 prev_node = '';
for i = 1:length(Out_solName)    
    
    split_sol = strsplit(Out_solName{i},','); % splitting 1;A,1-1 to give cell(2,1) with cell{1} = 1;A and cell{2} = 1-1
    
    if(~ ismember('e',Out_solName{i}))
        addnode = sprintf('V%s',split_sol{2}); % 1-1 becomes V1-1
    end
    
    if(~isequal(addnode, prev_node))  % checking if new node being added wasn't added in the previous step itself*********  % this was previous condition (~ismember(addnode, fin_sol))&& (~ ismember('e',Out_solName{i})
        
        fin_sol{sol_count} = addnode;
        
        sol_count = sol_count+1;
    end
    
    prev_node = addnode;
    
end

fin_rm_redunt(:) = {''};
finsol_count = 1;
solnode = cell(length(fin_sol),1);
solclus = cell(length(fin_sol),1);

for i = 1:length(fin_sol)
    
    if(i==1)
        fin_rm_redunt{finsol_count}  = strtok(fin_sol{i}, '-');
                 cur_node = strtok(fin_rm_redunt{finsol_count},'V');
        cur_clus = strsplit(fin_sol{i}, '-');
        solclus{finsol_count} = cur_clus{2};
        solnode{finsol_count} = cur_node;
        finsol_count = finsol_count+1;
    else
         cur_node  = strtok(strtok(fin_sol{i}, '-'),'V');
         cur_clus = strsplit(fin_sol{i}, '-');
         solclus{i} = cur_clus{2};       
         solnode{i} = cur_node;
         if(~((cur_clus{2} == solclus{i-1}) || (cur_node == solnode{i-1})))
             fin_rm_redunt{finsol_count}  = strtok(fin_sol{i}, '-');             
             finsol_count = finsol_count+1;
         end
         
    end
end


%%