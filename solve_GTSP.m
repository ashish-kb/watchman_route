 % V_adj =    [0     10    0     0;
%             10     0    15     0;
%             0     15     0     20;
%             0      0     20     0];
% V_adj nxn symmetric matrix where n is number of nodes and values of matrix indicate
% an edge of that weight, value '0' means no edge

% Cluster is cell of nx1 where n is number of nodes each entry of cell
% indicates the cluster relationship also can not take column matrix
% V_Cluster{1} = [3, 1, 4]; %node 1 belongs to cluster 3,1&4 
% V_Cluster{2} = [5, 2, 1];
% V_Cluster{3} = [3, 5, 6];
% V_Cluster{4} = [5, 3, 4, 6];
% 


function [fin_sol, fin_rm_redunt, Out_solName, Out_sol, G_init, G_gadget, G_gadget2, time_concorde_struct] = solve_GTSP(V_adj, V_Cluster) % the fomat is given above

    given_nodes = length(V_adj);
    node_name  = cell(1, given_nodes);
    
    for i = 1:given_nodes
               
        node_name{i} = sprintf('V%d',i);
    
    end
    
    

    G_init = digraph(V_adj, node_name);
 
    figure;
    P_init = plot(G_init,'EdgeLabel', G_init.Edges.Weight);
    title('Given Graph');


    [sz_r_G , sz_c_G] = size(V_adj);
    
    
%%
    % completing the graph
     V_comp = V_adj; 
 
%      for i= 1:sz_r_G
%          for j = 1:sz_c_G
%             if ((j~=i) && (V_comp(i,j)==0)) 
%                V_comp(i,j) =  distances(G_init, i,j);
%             end
% 
%          end
%      end

     V_comp =  distances(G_init, 1:sz_r_G,1:sz_c_G);


     G_comp = digraph(V_comp, node_name);
    

      figure;
      plot(G_comp, 'EdgeLabel', G_comp.Edges.Weight);
      title('Completed Graph');

%some plots in the plot_script
     G_comp.Nodes.Cluster = V_Cluster;
     
     ind_empty_clus = cellfun(@isempty, G_comp.Nodes.Cluster);  % this contains indices of graph where there is no cluster relationship means orphan nodes
     G_comp_tm = rmnode(G_comp, find(ind_empty_clus));
     
%%

    %add new nodes which belong to more than one cluster
    G_comp = G_comp_tm;
    IN_tab_Nodes = IN_transform_Nodes(G_comp);

    G_comp_temp = digraph([], []);
%     G_comp_temp1 = digraph([], []);
   

    G_comp_temp = addnode(G_comp_temp, IN_tab_Nodes); 
    


    figure;
    plot(G_comp);
    title('Given Completed Graph - without - duplicate nodes added and orphan nodes removed');
    
    given_nodes = G_comp.numnodes; 
    
    % add edges to the new nodes which are made in previous step

    [IN_Interedge_s, IN_Interedge_e, IN_Interedge_w] = IN_transform_Interedge(G_comp);
    
    
    IN_Interedge_s = IN_Interedge_s(~cellfun('isempty',IN_Interedge_s(:,1))); % remove the empty entries from the CELL
    IN_Interedge_e = IN_Interedge_e(~cellfun('isempty',IN_Interedge_e(:,1)));
    IN_Interedge_w = IN_Interedge_w(IN_Interedge_w>0);
    
    
%     IN_tab_Interedge = IN_tab_Interedge(~cellfun('isempty',IN_tab_Interedge.EndNodes(:,1)),:); % remove the empty entries from the table
    
%     G_comp_temp = addedge(G_comp_temp, IN_tab_Interedge);
    G_comp_temp = addedge(G_comp_temp, IN_Interedge_s, IN_Interedge_e, IN_Interedge_w);     % fin_sol = fin_sol(~cellfun('isempty',fin_sol)); % removing empty cells
    
    % % adding zero cost edges between the new nodes

    IN_tab_Zeroedge = IN_transform_Zeroedge(G_comp);

    G_comp_temp = addedge(G_comp_temp, IN_tab_Zeroedge);
    
    figure;
    title('Interedges, duplicate nodes and zero cost edges');
    plot(G_comp_temp, 'EdgeLabel', G_comp_temp.Edges.Weight);
    
    G_comp = G_comp_temp;

    G_gadget = digraph([], []);
    G_gadget1 = digraph([], []);
    
    [r_nodesize, ~] = size(G_comp.Nodes);

    exitvar =1;
    prenode_tab = G_comp.Nodes;
    
    add_edgetable =  table({}, [], 'VariableNames',{'EndNodes','Weight'});
    con_prevGadgetTab =  table({}, [], 'VariableNames',{'EndNodes','Weight'});
    con_lastToFirst =  table({}, [], 'VariableNames',{'EndNodes','Weight'});

    for k = 1:r_nodesize
        
        if( (prenode_tab.Cluster{k}~=0)) %
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
%                 addabc_table = table({str_a str_b str_c}', {cur_clus cur_clus cur_clus}', 'VariableNames', {'Name' 'Cluster'});
%                 addabc_table1 = [addabc_table1;table({str_a str_b str_c}', {cur_clus cur_clus cur_clus}', 'VariableNames', {'Name' 'Cluster'})];
                
%                 G_gadget = addnode(G_gadget, addabc_table);
                if(j==1)
                    str_e  = sprintf('e-%s', str_clus); % e node added only once for a cluster
%                     adde_table = table({str_e}', {cur_clus}', 'VariableNames', {'Name' 'Cluster'});
%                     adde_table1 = [adde_table1;table({str_e}', {cur_clus}', 'VariableNames', {'Name' 'Cluster'})];
%                     G_gadget = addnode(G_gadget, adde_table);
                    str_aFirst = str_a;
                    str_Oldc = str_c;           % the new c becomes old c for next cycle of abc...is used to connect cycles
                end
%                 add_edgetable = table({str_a str_b; str_b str_c; str_b str_a; str_c str_e; str_e str_b}, [0 0 0 0 0]', 'VariableNames',{'EndNodes','Weight'}); % adding a->b->c and b->a and c-> e and e-> b
                add_edgetable = [add_edgetable;table({str_a str_b; str_b str_c; str_b str_a; str_c str_e; str_e str_b}, [0 0 0 0 0]', 'VariableNames',{'EndNodes','Weight'})];
%                 G_gadget = addedge(G_gadget, add_edgetable);
                if(j>1)
%                     con_prevGadgetTab =  table({str_Oldc str_a}, 0, 'VariableNames',{'EndNodes','Weight'});% connect current abc to previous abc via an edge from previous c to current c
                    con_prevGadgetTab = [con_prevGadgetTab; table({str_Oldc str_a}, 0, 'VariableNames',{'EndNodes','Weight'})];
%                     G_gadget = addedge(G_gadget, con_prevGadgetTab); % connecting abc_old-> abd_new
                    str_Oldc = str_c;           % the new c becomes old c for next cycle
                end

                if(j==num_nodecurclus) % its the last node in current cluster then connect it like abc_last -> a_first
%                     con_lastToFirst = table({str_c str_aFirst}, 0, 'VariableNames',{'EndNodes','Weight'});
                    con_lastToFirst = [con_lastToFirst;table({str_c str_aFirst}, 0, 'VariableNames',{'EndNodes','Weight'})];
%                     G_gadget = addedge(G_gadget, con_lastToFirst); % connect it like abc_last -> a_first
                end

            end

            %declare the nodes and add edges abc->abc
            prenode_tab.Cluster(match_clus) = {0}; %this will ensure that matchclus is not run for same cluster again
        end




    end
%     addabc_table1
%     adde_table1
    G_gadget = addedge(G_gadget, [add_edgetable; con_prevGadgetTab; con_lastToFirst]);
%     add_edgetable1
%     con_prevGadgetTab1
%     con_lastToFirst1
      figure;
% % 
    G_gadget_temp = G_gadget;
     P_gadge = plot(G_gadget_temp, 'EdgeLabel', G_gadget_temp.Edges.Weight);

    %**************************** check the edge between V2-5 and V4-5 these
    %********should not exist after G-S(abc) transformation
    % adding interedge nodes
%     interEdge_table  = table({}, [], 'VariableNames',{'EndNodes','Weight'});
    init_Clus_vec = cell2mat(cellfun(@(x) str2double(x((regexp(x,'-','start')+1):end)) , G_comp.Edges.EndNodes(:,1),'uni',0));
    end_Clus_vec = cell2mat(cellfun(@(x) str2double(x((regexp(x,'-','start')+1):end)) , G_comp.Edges.EndNodes(:,2),'uni',0));
    
    init_node_vec = cell2mat(cellfun(@(x) str2double(x(2:(regexp(x,'-','start')-1))) , G_comp.Edges.EndNodes(:,1),'uni',0));
    end_node_vec  = cell2mat(cellfun(@(x) str2double(x(2:(regexp(x,'-','start')-1))) , G_comp.Edges.EndNodes(:,2),'uni',0));
    
    str_a_inter = arrayfun(@(x,y) sprintf('A,%d-%d', x, y), init_node_vec, init_Clus_vec,'uni', 0);
    str_c_inter = arrayfun(@(x,y) sprintf('C,%d-%d', x, y), end_node_vec, end_Clus_vec,'uni', 0);
    G_gadget = addedge(G_gadget, str_a_inter,str_c_inter,G_comp.Edges.Weight(:));
 %   for i = 1:size(G_comp.Edges.EndNodes,1)
%         init_NodeClus =
%         strsplit(cell2mat(strtok(G_comp.Edges.EndNodes(i,1), 'V')), '-'); 
%         end_NodeClus  = strsplit(cell2mat(strtok(G_comp.Edges.EndNodes(i,2), 'V')), '-');
%          init_NodeClus = strsplit(G_comp.Edges.EndNodes{i,1}(2:end), '-');
%          end_NodeClus = strsplit(G_comp.Edges.EndNodes{i,2}(2:end), '-');
%         init_NodeClus = init_NodeClus_vec{i};
%         end_NodeClus  = end_NodeClus_vec{i};
%        if(~isequal(init_NodeClus{2}, end_NodeClus{2}))     % condition
%        removed as such edges were already removed
%            str_a_inter = sprintf('A,%s-%s', init_NodeClus{1}, init_NodeClus{2}); % a node for intercluster edge  - init_NodeClus{1} is node name and init_Nodeclus{2} is cluster
%            str_c_inter = sprintf('C,%s-%s', end_NodeClus{1}, end_NodeClus{2}); % c node for intercluster edge
          

%            interEdge_table =  [interEdge_table;table({str_a_inter str_c_inter}, G_comp.Edges.Weight(i), 'VariableNames',{'EndNodes','Weight'})];
            
%        end
%    end
%str_a_inter = arrayfun(@(x,y) sprintf('A,%d-%d', x, y), )
     % 

     figure;
% 
     LWidths = 5*(G_gadget.Edges.Weight+1)/max(G_gadget.Edges.Weight);
     P_gadget = plot(G_gadget, 'XData',  P_gadge.XData, 'YData',  P_gadge.YData);
     %P_gadget = plot(G_gadget, 'Layout','circle','EdgeLabel', G_gadget.Edges.Weight, 'LineWidth',LWidths); %this is when we don't have XData and YData
     P_gadget.NodeColor = 'r';


    %%
    %*** directed to undirected graph

     % G_gadget2 will be converted to undirected its a copy of G_gadget but
     % with undirected edges
     temp_G = G_gadget;

     G_gadget2 = graph([], []);
     
     str_1node_intra = cellfun(@(x,y) sprintf('1;%s', x), G_gadget.Nodes.Name,'uni', 0);
     str_2node_intra = cellfun(@(x,y) sprintf('2;%s', x), G_gadget.Nodes.Name,'uni', 0);
     str_3node_intra = cellfun(@(x,y) sprintf('3;%s', x), G_gadget.Nodes.Name,'uni', 0);
     
     G_gadget2 = addedge(G_gadget2, [str_1node_intra;str_2node_intra], [str_2node_intra;str_3node_intra], zeros(2*G_gadget.numnodes,1));

%      [r_Gnodes c_Gnodes] = size(temp_G.Nodes);
%      [r_Gedge c_Gedge] = size(temp_G.Edges);
% 
%      undir_table = table(cell(r_Gnodes*3,1), 'VariableNames',{'Name'});
%      undir_edgetab = table(cell(r_Gnodes*2,2), zeros(r_Gnodes*2,1), 'VariableNames', {'EndNodes','Weight'}); % has only internal edges  1;A,1-1 -> 2;A,1-1 and 2;A,1-1 -> 3;A,1-1 of this type
%      % for every node we have 2 edges(to make it undirected) and all already
%      % existing edges in G_gadget we have a edge in G_gadget2
%      tab_counter = 0; % its actually tab node counter
%      tab_edgecounter = 1 ; % its tab edge
% 
%     % adding new nodes for undirected graph A,1-1 becomes 1;A,1-1   2;A,1-1 3;A,1-1
%     % also adding edges 1;A,1-1 -> 2;A,1-1 and 2;A,1-1 -> 3;A,1-1
% 
%     for i = 1:r_Gnodes
% 
% 
%         for j=1:3
%             tab_counter = tab_counter+1;
%             str_cur = sprintf('%d;%s', j, G_gadget.Nodes.Name{i}); % A,1-1 becomes 1;A,1-1   2;A,1-1 3;A,1-1
%             undir_table(tab_counter, 'Name') =  table({str_cur}, 'VariableNames',{'Name'});% making the undirected table and will add to the graph together as this saves time
% 
%         end
% 
% 
%         undir_edgetab(tab_edgecounter:(tab_edgecounter+1), :)  = table([undir_table{(tab_counter-2),'Name'} undir_table{(tab_counter-1),'Name'}...
%                                              ;undir_table{(tab_counter-1),'Name'} undir_table{tab_counter,'Name'}] , [0 0]',...
%                                              'VariableNames',{'EndNodes','Weight'});   % edges 1;A,1-1 -> 2;A,1-1 and 2;A,1-1 -> 3;A,1-1
% 
%         tab_edgecounter = tab_edgecounter+2; % as we are adding 2 rows in the table at once
% 
% 
%     end
    
 
%     G_gadget22 = graph([], []);
    
   % G_gadget22 = addnode(G_gadget22, [str_1node_intra;str_2node_intra;str_3node_intra]);
    
    %G_gadget22 = addedge(G_gadget22, , str_3node_intra);
%     G_gadget2  = G_gadget22;
    
 % need not  to add nodes adding edges does that automatically  %G_gadget2 = addnode(G_gadget2, undir_table); % we had made the tables above just putting those in the graph
%      G_gadget2 = addedge(G_gadget2, undir_edgetab); %we had made the tables above just putting those in the graph
    %%

    % ** we will add edges for an edge A,1-1 to C,1-1 the new edge becomes
    % 1;A,1-1 to 3;C,1-1

  %  undir_edgetab_ex =  table(cell(r_Gedge,2), zeros(r_Gedge,1), 'VariableNames', {'EndNodes','Weight'}); % external edge tab
    %str_a_inter = arrayfun(@(x,y) sprintf('A,%d-%d', x, y), init_node_vec, init_Clus_vec,'uni', 0);
    
    str_1node = cellfun(@(x,y) sprintf('1;%s', x), temp_G.Edges.EndNodes(:,1),'uni', 0);
    str_3node = cellfun(@(x,y) sprintf('3;%s', x), temp_G.Edges.EndNodes(:,2),'uni', 0);
    
    
    G_gadget2 = addedge(G_gadget2, str_1node,str_3node,temp_G.Edges.Weight(:)); %
%     for i = 1:r_Gedge
% 
% 
% 
%         str_1node = sprintf('%d;%s', 1, temp_G.Edges.EndNodes{i,1}); 
% 
%         str_3node = sprintf('%d;%s', 3, temp_G.Edges.EndNodes{i,2}); 
% 
% 
% 
%         undir_edgetab_ex(i, :) =   table({str_1node str_3node}, temp_G.Edges.Weight(i), 'VariableNames', {'EndNodes','Weight'});
% 
% 
% 
% 
%     end
    
    
    

 %   figure;
    

%     LWidths2 = 5*(G_gadget2.Edges.Weight+1)/max(G_gadget2.Edges.Weight);
%     p_gadget2 = plot(G_gadget2, 'Layout','circle','EdgeLabel', G_gadget2.Edges.Weight, 'LineWidth',LWidths2);
%     p_gadget2.NodeColor = 'r';
% 

    %% 
    %completed graph just so that we can feed it into tsp solver
% 
%     [s_gadget2 t_gadget2] = findedge(G_gadget2); % s - start -> t- target this is a list of edges in terms of index
% 
% 
%     tot_edgeweight = sum(G_gadget2.Edges.Weight);
%     gadget2_compMat = (tot_edgeweight+100)*ones(G_gadget2.numnodes, G_gadget2.numnodes) - diag((tot_edgeweight+99)*ones(1,G_gadget2.numnodes));
% 
% %     gadget2_compMat1 = gadget2_compMat;
% 
% %     for i = 1:G_gadget2.numedges
% % 
% %          gadget2_compMat(s_gadget2(i),t_gadget2(i))  =  G_gadget2.Edges.Weight(i);
% %          gadget2_compMat(t_gadget2(i), s_gadget2(i))  =  G_gadget2.Edges.Weight(i);
% % 
% % 
% % 
% %     end
%     ind_gadget2 = sub2ind(size(gadget2_compMat), s_gadget2(:), t_gadget2(:));
%     ind_gadget2 = [ind_gadget2;sub2ind(size(gadget2_compMat), t_gadget2(:), s_gadget2(:))];
%     
%     weight_feed_mat = G_gadget2.Edges.Weight(:); % will feed this below
%     gadget2_compMat(ind_gadget2(:)) = [weight_feed_mat; weight_feed_mat];
    


  %   [Out_sol, time_concorde_struct] = TSP_tour_EXPLICIT(gadget2_compMat,'/home/ashishkb/softwares/concorde/concorde/TSP/concorde');
            
    [Out_sol, time_concorde_struct] = TSP_tour_Dat(G_gadget2,'/home/ashishkb/softwares/concorde/concorde/TSP/concorde');
    
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
%    plot the tour in the graph
%     figure;
%     LWidths2 = 5*(G_gadget2.Edges.Weight+1)/max(G_gadget2.Edges.Weight);
%     p_gadget2 = plot(G_gadget2, 'Layout','Circle','EdgeLabel', G_gadget2.Edges.Weight, 'LineWidth',LWidths2);
%     p_gadget2.NodeColor = 'r';
%     highlight(p_gadget2,s_sol,t_sol, 'NodeColor','g','EdgeColor','r'); % if given edges solution is wrong this will not highlight the solution
%     %%
%     % plot spanning tree
%      figure;
% %    % LWidths2 = 5*(G_gadget2.Edges.Weight+1)/max(G_gadget2.Edges.Weight);
%      p_gadget2_T = plot(G_gadget2_T, 'NodeLabel',G_gadget2_T.Nodes.Name,'Layout','force','EdgeLabel', G_gadget2_T.Edges.Weight); %, 'LineWidth',LWidths2);
%      p_gadget2_T.NodeColor = 'r';
  %  highlight(p_gadget2,s_sol,t_sol, 'NodeColor','g','EdgeColor','r'); % if given edges solution is wrong this will not highlight the solution


    %%
    % get the final solution:



     fin_sol = cell(ceil(length(Out_solName)/9),1);  % final solution in a cell array divide by 9 because we are removing abc=3(nodes) and 123 = 3(used for undirected) therefore 3*3 = 9
     fin_rm_redunt = cell(length(given_nodes),1); % remove redundant and repeated nodes from fin_sol 
     sol_count = 1;
     prev_node = '';
    for i = 1:length(Out_solName)    

        split_sol = Out_solName{i}(5:end);% strsplit(Out_solName{i},','); % splitting 1;A,1-1 to give cell(2,1) with cell{1} = 1;A and cell{2} = 1-1

        if(~ ismember('e',Out_solName{i}))
            add_node = sprintf('V%s',split_sol); % 1-1 becomes V1-1
        end

        if(~isequal(add_node, prev_node))  % checking if new node being added wasn't added in the previous step itself*********  % this was previous condition (~ismember(addnode, fin_sol))&& (~ ismember('e',Out_solName{i})

            fin_sol{sol_count} = add_node;

            sol_count = sol_count+1;
        end

        prev_node = add_node;

    end
    
    
    fin_sol = fin_sol(~cellfun('isempty',fin_sol)); % removing empty cells
%     fin_rm_redunt(:) = {''};
%     finsol_count = 1;
%     solnode = cell(length(fin_sol),1);
%     solclus = cell(length(fin_sol),1);
% 
%     for i = 1:length(fin_sol)
% 
%         if(i==1)
%             fin_rm_redunt{finsol_count}  = strtok(fin_sol{i}, '-');
%             cur_node = strtok(fin_rm_redunt{finsol_count},'V');
%             cur_clus = strsplit(fin_sol{i}, '-');
%             solclus{finsol_count} = cur_clus{2};
%             solnode{finsol_count} = cur_node;
%             finsol_count = finsol_count+1;
%         else
%              cur_node  = strtok(strtok(fin_sol{i}, '-'),'V');
%              cur_clus = strsplit(fin_sol{i}, '-');
%              solclus{i} = cur_clus{2};       
%              solnode{i} = cur_node;
%              if(~((cur_clus{2} == solclus{i-1}) || (cur_node == solnode{i-1})))
%                  fin_rm_redunt{finsol_count}  = strtok(fin_sol{i}, '-');             
%                  finsol_count = finsol_count+1;
%              end
% 
%         end
%     end

    
    





end

