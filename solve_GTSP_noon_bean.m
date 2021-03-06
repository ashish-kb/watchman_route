function [ outfin_sol, outfin_cost,Out_solName, Out_sol, G_init, edges_totsp, nodes_totsp, time_concorde_struct] = solve_GTSP_noon_bean(V_adj, V_Cluster) % the fomat is given above

    given_nodes = length(V_adj);
    node_name  = cell(1, given_nodes);
    
    for i = 1:given_nodes
               
        node_name{i} = sprintf('V%d',i);
    
    end
    
    

    G_init = digraph(V_adj, node_name);
 
%     figure;
%     P_init = plot(G_init,'EdgeLabel', G_init.Edges.Weight);
%     title('Given Graph');


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
    
     G_alpha = graph(V_comp, node_name, 'upper','OmitSelfLoops');
     Cluster_to_node = arrayfun(@(i)find(cellfun(@(s)ismember(i,s), V_Cluster)), 1:max([V_Cluster{:}]) , 'UniformOutput', false); % reverse lookup for Cluster_cell %http://stackoverflow.com/questions/14934796/reverse-lookup-in-matlab-cell-array-of-indices
%       figure;
%       plot(G_comp, 'EdgeLabel', G_comp.Edges.Weight);
%       title('Completed Graph');
% extracting selective nodes to get a tour which would satisfy alpha and
% beta requirement for noon bean. instead of sum of edges be upperbound this
% will be the upperbound. And below we look for a tour to just any visit the
% node in a cluster but cover all clusters(set cover not optimal) this tour
% though valid but will cost more than optimal. hence we use it as alpha
     node_clus_mat = zeros(length(V_Cluster), length(Cluster_to_node)); % col is clus and row is node
     for i = 1:length(V_Cluster)
        node_clus_mat(i,V_Cluster{i}) = 1;
     end
     
     keep_node = [];
     for i = 1:length(Cluster_to_node)
        
        temp_keep = find(node_clus_mat(:,i)==1);
        keep_node = [keep_node;temp_keep(1)];
     end
     keep_node =  unique(keep_node);
    
    nod_num = 1:length(V_Cluster);
    rm_node = nod_num(~ismember(nod_num, keep_node));
    
    G_alpha = rmnode(G_alpha, rm_node);
    G_alpha_mst = minspantree(G_alpha);
    alpha_noon = 2*sum(G_alpha_mst.Edges.Weight(:)); % check p1 we want to add zero cost edges in the next step in p3 between (5b)and(5c) these edges should be preferred over a tour in p1 that is 2*weightMST
    beta_noon  = 2*2*sum(G_alpha_mst.Edges.Weight(:))*length(G_alpha_mst.Edges.Weight(:));% check edges in p4 pg 28 noon bean - the zero cost edges that will be added between (4b)and (3b) in p6 should be 
                                                                                        % prefered over a complete tour in p4 and a complete tour in p4 would not be bigger than alpha_noon*numberOfEdgesInG_alphaMST

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
    


%     figure;
%     plot(G_comp);
%     title('Given Completed Graph - without - duplicate nodes added and orphan nodes removed');
    
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
    
%     figure;
%     title('Interedges, duplicate nodes and zero cost edges');
%     plot(G_comp_temp, 'EdgeLabel', G_comp_temp.Edges.Weight);
    
%    pdata_right = [3.98084484112396,7.95533836363301;9.51039443389633,0.414757900564329;3.15039773815711,0.414757900564329;3.15722177054659,4.00932637212185;5.94136112806105,3.20478998863015;0.411811535288686,6.76618145828282; 5.94136112806098,0.414757900564329; 0.411811535288686,3.17981057099319]; 
    
    
%     figure;
%     plot(G_comp_temp, 'XData',  pdata_right(:,1), 'YData',  pdata_right(:,2), 'EdgeLabel', G_comp_temp.Edges.Weight);
%     title('I-N transformation ducplicate nodes');
    
    
    G_comp = G_comp_temp;
    
    [s_ t_] = findedge(G_comp);
    
    Adj_G_comp = full(adjacency(G_comp));
    Adj_G_comp_ind = sub2ind(size(Adj_G_comp), s_(:),t_(:));

    Adj_G_comp(Adj_G_comp_ind(:)) = G_comp.Edges.Weight(:);
    diag_ind = sub2ind(size(Adj_G_comp), 1:length(Adj_G_comp),1:length(Adj_G_comp));
    
   
    %convert to atsp 
    [atspAdjMatrix infcost]  = gtsp_to_atsp(Adj_G_comp, cell2mat(G_comp.Nodes.Cluster), alpha_noon, beta_noon, G_comp);
    
   % symTSP = atsp_to_tsp(atspAdjMatrix, infcost);
    %\todo make atspAdjMatrix diagonal =0 or omitself loops and I have made
    %0 cost edges to 0.01 cost, do something about that.
     [X_t Y_s] = meshgrid(1:length(G_comp.Nodes.Name), 1:length(G_comp.Nodes.Name));
     G_atsp = digraph(Y_s,X_t, atspAdjMatrix(:), G_comp.Nodes.Name,'OmitSelfLoops');  % the command below omits zero cost edges that is why it's done like this
    %G_atsp = digraph(atspAdjMatrix, G_comp.Nodes.Name);
     str_1node = cellfun(@(x,y) sprintf('1;%s', x), G_atsp.Edges.EndNodes(:,1),'uni', 0);
     str_3node = cellfun(@(x,y) sprintf('3;%s', x), G_atsp.Edges.EndNodes(:,2),'uni', 0);
    
     G_atsp2_tsp = graph([], []);
     
     str_1node_intra = cellfun(@(x,y) sprintf('1;%s', x), G_atsp.Nodes.Name,'uni', 0);
     str_2node_intra = cellfun(@(x,y) sprintf('2;%s', x), G_atsp.Nodes.Name,'uni', 0);
     str_3node_intra = cellfun(@(x,y) sprintf('3;%s', x), G_atsp.Nodes.Name,'uni', 0);
     
     G_atsp2_tsp = addedge(G_atsp2_tsp, [str_1node_intra;str_2node_intra], [str_2node_intra;str_3node_intra], zeros(2*G_atsp.numnodes,1));

     str_1node = cellfun(@(x,y) sprintf('1;%s', x), G_atsp.Edges.EndNodes(:,1),'uni', 0);
     str_3node = cellfun(@(x,y) sprintf('3;%s', x), G_atsp.Edges.EndNodes(:,2),'uni', 0);
     
     G_atsp2_tsp = addedge(G_atsp2_tsp, str_1node,str_3node,G_atsp.Edges.Weight(:)); %
     nodes_totsp = G_atsp2_tsp.numnodes;
     edges_totsp = G_atsp2_tsp.numedges;
    
     [Out_sol, time_concorde_struct] = TSP_tour_Dat(G_atsp2_tsp,'/home/ashishkb/softwares/concorde/concorde/TSP/concorde');
    
    Out_solName = G_atsp2_tsp.Nodes.Name(Out_sol);
    
    %vertexSequenceOrdered = get_concorde_result(symTSP,cell2mat(G_comp.Nodes.Cluster));
    
    Out_solName = Out_solName(cell2mat(cellfun(@(x) ismember('1',x(1)) , Out_solName,'uni',0)));
    cnter = 1;
    for i =1:length(Out_solName) 
        split_Out{cnter}= Out_solName{i}(3:end);
        cnter  = cnter+1;
    end
    
    %%
    % cycle so that first and last clusters are same
    
    cyclic = 1; % assumed not cyclic no split 
    cur_cycle_checked = 1;
   

    while (cyclic == 1)
    
    % repeat_check = find(~cellfun(@isempty,strfind(temp_fin_sol, temp_fin_sol{1}))); %test if there is a repeat
          conc_nod_first = strtok(strtok(split_Out{1}, '-'),'V');  %connected nod first 
          conc_nod_last = strtok(strtok(split_Out{end}, '-'),'V');  % connected last last

          conc_clus_first = strsplit(split_Out{1}, '-');  %connected cluster first 
          conc_clus_last = strsplit(split_Out{end}, '-');  % connected cluster last
          if((isequal(conc_clus_first{2},conc_clus_last{2}))|| (isequal(conc_nod_first,conc_nod_last)))
                            
             
             
             if ((cur_cycle_checked > length(split_Out)) && (isequal(conc_nod_first,conc_nod_last)) && (~isequal(conc_clus_first{2},conc_clus_last{2})))
                cyclic = 0;  % if first and last node are same and we have cycled more than once
             else
                 split_Out = circshift(split_Out,1,2);%cycle me 
                 cyclic = 1;
             end
          else
               cyclic = 0; % it is cyclic already
          end
          cur_cycle_checked = cur_cycle_checked+1;
          
          
          
    end
    
    
    
    
    
    
    
    
    
    
    
    %%
    
    
    
   % split_Out{cnter} = split_Out{1};
    %seq_clus = [1,3,5,2,4,4,6,6];
    %seq_forw = true(1,length(seq_clus));
    node_num_forw = [];
    test_prev_forw = 'o';
    cnter = 1;
    % checking forward set for solution as we can have solution in both
    % direction of split_Out
   for i = 1:length(split_Out) % random initialisation
        if(isequal(test_prev_forw, split_Out{i}(regexp(split_Out{i},'-','start'):end)))
            %seq_forw(i) = false;
        else
            node_num_forw(cnter) =  str2num(split_Out{i}(2:(regexp(split_Out{i},'-','start')-1)));
            cnter = cnter + 1;
        end   
        
        test_prev_forw = split_Out{i}(regexp(split_Out{i},'-','start'):end);
        
   end
    % checking for backward set
   cnter = 1;
   node_num_back = [];
   test_prev_back = 't'; % random initialisation
   for i = length(split_Out):-1:1
        if(isequal(test_prev_back,split_Out{i}(regexp(split_Out{i},'-','start'):end)))
            %seq_forw(i) = false;
        else
            node_num_back(cnter) = str2num( split_Out{i}(2:(regexp(split_Out{i},'-','start')-1)));
            cnter = cnter + 1;
        end   
        
        test_prev_back = split_Out{i}(regexp(split_Out{i},'-','start'):end);
        
   end
   
   
   
   s_forw(1:(length(node_num_forw)-1)) = node_num_forw(1:(end-1));
   t_forw(1:(length(node_num_forw)-1)) = node_num_forw(2:end);
   s_forw(length(node_num_forw)) = node_num_forw(end);
   t_forw(length(node_num_forw)) = node_num_forw(1);
   
   G_comp_allnode = digraph(V_comp, node_name);
   
   edges_forw = findedge(G_comp_allnode, s_forw,t_forw);
   forw_cost = sum(G_comp_allnode.Edges.Weight(edges_forw(edges_forw~=0)));
   
   s_back(1:(length(node_num_back)-1)) = node_num_back(1:(end-1));
   t_back(1:(length(node_num_back)-1)) = node_num_back(2:end);
   s_back(length(node_num_back)) = node_num_back(end);
   t_back(length(node_num_back)) = node_num_back(1);
   
   edges_back = findedge(G_comp_allnode, s_back,t_back);
   back_cost = sum(G_comp_allnode.Edges.Weight(edges_back(edges_back~=0)));
    
   
   split_Out_back = flip(split_Out);
   
   if(forw_cost<=back_cost)
        outfin_sol = node_num_forw;
        outfin_cost = forw_cost;
   else
        outfin_sol = node_num_back;
        outfin_cost = back_cost;
   end
    %fin_rm_redunt = fin_rm_redunt(~cellfun('isempty',fin_rm_redunt)); % removing empty cells
%     fin_rm_redunt = {'V20-9','V19-4','V5-6','V18-8','V14-5','V12-10','V17-7','V16-2','V8-3','V2-1'};%split_Out_back;
%     whole_path = {};
%     total_cost = 0;
%     for i = 1:(length(fin_rm_redunt)-1)
% 
%        %  cur_path = shortestpath(G_init, fin_rm_redunt(i), fin_rm_redunt(i+1)); % path between currently selected nodes
%          total_cost = total_cost + distances(G_comp, findnode(G_comp,fin_rm_redunt{i}), findnode(G_comp,fin_rm_redunt{i+1}));
%          if(i==1)
%             % whole_path = G_init.Nodes.Name{fin_rm_redunt(i)};
% 
%          end
%        %  whole_path =  [whole_path;  G_init.Nodes.Name{cur_path(2:end)'}]; %adding cur_path to whole_path
% 
%     end
%     
%     if( length(fin_rm_redunt) ==1) % won't go into previous loop as length is 1
%       % whole_path = fin_rm_redunt(1);
%     end
% 
%       %  total_cost = total_cost + distances(G_init, findnode(G_init,whole_path{end}), findnode(G_init,whole_path{1}));
       % end_path = shortestpath(G_init, whole_path{end}, whole_path{1}); % path between last node and first node
      %  whole_path = [whole_path; end_path(2:end)'];
    
    
end