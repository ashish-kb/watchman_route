function [ outfin_sol, outfin_cost, Out_solName, whole_path_nodes,G_init, G_nodebot, edges_totsp, nodes_totsp, time_concorde_struct] = solveMulti_GTSP_noon_bean(V_adj, V_Cluster, V_adj_bot) % the fomat is given above
    
     given_nodes = length(V_adj);
    node_name  = cell(1, given_nodes);
    
    for i = 1:given_nodes
               
        node_name{i} = sprintf('V%d',i);
    
    end
    
    

    G_init = digraph(V_adj, node_name);
 
%      figure;
%      P_init = plot(G_init,'EdgeLabel', G_init.Edges.Weight);
%      title('Given Graph');


    [sz_r_G , sz_c_G] = size(V_adj);
%%
%side track calculating distances of from bots

    nodebot_name = cell(1, length(V_adj_bot));
    num_bots = length((given_nodes+1):length(V_adj_bot));

    for i = 1:given_nodes
              
       nodebot_name{i} = sprintf('V%d',i);
   
    end
     
     
     for i = (given_nodes+1):length(V_adj_bot)
         nodebot_name{i} = sprintf('B%d',(i-given_nodes));
     end
     
     G_nodebot = digraph(V_adj_bot, nodebot_name);
     
     V_nodebot_comp = distances(G_nodebot, 1:length(V_adj_bot), 1:length(V_adj_bot));
     
     G_nodebot_comp = digraph(V_nodebot_comp, nodebot_name);
     
%check if this has all the data we need 
    
%%
    % completing the graph
     V_comp = V_adj; 
 
     V_comp =  distances(G_init, 1:sz_r_G,1:sz_c_G);


     G_comp = digraph(V_comp, node_name);
    
     G_alpha = graph(V_comp, node_name, 'upper','OmitSelfLoops');
     Cluster_to_node = arrayfun(@(i)find(cellfun(@(s)ismember(i,s), V_Cluster)), 1:max([V_Cluster{:}]) , 'UniformOutput', false); % reverse lookup for Cluster_cell %http://stackoverflow.com/questions/14934796/reverse-lookup-in-matlab-cell-array-of-indices
%     figure;
%     plot(G_comp, 'EdgeLabel', G_comp.Edges.Weight);
%     title('Completed Graph');

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
    alpha_noon = 2*sum(G_alpha_mst.Edges.Weight(:)); % check p1(noon bean paper) we want to add zero cost edges in the next step in p3 between (5b)and(5c) these edges should be preferred over a tour in p1 that is 2*weightMST
    % verfied these alpha and beta costs works
    alpha_noon_botadd = alpha_noon + 2*num_bots*max(G_nodebot_comp.Edges.Weight);
    beta_noon  = 2*alpha_noon_botadd*(1+length(G_alpha_mst.Edges.Weight(:))); % check edges in p4 pg 28 noon bean - the zero cost edges that will be added between (4b)and (3b) in p6 should be 
                                                                                         % prefered over a complete tour in p4 and a complete tour in p4 would not be bigger than alpha_noon*numberOfEdgesInG_alphaMST
    beta_noon_botadd = beta_noon + 2*num_bots*(1+alpha_noon_botadd);
    %****my addition ***********%
    %beta_noon_botadd = 2*beta_noon_botadd*length(G_alpha_mst.Edges.Weight(:)) + 2*num_bots*beta_noon_botadd;
%some plots in the plot_script
     G_comp.Nodes.Cluster = V_Cluster;
     
     ind_empty_clus = cellfun(@isempty, G_comp.Nodes.Cluster);  % this contains indices of graph where there is no cluster relationship means orphan nodes
     G_comp_tm = rmnode(G_comp, find(ind_empty_clus));
     
%%

    %add new nodes which belong to more than one cluster
    G_comp = G_comp_tm;
    IN_tab_Nodes = IN_transform_Nodes(G_comp); % intersecting to non intersecting clusters nodes

    G_comp_temp = digraph([], []);
%     G_comp_temp1 = digraph([], []);
   

    G_comp_temp = addnode(G_comp_temp, IN_tab_Nodes); 
    

% 
%      figure;
%      plot(G_comp, 'EdgeLabel', G_comp.Edges.Weight);
%      title('Given Completed Graph - without - duplicate nodes added and orphan nodes removed');
    
    given_nodes = G_comp.numnodes; 
    
    % add edges to the new nodes which are made in previous step

    [IN_Interedge_s, IN_Interedge_e, IN_Interedge_w] = IN_transform_Interedge(G_comp); % intersecting to non intersecting clusters edges
    
    
    IN_Interedge_s = IN_Interedge_s(~cellfun('isempty',IN_Interedge_s(:,1))); % remove the empty entries from the CELL
    IN_Interedge_e = IN_Interedge_e(~cellfun('isempty',IN_Interedge_e(:,1)));
    IN_Interedge_w = IN_Interedge_w(IN_Interedge_w>0);
    
    
%   IN_tab_Interedge = IN_tab_Interedge(~cellfun('isempty',IN_tab_Interedge.EndNodes(:,1)),:); % remove the empty entries from the table
    
%   G_comp_temp = addedge(G_comp_temp, IN_tab_Interedge);
    G_comp_temp = addedge(G_comp_temp, IN_Interedge_s, IN_Interedge_e, IN_Interedge_w);  % fin_sol = fin_sol(~cellfun('isempty',fin_sol)); % removing empty cells
    
%   adding zero cost edges between the new nodes

    IN_tab_Zeroedge = IN_transform_Zeroedge(G_comp);

    G_comp_temp = addedge(G_comp_temp, IN_tab_Zeroedge);
%     
%   figure;     
%   p_interedge= plot(G_comp_temp, 'EdgeLabel', G_comp_temp.Edges.Weight);
%   title('Interedges, duplicate nodes and zero cost edges');
    
%   pdata_right = [3.98084484112396,7.95533836363301;9.51039443389633,0.414757900564329;3.15039773815711,0.414757900564329;3.15722177054659,4.00932637212185;5.94136112806105,3.20478998863015;0.411811535288686,6.76618145828282; 5.94136112806098,0.414757900564329; 0.411811535288686,3.17981057099319]; 
    
    
%   figure;
%   plot(G_comp_temp, 'XData',  pdata_right(:,1), 'YData',  pdata_right(:,2), 'EdgeLabel', G_comp_temp.Edges.Weight);
%   title('I-N transformation ducplicate nodes');
    
    
    G_comp = G_comp_temp;
    
    [s_ t_] = findedge(G_comp);
    
    Adj_G_comp = full(adjacency(G_comp));
    Adj_G_comp_ind = sub2ind(size(Adj_G_comp), s_(:),t_(:));

    Adj_G_comp(Adj_G_comp_ind(:)) = G_comp.Edges.Weight(:);
    diag_ind = sub2ind(size(Adj_G_comp), 1:length(Adj_G_comp),1:length(Adj_G_comp));
    
   
    %convert to atsp code used from https://www.mathworks.com/matlabcentral/fileexchange/44467-noon-bean-transformation
    %have made significant changes to their implementation
    [atspAdjMatrix infcost]  = gtsp_to_atsp(Adj_G_comp, cell2mat(G_comp.Nodes.Cluster), alpha_noon_botadd, beta_noon_botadd, G_comp);
    
   % symTSP = atsp_to_tsp(atspAdjMatrix, infcost);
    %\todo make atspAdjMatrix diagonal =0 or omitself loops and I have made
    %0 cost edges to 0.01 cost, do something about that.
     [X_t Y_s] = meshgrid(1:length(G_comp.Nodes.Name), 1:length(G_comp.Nodes.Name));
     G_atsp = digraph(Y_s,X_t, atspAdjMatrix(:), G_comp.Nodes.Name,'OmitSelfLoops');  % the command below omits zero cost edges that is why it's done like this
     
%      figure;
     %[G_atsp_s G_atsp_t] = findedge(G_atsp);
%      p_noonbean = plot(G_atsp, 'XData',  p_interedge.XData, 'YData', p_interedge.YData, 'EdgeLabel', G_atsp.Edges.Weight);
%      highlight(p_noonbean, [2;3], [3;2], 'EdgeColor', 'k', 'LineStyle', ':', 'LineWidth', 1);
%      highlight(p_noonbean, [2;2;3;3], [1;4;1;4], 'EdgeColor', 'r', 'LineStyle', '-', 'LineWidth', 0.8);
%      title('tail shifted noon bean done');
%      axis equal;
     
%%
     % adding the bot nodes  to the G_atsp graph
 

     bot_Nodes =  table(cell(num_bots,1), cell(num_bots,1), 'VariableNames', {'B_d', 'B_f'}); % external edge tab
     
     depot_rep = repmat(1:num_bots, num_bots, 1);
     finish_rep = repmat([1:num_bots]', num_bots, 1);
     
     str_d = arrayfun(@(i) sprintf('B%d-d', i), depot_rep(:), 'UniformOutput', false);
     str_f = arrayfun(@(i) sprintf('B%d-f', i), finish_rep, 'UniformOutput', false);
         
     bot_Nodes = table(str_d, str_f, 'VariableNames', {'B_d', 'B_f'});




     % adding bot edges to the G_atsp graph     
     
     %******************the following logic is currently coded*******%
     % or make a table with all the edges we need B1_d to V1-1 and also
     % V1-1 and B1_f then search first 2 elements of these cells(cellfun(@(x) x(1:2) , G_atsp.Nodes.Name(:,1),'uni',0)) using
     % findedge in G_nodebot use this list to get the column of weights
     %
     %
    str_d_vec = cell(length(G_atsp.Nodes.Name(:))*num_bots, 1);
    str_V_vec = cell(length(G_atsp.Nodes.Name(:))*num_bots, 1);
    str_f_vec = cell(length(G_atsp.Nodes.Name(:))*num_bots, 1);
    weight_vec = zeros(length(G_atsp.Nodes.Name(:))*num_bots, 1);
    cell_jumpcount = 1; % this counter jumps from 0 to length(G_atsp.Nodes.Name(:)) and then 2*length(G_atsp.Nodes.Name(:)) and then fills up the array

    for i = 1:num_bots
       str_d_vec(cell_jumpcount:i*length(G_atsp.Nodes.Name(:))) = cellfun(@(x) sprintf('B%d-d', i), G_atsp.Nodes.Name(:),'uni', 0);
       str_f_vec(cell_jumpcount:i*length(G_atsp.Nodes.Name(:))) = cellfun(@(x) sprintf('B%d-f', i), G_atsp.Nodes.Name(:),'uni', 0);
       str_V_vec(cell_jumpcount:i*length(G_atsp.Nodes.Name(:))) = G_atsp.Nodes.Name(:);
       weight_vec(cell_jumpcount:i*length(G_atsp.Nodes.Name(:)))...
           = G_nodebot_comp.Edges.Weight(findedge(G_nodebot_comp, cellfun(@(x) sprintf('B%d', i), G_atsp.Nodes.Name(:),'uni', 0), cellfun(@(x) x(1:(regexp(x,'-','start')-1)) , G_atsp.Nodes.Name(:,1),'uni',0)));
       cell_jumpcount = cell_jumpcount + length(G_atsp.Nodes.Name(:));
    end
    
    
    
    G_atsp = addedge(G_atsp, str_d_vec, str_V_vec, [weight_vec+(alpha_noon_botadd+beta_noon_botadd)*ones(size(weight_vec))]); % departure nodes are added penalty
%   G_atsp = addedge(G_atsp, str_V_vec, str_f_vec, [weight_vec]); 
%  ******* incoming node without penalty % Adding tail switching code
    
    edge_table = G_atsp.Edges;
    zero_weight_ind = edge_table.Weight==0;

    cur_clus_cell = cell(1, length(Cluster_to_node));
    weight_vec_cell = cell(num_bots, 1);
    str_f_cell = cell(num_bots,1);

    for i = 1: length(Cluster_to_node)
        log_ind_clus(:,i) = cell2mat(cellfun(@(x) isequal(sprintf('-%d', i),x(regexp(x,'-','start'):end)), edge_table.EndNodes(:,2), 'uni', 0 ));
        cur_clus_ind = log_ind_clus(:,i) & zero_weight_ind;
        cur_clus_cell{1,i} = edge_table.EndNodes(cur_clus_ind==1, 1); % cell containing names of all the nodes in cluster i connected in a cycle
        
        for j = 1:num_bots
            
            weight_vec_cell{j,1} = ...
            circshift(G_nodebot_comp.Edges.Weight(findedge(G_nodebot_comp, cellfun(@(x) sprintf('B%d', j), cur_clus_cell{1,i},'uni', 0), cellfun(@(x) x(1:(regexp(x,'-','start')-1)) , cur_clus_cell{1,i},'uni',0))),-1); 
            
            str_f_cell{j, 1} = cellfun(@(x) sprintf('B%d-f', j), cur_clus_cell{1,i},'uni', 0);


        end

        concat_f_cell = [str_f_cell{:,:}];
        concat_w_cell = [weight_vec_cell{:,:}];
        G_atsp = addedge(G_atsp, repmat(cur_clus_cell{1,i}, num_bots, 1), concat_f_cell(:), [concat_w_cell(:)]); 

    end

% **********end of tail shifting code for edges from nodes to finish depots
   
   
   
    G_atsp = addedge(G_atsp, bot_Nodes.B_d, bot_Nodes.B_f, zeros(length(bot_Nodes.B_d), 1));
    G_atsp = addedge(G_atsp, bot_Nodes.B_f, bot_Nodes.B_d, zeros(length(bot_Nodes.B_d), 1));
    
%     figure;
    %plot(G_atsp, 'EdgeLabel', G_atsp.Edges.Weight);
%     p_final = plot(G_atsp, 'XData',  [p_interedge.XData 0 1 1.3 2], 'YData', [p_interedge.YData -2.5  -3.5 -2.5  -3.5], 'EdgeLabel', G_atsp.Edges.Weight);
%     highlight(p_final,[7 8],'NodeColor','k');
%     highlight(p_final,[5 6],'NodeColor','g');
%     highlight(p_final, [2;3;5;5;6;6;7;7;8;8], [3;2;7;8;7;8;5;6;5;6], 'EdgeColor', 'k', 'LineStyle', ':', 'LineWidth', 1);
%     highlight(p_final, [5;5;5;5;6;6;6;6], [1;2;3;4;1;2;3;4], 'EdgeColor', 'g', 'LineStyle', '-', 'LineWidth', 0.8);
%     highlight(p_final, [2;2;3;3], [1;4;1;4], 'EdgeColor', 'r', 'LineStyle', '-', 'LineWidth', 0.8);
%     
    
    
    
%     p_final.ArrowSize = 11
%     
%     axis equal;
%     title('bot edges added to the noonbean-ed graph');
    
    %%
    %preparing the directed tsp to normal tsp to send it to concorde
     
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
    %following is just a logic to extract the final solution from
    %split_Out, you can have a look at it and understand the answer or node
    %is same, 
    
    
    
    % cycle when first and last clusters are same
    % stops and when both first and last clus are different and nodes are different
    cyclic = 1; % assumed not cyclic no split 
    cur_cycle_checked = 1;
   

    while (cyclic == 1)
    
    % repeat_check = find(~cellfun(@isempty,strfind(temp_fin_sol, temp_fin_sol{1}))); %test if there is a repeat
          conc_nod_first = strtok(strtok(split_Out{1}, '-'),'V');  %connected nod first 
          conc_nod_last = strtok(strtok(split_Out{end}, '-'),'V');  % connected last last

          conc_clus_first = strsplit(split_Out{1}, '-');  %connected cluster first 
          conc_clus_last = strsplit(split_Out{end}, '-');  % connected cluster last
          if((isequal(conc_clus_first{end},conc_clus_last{end}))|| (isequal(conc_nod_first,conc_nod_last)))
                            
             
             
             if ((cur_cycle_checked > length(split_Out)) && (isequal(conc_nod_first,conc_nod_last)) && (~isequal(conc_clus_first{2},conc_clus_last{2})))
                cyclic = 0;  % if first and last node are same and we have cycled more than once this happens very rarely - i remember when there are less than 3 orig nodes 
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
    %extracting individual robot veritices 
    
    bot_path = cell(length(split_Out), num_bots);

    for i=1:num_bots

       b_start = sprintf('B%d-d', i);
       %b_end = sprintf('B%d-f', i);
       b_ind_start =  find(cell2mat(cellfun(@(x) ismember(x, b_start), {split_Out'},'uni',0)));
      % b_ind_end =  find(cell2mat(cellfun(@(x) ismember(x, b_end), {split_Out'},'uni',0)));
       b_ind_next = b_ind_start + 1;
       b_ind_prev = b_ind_start - 1;
       if(b_ind_next > length(split_Out))
           b_ind_next = 1;
       end
       if(b_ind_prev < 1)
           b_ind_prev = length(split_Out);
       end

       if (isequal(split_Out{b_ind_next}(1), 'V'))
            split_circ = split_Out';                 
            split_circ = circshift(split_circ, -(b_ind_start-1));
            depot_f_vec = find(cell2mat(cellfun(@(x) ismember(x(end), 'f'), split_circ, 'uni',0))); % B1-f type of nodes i.e. finish nodes 
            b_ind_end = depot_f_vec(1);
            b_end = split_circ(b_ind_end);
            
           % b_ind_end =  find(cell2mat(cellfun(@(x) ismember(x, b_end), {split_circ},'uni',0)));
           
            bot_path(1:b_ind_end, i) = split_circ(1:b_ind_end)';
            
       elseif(isequal(split_Out{b_ind_prev}(1), 'V'))
           split_circ = flip(split_Out');
           b_ind_start = length(split_Out) - b_ind_start + 1;           
           split_circ = circshift(split_circ, -(b_ind_start-1));
           depot_f_vec =  find(cell2mat(cellfun(@(x) ismember(x(end), 'f'), split_circ, 'uni',0))); % B1-f type of nodes i.e. finish nodes 
           b_ind_end = depot_f_vec(1);
           b_end = split_circ(b_ind_end);
%            b_ind_end = length(split_Out) - b_ind_end + 1;
%            b_ind_end =  find(cell2mat(cellfun(@(x) ismember(x, b_end), {split_circ},'uni',0)));
          
           bot_path(1:b_ind_end, i) = split_circ(1:b_ind_end)';
          
       else
           % do nothing if both sides don't have a node     
       end


    end
    
    % extracting noden_num for all the bots
    
    node_num_forw = [];
    test_prev_forw = 'o';
    cnter = 2;  % counter
    
    node_num_forw_cell = cell(num_bots,1);
    % checking forward set for solution as we can have solution in both
    % direction of split_Out
    
   for j = 1:num_bots 
        
        bot_path_trim = bot_path(:,j);
        bot_path_trim = bot_path_trim(~cellfun('isempty',bot_path_trim));
        
        if(~isempty(bot_path_trim))
       
            for i = 2:(length(bot_path_trim)-1) % random initialisation
                if(isequal(test_prev_forw, bot_path_trim{i}(regexp(bot_path_trim{i},'-','start'):end)))
                    %seq_forw(i) = false;
                else
                    node_num_forw(cnter) =  str2num(bot_path_trim{i}(2:(regexp(bot_path_trim{i},'-','start')-1)));
                    cnter = cnter + 1;
                end   

                test_prev_forw = bot_path_trim{i}(regexp(bot_path_trim{i},'-','start'):end);

            end

            node_num_forw(1) = -1*str2num(bot_path_trim{1}(2:(regexp(bot_path_trim{1},'-','start')-1)));
            node_num_forw(cnter) = -1*str2num(bot_path_trim{end}(2:(regexp(bot_path_trim{end},'-','start')-1)));

            node_num_forw_cell{j,1} = node_num_forw; %unique(node_num_forw, 'stable');
            node_num_forw = [];
            test_prev_forw = 'o';
            cnter = 2;
        end
   end
    
   % will make a graph out of these to figure out the cost
   forw_cost = zeros(num_bots, 1);
   
    for j = 1:num_bots
        if (~isempty(node_num_forw_cell{j, 1}))        
           s_forw(1:(length(node_num_forw_cell{j, 1})-3)) = node_num_forw_cell{j, 1}(2:(end-2));
           t_forw(1:(length(node_num_forw_cell{j, 1})-3)) = node_num_forw_cell{j, 1}(3:end-1);
           %s_forw(length(node_num_forw_cell{j, 1})-1) = ;
           %t_forw(length(node_num_forw_cell{j, 1})-1) = ;

           %s_forw(1) = ;
           %t_forw(1) = ;
           %for edges betwwen nodes
           G_comp_allnode = digraph(V_comp, node_name);

           edges_forw = findedge(G_comp_allnode, s_forw,t_forw);

           %for edges between bots to nodes and nodes to bots
           str_s_bot2node = sprintf('B%d', -1*node_num_forw_cell{j,1}(1));
           str_t_bot2node = sprintf('V%d', node_num_forw_cell{j,1}(2));
           edge_bot2node = findedge(G_nodebot_comp, str_s_bot2node,str_t_bot2node);
           cost_edge_bot2node = G_nodebot_comp.Edges.Weight(edge_bot2node);


           str_s_node2bot = sprintf('V%d', node_num_forw_cell{j,1}(end-1));
           str_t_node2bot = sprintf('B%d', -1*node_num_forw_cell{j,1}(end));
           edge_node2bot = findedge(G_nodebot_comp, str_s_node2bot,str_t_node2bot);
           cost_edge_node2bot = G_nodebot_comp.Edges.Weight(edge_node2bot);


           forw_cost(j,1) = sum(G_comp_allnode.Edges.Weight(edges_forw(edges_forw~=0))) + cost_edge_bot2node + cost_edge_node2bot;


           s_forw = [];
           t_forw = [];
        end
    end
   
    outfin_cost = forw_cost;
    outfin_sol = node_num_forw_cell;
    
    whole_path_nodes = cell(num_bots,1);

    for i = 1:num_bots
        for j = 1:length(node_num_forw_cell{i})
           if(j==1)
               str_s_bot2node = sprintf('B%d', -1*node_num_forw_cell{i,1}(1));
               str_t_bot2node = sprintf('V%d', node_num_forw_cell{i,1}(2));
               whole_path_nodes{i} =  [whole_path_nodes{i}; shortestpath(G_nodebot, str_s_bot2node, str_t_bot2node)];  
           elseif(j>1 && j< length(node_num_forw_cell{i})-1)
               str_s_node2node = sprintf('V%d', node_num_forw_cell{i,1}(j));
               str_t_node2node = sprintf('V%d', node_num_forw_cell{i,1}(j+1));
               whole_path_nodes{i} =  [whole_path_nodes{i}(1:end-1) shortestpath(G_nodebot, str_s_node2node, str_t_node2node)]; 
           elseif (j==length(node_num_forw_cell{i})-1)

               str_s_node2bot = sprintf('V%d', node_num_forw_cell{i,1}(end-1));
               str_t_node2bot = sprintf('B%d', -1*node_num_forw_cell{i,1}(end));
               whole_path_nodes{i} =  [whole_path_nodes{i}(1:end-1) shortestpath(G_nodebot, str_s_node2bot, str_t_node2bot)]; 
           end



        end
    end
    % 
    
 
