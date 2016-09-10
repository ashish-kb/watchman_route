function [ outfin_sol, outfin_cost, Out_solName, G_init, G_nodebot, edges_totsp, nodes_totsp, time_concorde_struct] = solveMulti_GTSP_noon_bean(V_adj, V_Cluster, V_adj_bot) % the fomat is given above
    
    V_adj  = 10*V_adj;
    V_adj_bot = 10*V_adj_bot;
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
    alpha_noon = 2*sum(G_alpha_mst.Edges.Weight(:)); % check p1(noon bean paper) we want to add zero cost edges in the next step in p3 between (5b)and(5c) these edges should be preferred over a tour in p1 that is 2*weightMST
    %major additions that are a bit doubtful
    alpha_noon_botadd = alpha_noon + 2*num_bots*max(G_nodebot_comp.Edges.Weight);
    beta_noon  = 2*alpha_noon_botadd*length(G_alpha_mst.Edges.Weight(:)); % check edges in p4 pg 28 noon bean - the zero cost edges that will be added between (4b)and (3b) in p6 should be 
                                                                                         % prefered over a complete tour in p4 and a complete tour in p4 would not be bigger than alpha_noon*numberOfEdgesInG_alphaMST
    beta_noon_botadd = beta_noon + 2*num_bots*alpha_noon_botadd;
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
    


%     figure;
%     plot(G_comp);
%     title('Given Completed Graph - without - duplicate nodes added and orphan nodes removed');
    
    given_nodes = G_comp.numnodes; 
    
    % add edges to the new nodes which are made in previous step

    [IN_Interedge_s, IN_Interedge_e, IN_Interedge_w] = IN_transform_Interedge(G_comp); % intersecting to non intersecting clusters edges
    
    
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
    [atspAdjMatrix infcost]  = gtsp_to_atsp(Adj_G_comp, cell2mat(G_comp.Nodes.Cluster), alpha_noon_botadd, beta_noon_botadd, G_comp);
    
   % symTSP = atsp_to_tsp(atspAdjMatrix, infcost);
    %\todo make atspAdjMatrix diagonal =0 or omitself loops and I have made
    %0 cost edges to 0.01 cost, do something about that.
     [X_t Y_s] = meshgrid(1:length(G_comp.Nodes.Name), 1:length(G_comp.Nodes.Name));
     G_atsp = digraph(Y_s,X_t, atspAdjMatrix(:), G_comp.Nodes.Name,'OmitSelfLoops');  % the command below omits zero cost edges that is why it's done like this
     %%
     % adding the bot nodes  to the G_atsp graph
     
     % (a) add node b1 and b2 (b) look for V1 in G_atsp.Nodes and then look
     % for b1 to v1 in G_nodebot.Edges....make 2 edges outgoing and
     % incoming with the right costs... then .....then do the same with b2
     % to v1 and so on.
     
     bot_Nodes =  table(cell(num_bots,1), cell(num_bots,1), 'VariableNames', {'B_d', 'B_f'}); % external edge tab
     
     depot_rep = repmat(1:num_bots, num_bots, 1);
     finish_rep = repmat([1:num_bots]', num_bots, 1);
     
     str_d = arrayfun(@(i) sprintf('B%d-d', i), depot_rep(:), 'UniformOutput', false);
     str_f = arrayfun(@(i) sprintf('B%d-f', i), finish_rep, 'UniformOutput', false);
         
     bot_Nodes = table(str_d, str_f, 'VariableNames', {'B_d', 'B_f'});
                                 
         
%      str_d_comp_from = 0;mat2cell(cell2mat(cellfun(@(x) repmat(x, length(bot_Nodes.B_d),1), bot_Nodes.B_d(:),'uni',0)), ones(1,length(bot_Nodes.B_d)^2), length(bot_Nodes.B_d{}));
%      str_f_comp_to   = remat(bot_Nodes.finish, length(bot_Nodes.B_f),1);
     
     
     % adding bot edges to the G_atsp graph
     
     
     % index of B1 edges 
     %cell2mat(cellfun(@(x) ismember('B1',{x(1:2)}) , G_nodebot_comp.Edges.EndNodes(:,1),'uni',0))
     %2. cell2mat(cellfun(@(x,y) ismember('B1',{x(1:2)})&ismember('V1',{y(1:2)}), G_nodebot_comp.Edges.EndNodes(:,1), G_nodebot_comp.Edges.EndNodes(:,2),'uni',0))
     % 3. G_atsp.Nodes.Name(find(cell2mat(cellfun(@(x) ismember('V6',{x(1:2)}) , G_atsp.Nodes.Name(:,1),'uni',0))))
     %findedge(G_nodebot_comp, 'B1', G_nodebot_comp.Nodes.Name)
     % G_atsp.Nodes.Name{1,:}(1:2)
     
     % logic
     % from G_nodebot_comp get the entries corresponding to B1 then add
     % start search for 'V1' in G_atsp using (2.) above then use the weight
     % to make new edges str_fv V1-1 (3.) -> str_fb B1_f and str_db B1_d to str_dv V1-1 then 
     % finally add these edges to the graph 
     
     % or in G_atsp start making nodes between B1_d to all named nodes and
     % search for the associated weight from (2.) add beta1+beta2 to nodes
     % from depots
     
     
     
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
       weight_vec(cell_jumpcount:i*length(G_atsp.Nodes.Name(:))) = G_nodebot_comp.Edges.Weight(findedge(G_nodebot_comp, cellfun(@(x) sprintf('B%d', i), G_atsp.Nodes.Name(:),'uni', 0), cellfun(@(x) x(1:(regexp(x,'-','start')-1)) , G_atsp.Nodes.Name(:,1),'uni',0)));
       cell_jumpcount = cell_jumpcount + length(G_atsp.Nodes.Name(:));
    end
    
    
    
    G_atsp = addedge(G_atsp, str_d_vec, str_V_vec, [weight_vec+(alpha_noon_botadd+beta_noon_botadd)*ones(size(weight_vec))]); % departure nodes are added penalty
    G_atsp = addedge(G_atsp, str_V_vec, str_f_vec, [weight_vec*0]); % incoming node trying with and without penalty penalty = beta_noon_botadd*ones(size(weight_vec)) + alpha_noon_botadd*ones(size(weight_vec))
    G_atsp = addedge(G_atsp, bot_Nodes.B_d, bot_Nodes.B_f, zeros(length(bot_Nodes.B_d), 1));
    G_atsp = addedge(G_atsp, bot_Nodes.B_f, bot_Nodes.B_d, zeros(length(bot_Nodes.B_d), 1));
    %bot_order_changed = {'B2-d';'B2-f';'B1-d';'B1-f';'B3-d';'B3-f'};
  %  bot_order_changed = {'B2-d';'B2-f';'B3-d';'B3-f';'B1-d';'B1-f'};
   % G_atsp = addedge(G_atsp,  bot_order_changed, circshift(bot_order_changed, -1), zeros(length(bot_Nodes.Name), 1));

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

       if (isequal(split_Out{b_ind_start+1}(1), 'V'))
            split_circ = split_Out';                 
            split_circ = circshift(split_circ, -(b_ind_start-1));
            depot_f_vec = find(cell2mat(cellfun(@(x) ismember(x(4), 'f'), split_circ, 'uni',0))); % B1-f type of nodes i.e. finish nodes 
            b_ind_end = depot_f_vec(1);
            b_end = split_circ(b_ind_end);
            
           % b_ind_end =  find(cell2mat(cellfun(@(x) ismember(x, b_end), {split_circ},'uni',0)));
           
            bot_path(1:b_ind_end, i) = split_circ(1:b_ind_end)';
            
       elseif(isequal(split_Out{b_ind_start-1}(1), 'V'))
           split_circ = flip(split_Out');
           b_ind_start = length(split_Out) - b_ind_start + 1;           
           split_circ = circshift(split_circ, -(b_ind_start-1));
           depot_f_vec =  find(cell2mat(cellfun(@(x) ismember(x(4), 'f'), split_circ, 'uni',0))); % B1-f type of nodes i.e. finish nodes 
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
    
    % 
    
    % a hand made tour which costs less than the optimal here 
    % check the costs in G_nodebot_comp.Edges
    % b2 -> 9(7) -> 11(4,8) -> 13(5,6,7) -> b2 ::: 78.372814269112240 + 1.103684171840059e+02 + 1.425669603284892e+02 + 2.277035256767233e+02 
    %  
    % b3 -> v3(2, 1) -> v1(3) -> b3 ::: 1.275029275866113e+02 + 1.650546691729253e+02 + 40.933535537731274
    %  
    % b1 -> v18(19) -> v19(10*****) -> v22(11) -> b1   :::  85.316768751418860 + 58.366053168156080 + 1.508696283011450e+02 + 1.897923321106173e+02
    
    % tsp of the above checking only the actual visited nodes and all the
    % clusters they visit not checking intracluster edges
    % 3;2;1b2-d -> 3;2;1;V9-7 -> 3:2:1V11-4 -> 3;2;1;V11-8 -> 3:2:1V13-5-> 3;2;1V13-6 ->      ..
    
    
    % mostly the problem is with these edges -> {'1;B3-d','3;V1-3';'1;B3-d','3;V11-4'}
    
    
%     %%
%     
%     
%     
%    % split_Out{cnter} = split_Out{1};
%     %seq_clus = [1,3,5,2,4,4,6,6];
%     %seq_forw = true(1,length(seq_clus));
%     node_num_forw = [];
%     test_prev_forw = 'o';
%     cnter = 1;  % counter
%     % checking forward set for solution as we can have solution in both
%     % direction of split_Out
%    for i = 1:length(split_Out) % random initialisation
%         if(isequal(test_prev_forw, split_Out{i}(regexp(split_Out{i},'-','start'):end)))
%             %seq_forw(i) = false;
%         else
%             node_num_forw(cnter) =  str2num(split_Out{i}(2:(regexp(split_Out{i},'-','start')-1)));
%             cnter = cnter + 1;
%         end   
%         
%         test_prev_forw = split_Out{i}(regexp(split_Out{i},'-','start'):end);
%         
%    end
%     % checking for backward set
%    cnter = 1;
%    node_num_back = [];
%    test_prev_back = 't'; % random initialisation
%    for i = length(split_Out):-1:1
%         if(isequal(test_prev_back,split_Out{i}(regexp(split_Out{i},'-','start'):end)))
%             %seq_forw(i) = false;
%         else
%             node_num_back(cnter) = str2num( split_Out{i}(2:(regexp(split_Out{i},'-','start')-1)));
%             cnter = cnter + 1;
%         end   
%         
%         test_prev_back = split_Out{i}(regexp(split_Out{i},'-','start'):end);
%         
%    end
%    
%    
%    
%    s_forw(1:(length(node_num_forw)-1)) = node_num_forw(1:(end-1));
%    t_forw(1:(length(node_num_forw)-1)) = node_num_forw(2:end);
%    s_forw(length(node_num_forw)) = node_num_forw(end);
%    t_forw(length(node_num_forw)) = node_num_forw(1);
%    
%    G_comp_allnode = digraph(V_comp, node_name);
%    
%    edges_forw = findedge(G_comp_allnode, s_forw,t_forw);
%    forw_cost = sum(G_comp_allnode.Edges.Weight(edges_forw(edges_forw~=0)));
%    
%    s_back(1:(length(node_num_back)-1)) = node_num_back(1:(end-1));
%    t_back(1:(length(node_num_back)-1)) = node_num_back(2:end);
%    s_back(length(node_num_back)) = node_num_back(end);
%    t_back(length(node_num_back)) = node_num_back(1);
%    
%    edges_back = findedge(G_comp_allnode, s_back,t_back);
%    back_cost = sum(G_comp_allnode.Edges.Weight(edges_back(edges_back~=0)));
%     
%    
%    split_Out_back = flip(split_Out);
%    
%    if(forw_cost<=back_cost)
%         outfin_sol = node_num_forw;
%         outfin_cost = forw_cost;
%    else
%         outfin_sol = node_num_back;
%         outfin_cost = back_cost;
%    end
%     %fin_rm_redunt = fin_rm_redunt(~cellfun('isempty',fin_rm_redunt)); % removing empty cells
% %     fin_rm_redunt = {'V20-9','V19-4','V5-6','V18-8','V14-5','V12-10','V17-7','V16-2','V8-3','V2-1'};%split_Out_back;
% %     whole_path = {};
% %     total_cost = 0;
% %     for i = 1:(length(fin_rm_redunt)-1)
% % 
% %        %  cur_path = shortestpath(G_init, fin_rm_redunt(i), fin_rm_redunt(i+1)); % path between currently selected nodes
% %          total_cost = total_cost + distances(G_comp, findnode(G_comp,fin_rm_redunt{i}), findnode(G_comp,fin_rm_redunt{i+1}));
% %          if(i==1)
% %             % whole_path = G_init.Nodes.Name{fin_rm_redunt(i)};
% % 
% %          end
% %        %  whole_path =  [whole_path;  G_init.Nodes.Name{cur_path(2:end)'}]; %adding cur_path to whole_path
% % 
% %     end
% %     
% %     if( length(fin_rm_redunt) ==1) % won't go into previous loop as length is 1
% %       % whole_path = fin_rm_redunt(1);
% %     end
% % 
% %       %  total_cost = total_cost + distances(G_init, findnode(G_init,whole_path{end}), findnode(G_init,whole_path{1}));
%        % end_path = shortestpath(G_init, whole_path{end}, whole_path{1}); % path between last node and first node
%       %  whole_path = [whole_path; end_path(2:end)'];
%     
    
end