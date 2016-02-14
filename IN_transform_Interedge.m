function [IN_Interedge_s, IN_Interedge_e, IN_Interedge_w]= IN_transform_Interedge(graph) % the graph should contain Graph.Nodes.Cluster values....outputs a table - its part of IN transformation
    
    tot_clus = cellfun(@numel, graph.Nodes.Cluster); % matrix containing total of each cluster

    int_clus_tot = 0; % number of edges in the new graph

    for i = 1:length(tot_clus) %  for a,b,c = a(a+b+c-a) + b(a+b+c-b) + c(a+b+c-c) = total number of Intercluster edges after addition of new nodes

        int_clus_tot = int_clus_tot + tot_clus(i)*(sum(tot_clus)-tot_clus(i));
    end
  
    %table(cell(int_clus_tot,2), zeros(int_clus_tot,1), 'VariableNames', {'EndNodes','Weight'}); 
    IN_Interedge_s = cell(int_clus_tot,1);
    IN_Interedge_e = cell(int_clus_tot,1);
    IN_Interedge_w = zeros(int_clus_tot,1);
    
    
    edge_count = 1;
    %IN_tab_Interedge = table({}, [], 'VariableNames', {'EndNodes','Weight'});
    
    start_idx =  findnode(graph, graph.Edges.EndNodes(:,1));
    end_idx   =  findnode(graph, graph.Edges.EndNodes(:,2));
    
    clus_start_vec = graph.Nodes.Cluster(start_idx(:));
    clus_end_vec = graph.Nodes.Cluster(end_idx(:));
    
    node_nam_start_vec = graph.Nodes.Name(start_idx(:));
    node_nam_end_vec = graph.Nodes.Name(end_idx(:));
    
    for i = 1:graph.numedges
    %     start_idx = findnode(graph, graph.Edges.EndNodes{i,1});
    %     end_idx   = findnode(graph, graph.Edges.EndNodes{i,2});
         cur_clus_start_vec = clus_start_vec{i};
         cur_clus_end_vec   = clus_end_vec{i};
         
         [r_rep_s, c_rep_s]= size(cur_clus_start_vec); % check the cluster size at start idx 
         [r_rep_e, c_rep_e]= size(cur_clus_end_vec); % check the cluster size at end_idx 
         
         
         
         for j = 1:c_rep_s
             
             for k = 1:c_rep_e
                 cur_clus_end_vecKth =  cur_clus_end_vec(k);% chosen cluster for end node
                 cur_clus_start_vecJth = cur_clus_start_vec(j);
                 
                 if(cur_clus_end_vecKth ~= cur_clus_start_vecJth) % no nodes between same cluster
                     str_s =   sprintf('%s-%d', node_nam_start_vec{i}, cur_clus_start_vecJth); % string without concatenation
                     str_e =   sprintf('%s-%d', node_nam_end_vec{i}, cur_clus_end_vecKth);
                     
                     %IN_tab_Interedge = [IN_tab_Interedge;table({str_s str_e}, [ graph.Edges.Weight(i)], 'VariableNames', {'EndNodes','Weight'})];
                     IN_Interedge_s(edge_count) = str_s;
                     IN_Interedge_e(edge_count) = str_e;
                     IN_Interedge_w(edge_count) = graph.Edges.Weight(i);
%                      IN_tab_Interedge(edge_count,:) = table({str_s str_e}, [graph.Edges.Weight(i)], 'VariableNames', {'EndNodes','Weight'});
                     edge_count = edge_count+1;
                end
             end
             
             
         end
         
         
         
         
    end
    
    

end

%undir_edgetab = table(cell(r_Gnodes*2,2), zeros(r_Gnodes*2,1), 'VariableNames', {'EndNodes','Weight'}); 