function IN_tab_Interedge = IN_transform_Interedge(graph) % the graph should contain Graph.Nodes.Cluster values....outputs a table - its part of in transformation
    
    tot_clus = cellfun(@numel, graph.Nodes.Cluster); % matrix containing total of each cluster

    int_clus_tot = 0; % number of edges in the new graph

    for i = 1:length(tot_clus) % what this is doing is for a,b,c = a(a+b+c-a) + b(a+b+c-b) + c(a+b+c-c) = total number of Intercluster edges after addition of new nodes

        int_clus_tot = int_clus_tot + tot_clus(i)*(sum(tot_clus)-tot_clus(i));
    end
  
    IN_tab_Interedge = table(cell(int_clus_tot,2), zeros(int_clus_tot,1), 'VariableNames', {'EndNodes','Weight'}); 
    edge_count = 1;
    %IN_tab_Interedge = table({}, [], 'VariableNames', {'EndNodes','Weight'});
    
    for i = 1:graph.numedges
         start_idx = findnode(graph, graph.Edges.EndNodes{i,1});
         end_idx   = findnode(graph, graph.Edges.EndNodes{i,2});
         
         [r_rep_s, c_rep_s]= size(graph.Nodes.Cluster{start_idx}); % check the cluster size at start idx 
         [r_rep_e, c_rep_e]= size(graph.Nodes.Cluster{end_idx}); % check the cluster size at end_idx 
         
         for j = 1:c_rep_s
             
             for k = 1:c_rep_e
                 str_s =   sprintf('%s-%d', graph.Nodes.Name{start_idx}, graph.Nodes.Cluster{start_idx}(j)); % string without concatenation
                 str_e =   sprintf('%s-%d', graph.Nodes.Name{end_idx}, graph.Nodes.Cluster{end_idx}(k));
                 
                 %IN_tab_Interedge = [IN_tab_Interedge;table({str_s str_e}, [ graph.Edges.Weight(i)], 'VariableNames', {'EndNodes','Weight'})];
                 IN_tab_Interedge(edge_count,:) = table({str_s str_e}, [ graph.Edges.Weight(i)], 'VariableNames', {'EndNodes','Weight'});
                 edge_count = edge_count+1;
             end
             
             
         end
         
         
         
         
    end
    
    

end

%undir_edgetab = table(cell(r_Gnodes*2,2), zeros(r_Gnodes*2,1), 'VariableNames', {'EndNodes','Weight'}); 