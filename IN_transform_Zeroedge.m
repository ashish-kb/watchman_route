function IN_tab_Zeroedge = IN_transform_Zeroedge(graph) % the graph should contain Graph.Nodes.Cluster values....outputs a table - its part of in transformation
    

    node_table = graph.Nodes;
    [sz_r_clst, sz_c_clust]= size(graph.Nodes.Cluster);
    num_edges =0;
    for i = 1:length(graph.Nodes.Cluster)
        if(length(graph.Nodes.Cluster{i})>1)
            num_edges = num_edges + length(nchoosek(graph.Nodes.Cluster{i},2)); % edge will only be there if 2 nodes 
        end
    end
    
    %num_edges = sum(factorial(cellfun(@numel, graph.Nodes.Cluster))); % number of zero edges are equal to sum of factorial of clusters in it
    IN_tab_Zeroedge = table(cell(num_edges,2), zeros(num_edges,1), 'VariableNames', {'EndNodes','Weight'}); 
    
    edge_count = 1;
   
    for i = 1:sz_r_clst

       [r_rep, c_rep]= size(graph.Nodes.Cluster{i}); % checking if repeating nodes are there
        for j = 1:c_rep
            for k = 1:c_rep
                if (k~=j)
                   %G_comp = addedge(G_comp, copy_node(j), copy_node(k), 0); % adding zero cost edges 
                   % ***step will take a lot of time  can be altered               
                   str1 = graph.Nodes.Name{i};  % choosing the previous node 'V1'
                   str_n1  = sprintf('%s-%d', str1, graph.Nodes.Cluster{i}(j));  % concatenate with current cluster - node 1 of the edge we create
                   str_n2  = sprintf('%s-%d', str1, graph.Nodes.Cluster{i}(k));  % concatenate with current cluster - node 2 of the edge we create
                   
                   IN_tab_Zeroedge(edge_count, :) = table({str_n1 str_n2}, [0]',  'VariableNames', {'EndNodes','Weight'}); 
                   edge_count = edge_count + 1;
                end
            end
        end

    end
    
    
    
    

end

%undir_edgetab = table(cell(r_Gnodes*2,2), zeros(r_Gnodes*2,1), 'VariableNames', {'EndNodes','Weight'}); 