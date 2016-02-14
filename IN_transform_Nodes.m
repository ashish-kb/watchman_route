function IN_tab_Nodes = IN_transform_Nodes(graph) % the graph should contain Graph.Nodes.Cluster values....outputs a table
   
    [sz_r_clst, sz_c_clust]= size(graph.Nodes.Cluster);
    node_count = 1;
    
    num_new = sum(cellfun(@numel, graph.Nodes.Cluster)); % number of new nodes
    
    IN_tab_Nodes =  table(cell(num_new,1), cell(num_new,1), 'VariableNames', {'Name','Cluster'}); % external edge tab
    
    
    for i = 1:sz_r_clst
        
        [r_rep, c_rep]= size(graph.Nodes.Cluster{i}); % checking if repeating nodes are there
        for j= 1:c_rep

                str1 = graph.Nodes.Name{i};  % choosing the previous node 'V1'
                str  = sprintf('%s-%d', str1, graph.Nodes.Cluster{i}(j));  % concatenate with current cluster
                cur_Name = str;
                cur_Cluster = graph.Nodes.Cluster{i}(j);
                IN_tab_Nodes(node_count,:) = table({cur_Name}', {cur_Cluster}', 'VariableNames', {'Name' 'Cluster'});
                node_count = node_count + 1;
          

        end
        %G_comp.Nodes.Name{2} = 'V20';
        
        
        
        
    end
    
    
    
    

end


%addabc_table = table({str_a str_b str_c}', {cur_clus cur_clus cur_clus}', 'VariableNames', {'Name' 'Cluster'});