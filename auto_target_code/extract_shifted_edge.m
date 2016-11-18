edge_table = G_atsp.Edges;
zero_weight_ind = edge_table.Weight==0;

cur_clus_cell = cell(1, length(Cluster_to_node));
weight_vec_cell = cell(num_bots, 1);
str_f_cell = cell(num_bots,1);

for i = 1: length(Cluster_to_node)
    for j = 1:num_bots
        log_ind_clus(:,i) = cell2mat(cellfun(@(x) isequal(sprintf('-%d', i),x(regexp(x,'-','start'):end)), edge_table.EndNodes(:,2), 'uni', 0 ));
        cur_clus_ind = log_ind_clus(:,i) & zero_weight_ind;
        cur_clus_cell{1,i} = edge_table.EndNodes(cur_clus_ind==1, 1); % cell containing names of all the nodes in cluster i connected in a cycle
        
        weight_vec_cell{j,1} = ...
        circshift(G_nodebot_comp.Edges.Weight(findedge(G_nodebot_comp, cellfun(@(x) sprintf('B%d', j), cur_clus_cell{1,i},'uni', 0), cellfun(@(x) x(1:(regexp(x,'-','start')-1)) , cur_clus_cell{1,i},'uni',0))),-1); 
        
        str_f_cell{j, 1} = cellfun(@(x) sprintf('B%d-f', j), cur_clus_cell{1,i},'uni', 0);
        
        
    end
    
    concat_f_cell = [str_f_cell{:,:}];
    concat_w_cell = [weight_vec_cell{:,:}];
    G_atsp = addedge(G_atsp, repmat(cur_clus_cell{1,i}, num_bots, 1), concat_f_cell(:), [concat_w_cell(:)]); 
    
end


%G_atsp = addedge(G_atsp, repmat(cur_clus_cell{1,i}, 2, 1), concat_f_cell(:), [concat_w_cell(:)]); 
%x{1}(regexp(x{1},'-','start')+1:end)

%saasd = log_ind_clus(:,i) & zero_weight_ind;
%edge_table.EndNodes(find(saasd), :)


% sada = cell2mat(cellfun(@(x) isequal(sprintf('-f'),x(regexp(x,'-','start'):end)), edge_table.EndNodes(:,2), 'uni', 0 ));
% adsss = edge_table.EndNodes(find(sada), :);