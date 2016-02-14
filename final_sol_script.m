fin_rm_redunt = fin_rm_redunt(~cellfun('isempty',fin_rm_redunt)); % removing empty cells
whole_path = {};
total_cost = 0;
for i = 1:(length(fin_rm_redunt)-1)
    
     cur_path = shortestpath(G_init, fin_rm_redunt{i}, fin_rm_redunt{i+1}); % path between currently selected nodes
     total_cost = total_cost + distances(G_init, findnode(G_init,fin_rm_redunt{i}), findnode(G_init,fin_rm_redunt{i+1}));
     if(i==1)
         whole_path = fin_rm_redunt(i);
         
     end
     whole_path =  [whole_path; cur_path(2:end)']; %adding cur_path to whole_path
    
end
    total_cost = total_cost + distances(G_init, findnode(G_init,fin_rm_redunt{length(fin_rm_redunt)}), findnode(G_init,fin_rm_redunt{1}));
    whole_path = [whole_path; fin_rm_redunt{1}];