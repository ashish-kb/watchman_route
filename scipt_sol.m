% nm_nodes = 5; % coming from Visilibility
% nm_targets = 6;
%run script_complete.m after this if not there then 
target_mat = zeros(length(targets_x), length(guards_x));

for i = 1:length(targets_x)

   target_mat(i,:)  = (visibility_adjacency_matrix((length(guards_x)+i) , (1:length(guards_x)))==1)*i;
   

end

V_Cluster = cell(length(guards_x), 1);

for i = 1:length(guards_x)
    V_Cluster{i} = target_mat(find(target_mat(:,i)), i)';
    
end

% target_mat(find(target_mat(:,12)), 12)
%b - diag(diag(b))
% a = rand(5);
% b = triu(a) + triu(a,1)'


%%
% Eucledean Distance
guard_mat = visibility_adjacency_matrix(1:length(guards_x),1:length(guards_x));

guard_mat = guard_mat - diag(diag(guard_mat));

guard_mat_weight = zeros(size(guard_mat));

[r_vis c_vis]= find(guard_mat ~=0); % row column of visiblility 

for i = 1:length(r_vis)
    
    guard_mat_weight(r_vis(i), c_vis(i)) = pdist2( [guards_x(r_vis(i)), guards_y(r_vis(i))], [guards_x(c_vis(i)), guards_y(c_vis(i))]);    
    
end


% Concorde rounds the distances to integers
% So if points have distances less than 0.5 units, they will become 0
% To avoid this, I am going to scale all the distances
% Right now, I am just going to scale by 100
% TODO a better way of scaling, do this only if necessary
scaled_guard_mat_weight = round(100*guard_mat_weight);


V_adj  = double(scaled_guard_mat_weight);

%%
% random assignment

% rand_mat = 10*randi([1 10], nm_nodes, 'uint8');
% 
% rand_mat(guard_mat~=1) = 0; % replace with zeros where there are no edges
% 
% guard_mat_weight = rand_mat;
% 
% guard_mat_weight = guard_mat_weight - diag(diag(guard_mat_weight));
% guard_mat_temp = tril(guard_mat_weight); % making it symmtric as we have symmetric edges
% guard_mat_weight  = tril(guard_mat_temp,-1)'+ guard_mat_temp;
% V_adj  = double(guard_mat_weight);

%%


[fin_sol, fin_rm_redunt, Out_solName, Out_sol, G_init, G_gadget, G_gadget2] = solve_GTSP(V_adj, V_Cluster);
% 
% 
% 
% givn_nodes = length(V_adj);
% nod_name  = cell(1, givn_nodes);
% 
% for i = 1:given_nodes
% 
%     nod_name{i} = sprintf('V%d',i);
% 
% end
% 
%  G_init_truew = digraph(V_adj, node_name); % initial graph with true weights
    

