function [outfin_sol, outfin_cost, Out_solName, whole_path_nodes, G_init, G_nodebot, edges_totsp, nodes_totsp, time_concorde_struct] = auto_for_content_Multi_noon_bean(visibility_adjacency_matrix, guard_target_struct, counter_struct)
    %[outfin_sol, outfin_cost,Out_solName, Out_sol, G_init, edges_totsp, nodes_totsp, time_concorde_struct] = solve_GTSP_noon_bean(V_adj, V_Cluster);
    %filepath_guards = 'ginput2.guards'; % where u want to store that data
    %filepath_nodetarget = 'ginput2.nodetarget'; % just number of nodes and targets
    time_auto_for_struct = struct('total_time',{},'concorde_time',{});
    
    
    



    %%
    % write data to files


%     fid = fopen(filepath_guards,'W');
%     if fid < 0
%         error('Cannot create  file');
%         return;
%     end
% 
%     fprintf(fid,'//Guard x-y coordinates \n');
% 
%     for i = 1:length(guards)
%         fprintf(fid,'%d %d\n', guards(i,1), guards(i,2));
%     end
% 
%     fclose(fid);
% 
% 
%     fid = fopen(filepath_nodetarget,'W');
%     if fid < 0
%         error('Cannot create file');
%         return;
%     end
%     fprintf(fid, '%d  %d', length(guards_x), length(targets_x));
% 
%     fclose(fid);


    %%
    %
    
    time_auto_for_struct(1).total_time = tic;
    
    targets_x = guard_target_struct.targets_x; % 
    targets_y = guard_target_struct.targets_y; % 
    
    
    guards_x = guard_target_struct.guards_x;
    guards_y = guard_target_struct.guards_y;
    
    guards = guard_target_struct.guards;
    bot_x = guard_target_struct.bot_x;
    bot_y = guard_target_struct.bot_y;
    target_mat = zeros(length(targets_x), length(guards_x));
    %target_mat is 
    for i = 1:length(targets_x)

       target_mat(i,:)  = (visibility_adjacency_matrix((length(guards_x)+i) , (1:length(guards_x)))==1)*i;


    end

    V_Cluster = cell(length(guards_x), 1);

    for i = 1:length(guards_x)
        V_Cluster{i} = target_mat(find(target_mat(:,i)), i)';

    end

    %%
    % Eucledean Distance
    guard_mat = visibility_adjacency_matrix(1:length(guards_x),1:length(guards_x));

    guard_mat = guard_mat - diag(diag(guard_mat));

    guard_mat_weight = zeros(size(guard_mat));

    [r_vis c_vis]= find(guard_mat ~=0); % row column of visiblility 

    for i = 1:length(r_vis)

        guard_mat_weight(r_vis(i), c_vis(i)) = pdist2([guards_x(r_vis(i)), guards_y(r_vis(i))], [guards_x(c_vis(i)), guards_y(c_vis(i))]);    

    end


    % Concorde rounds the distances to integers
    % So if points have distances less than 0.5 units, they will become 0
    % To avoid this, I am going to scale all the distances
    % Right now, I am just going to scale by 100
    % TODO a better way of scaling, do this only if necessary
    scaled_guard_mat_weight = guard_mat_weight;%round(10*guard_mat_weight);


    V_adj  = double(scaled_guard_mat_weight);
    
    %%
    %new V_adj for bots
    
    %getbot_col = visibility_adjacency_matrix(1:length(guards_x), (length(visibility_adjacency_matrix) - length(bot_x)+1):end); %
    %getbot_row = getbot_col';
    %cat_interbot = visibility_adjacency_matrix((length(visibility_adjacency_matrix) - length(bot_x)+1):end, (length(visibility_adjacency_matrix) - length(bot_x)+1):end);
    %bot_mat = vertcat(horzcat(visibility_adjacency_matrix(1:length(guards_x),1:length(guards_x)), getbot_col), [getbot_row cat_interbot]);
    
    visility_adj_guardbot = zeros(size(visibility_adjacency_matrix));
    bot_locat_invisilibility = (length(visibility_adjacency_matrix) - length(bot_x)+1):length(visibility_adjacency_matrix);
    visility_adj_guardbot(1:length(guards_x), 1:length(guards_x)) = visibility_adjacency_matrix(1:length(guards_x),1:length(guards_x));
    visility_adj_guardbot(bot_locat_invisilibility, bot_locat_invisilibility) = visibility_adjacency_matrix(bot_locat_invisilibility, bot_locat_invisilibility);
    visility_adj_guardbot(1:length(guards_x), bot_locat_invisilibility) = visibility_adjacency_matrix(1:length(guards_x), bot_locat_invisilibility);
    visility_adj_guardbot(bot_locat_invisilibility, 1:length(guards_x)) = visibility_adjacency_matrix(bot_locat_invisilibility, 1:length(guards_x));
    
    visility_adj_guardbot = visility_adj_guardbot - diag(diag(visility_adj_guardbot));
    visility_adj_guardbot_weight = zeros(size(guard_mat));
        
    [r_visbot c_visbot]= find(visility_adj_guardbot ~=0);
    
     for i = 1:length(r_visbot)

        visility_adj_guardbot_weight(r_visbot(i), c_visbot(i)) = pdist2([guards(r_visbot(i), 1), guards(r_visbot(i), 2)], [guards(c_visbot(i), 1), guards(c_visbot(i), 2)]);    

     end
    
    V_adj_bot = double(visility_adj_guardbot_weight);
    
    %removing target rows and cols 
    V_adj_bot = V_adj_bot([1:length(guards_x) (length(guards_x)+length(targets_x)+1):length(V_adj_bot)], [1:length(guards_x) (length(guards_x)+length(targets_x)+1):length(V_adj_bot)]);

    
    % I think putting zeros where we have targets and then extracting
    % guard+bot_mat_weight is better... after that remove the zeros(where targets are) while making the
    % graph in solveMulti_GTSP_noon_bean, then use distances and then use
    % these distance to make edges to all nodes finally
    
    
    

    %%

    [ outfin_sol, outfin_cost, Out_solName, whole_path_nodes,G_init, G_nodebot, edges_totsp, nodes_totsp, time_concorde_struct] = solveMulti_GTSP_noon_bean(V_adj, V_Cluster, V_adj_bot);
    

    time_auto_for_struct(1).concorde_time = time_concorde_struct(1).concorde_time;
     
    time_auto_for_struct(1).total_time = toc(time_auto_for_struct(1).total_time);
    
end