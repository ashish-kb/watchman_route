function [fin_sol, fin_rm_redunt, G_init, G_gadget2, nodes_totsp, total_cost, whole_path, time_auto_for_struct] = auto_for_content(visibility_adjacency_matrix, guard_target_struct, counter_struct)
    
%     filepath_guards = 'ginput2.guards'; % where u want to store that data
%     filepath_nodetarget = 'ginput2.nodetarget'; % just number of nodes and targets
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
    %script_sol take the generated data and solve the problem
    
    time_auto_for_struct(1).total_time = tic;
    
    targets_x = guard_target_struct.targets_x; % 
    targets_y = guard_target_struct.targets_y; % 
    
    
    guards_x = guard_target_struct.guards_x;
    guards_y = guard_target_struct.guards_y;
    target_mat = zeros(length(targets_x), length(guards_x));

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

        guard_mat_weight(r_vis(i), c_vis(i)) = pdist2( [guards_x(r_vis(i)), guards_y(r_vis(i))], [guards_x(c_vis(i)), guards_y(c_vis(i))]);    

    end


    % Concorde rounds the distances to integers
    % So if points have distances less than 0.5 units, they will become 0
    % To avoid this, I am going to scale all the distances
    % Right now, I am just going to scale by 100
    % TODO a better way of scaling, do this only if necessary
    scaled_guard_mat_weight = round(10*guard_mat_weight);


    V_adj  = double(scaled_guard_mat_weight);

    %%

    [fin_sol, fin_rm_redunt, ~, ~, G_init, G_gadget, G_gadget2, time_concorde_struct] = solve_GTSP(V_adj, V_Cluster);
    nodes_totsp = G_gadget2.numnodes;

    time_auto_for_struct(1).concorde_time = time_concorde_struct(1).concorde_time;
    %%
    %script auto

    temp_fin_sol = fin_sol;
    cyclic = 1; % assumed not cyclic no split 
    cur_cycle_checked = 1;
    check_cyc_cnt = 2;
    last_node = 1;
    

    while (cyclic == 1)
    
    % repeat_check = find(~cellfun(@isempty,strfind(temp_fin_sol, temp_fin_sol{1}))); %test if there is a repeat
          conc_nod_first = strtok(strtok(temp_fin_sol{1}, '-'),'V');  %connected nod first 
          conc_nod_last = strtok(strtok(temp_fin_sol{end}, '-'),'V');  % connected last last

          conc_clus_first = strsplit(temp_fin_sol{1}, '-');  %connected cluster first 
          conc_clus_last = strsplit(temp_fin_sol{end}, '-');  % connected cluster last
          if((isequal(conc_clus_first{2},conc_clus_last{2}))|| (isequal(conc_nod_first,conc_nod_last)))
                            
             
             
             if ((cur_cycle_checked > length(temp_fin_sol)) && (isequal(conc_nod_first,conc_nod_last)) && (~isequal(conc_clus_first{2},conc_clus_last{2})))
                cyclic = 0;  % only one node so already cyclic and we have pushed the cycle more than once
             else
                 temp_fin_sol = circshift(temp_fin_sol,1);%cycle me 
                 cyclic = 1;
             end
          else
               cyclic = 0; % it is cyclic already
          end
          cur_cycle_checked = cur_cycle_checked+1;
          
          
          
     end


    fin_sol = temp_fin_sol;


    fin_rm_redunt(:) = {''};
    finsol_count = 1;
    solnode = cell(length(fin_sol),1);
    solclus = cell(length(fin_sol),1);

    for i = 1:length(fin_sol)

        if(i==1)
            fin_rm_redunt{finsol_count}  = strtok(fin_sol{i}, '-');
            cur_node = strtok(fin_rm_redunt{finsol_count},'V');
            cur_clus = strsplit(fin_sol{i}, '-');
            solclus{finsol_count} = cur_clus{2};
            solnode{finsol_count} = cur_node;
            finsol_count = finsol_count+1;
        else
             cur_node  = strtok(strtok(fin_sol{i}, '-'),'V');
             cur_clus =  fin_sol{i}((regexp(fin_sol{i},'-','start')+1):end)   ;%strsplit(fin_sol{i}, '-');
             solclus{i} = cur_clus;       
             solnode{i} = cur_node;
             if(~(isequal(cur_clus,solclus{i-1}) || isequal(cur_node,solnode{i-1})))
                 fin_rm_redunt{finsol_count}  = strtok(fin_sol{i}, '-');             
                 finsol_count = finsol_count+1;
             end

        end
    end
    %%
    % whole_path = cell(G_init.numnodes, 1);% = length(cur_path)
    % counter_whlpath = 1; % goes through whole path
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
    
    if( length(fin_rm_redunt) ==1) % won't go into previous loop as length is 1
       whole_path = fin_rm_redunt(1);
    end

        total_cost = total_cost + distances(G_init, findnode(G_init,whole_path{end}), findnode(G_init,whole_path{1}));
        end_path = shortestpath(G_init, whole_path{end}, whole_path{1}); % path between last node and first node
        whole_path = [whole_path; end_path(2:end)'];
%[sdf';sas(2:end)']



%%
  % following we are doing to highlight the solution
%     s_whole = zeros((length(whole_path)-1),1);
%     t_whole = zeros((length(whole_path)-1),1);
% 
%     for i = 2:(length(whole_path))
%         s_whole(i-1) = findnode(G_init,whole_path{i-1});%Out_sol(i-1);
%         t_whole(i-1) = findnode(G_init,whole_path{i});%Out_sol(i);
%     end

%     s_whole(G_init.numnodes) = t_whole(end-1);
%     t_whole(G_init.numnodes) = s_whole(1);
%%
    %plot the tour in the graph
%     figure(counter_struct+200);
%     LWidths2 = 5*(G_init.Edges.Weight+1)/max(G_init.Edges.Weight);
%     p_whole = plot(G_init,'EdgeLabel', G_init.Edges.Weight, 'LineWidth',LWidths2);
%     p_whole.NodeColor = 'r';
%     highlight(p_whole,s_whole,t_whole, 'NodeColor','g','EdgeColor','r'); % if given edges solution is wrong this will not highlight the solution
%     
%     
%     figure(counter_struct);
%     % plot3( [guards(i,1) guards(j,1)], [guards(i,2) guards(j,2)], 0.3*[1 1], 'k', 'LineWidth', 0.5 , 'LineStyle' , '-' );
% 
%     for i = 1:length(s_whole)
%         plot3( [guards(s_whole(i),1) guards(t_whole(i),1)], [guards(s_whole(i),2) guards(t_whole(i),2)], 0.3*[1 1], 'b', 'LineWidth', 3 , 'LineStyle' , '-' );
%         
%     end
    
    
     time_auto_for_struct(1).total_time = toc(time_auto_for_struct(1).total_time);
    
end