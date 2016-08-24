        

counter_struct =1;
        
        


        [fin_sol, fin_rm_redunt, G_init, G_gadget2, nodes_totsp, total_cost, whole_path, time_auto_for_struct] = auto_for_content(visibility_adjacency_matrix, guard_target_struct, counter_struct);
        
        %%
        
         [fin_sol, fin_rm_redunt, Out_solName, Out_sol, G_init, G_gadget, G_gadget2, time_concorde_struct] = solve_GTSP(V_adj, V_Cluster) ;
         
         %%
         %********noon-bean************
         counter_struct =1;
        
        %currently cost I
        %
        % <http://www.mathworks.com N gtsp_to_atsp is changed ..> .we are using max max
        % 
        % * ITEM1
        % * ITEM2
        % 
        %instead OF SUMSUM


        [outfin_sol, outfin_cost,Out_solName, Out_sol, G_init, edges_totsp, nodes_totsp, time_auto_for_struct] = auto_for_content_noon_bean(visibility_adjacency_matrix, guard_target_struct, counter_struct);
        %%
        outfin_sol_unique = [unique(outfin_sol,'stable') outfin_sol(1)];
        
        fid = fopen('goals.txt','W');
        if fid < 0
            error('Cannot create  file');
            return;
        end

        %fprintf(fid,'//goals or outfin_sol_unique x-y coordinates \n');
        %goals or outfin_sol_unique x-y coordinates
        for i = 1:length(outfin_sol_unique)
            fprintf(fid,'%f\t%f\n', guards_x(outfin_sol_unique(i)), guards_y(outfin_sol_unique(i)));
        end

        fclose(fid);
        
        %currently cost IN gtsp_to_atsp is changed ...we are using max max
        %instead OF SUMSUM
         
         %%
%          tic;
%          for i =1:100000
%              x_ss = ismember(X,[1 2 3; 2 3 4]);
%          end
%          toc;

%%

%Adj_G_comp = full(adjacency(G_comp));
[s_ t_]=findedge(G_comp);

 Adj_G_comp(:,:) = sum(G_comp.Edges.Weight(:));
 Adj_G_comp_ind = sub2ind(size(Adj_G_comp), s_(:),t_(:));
 
 Adj_G_comp(Adj_G_comp_ind(:)) = G_comp.Edges.Weight(:);
 diag_ind = sub2ind(size(Adj_G_comp), 1:length(Adj_G_comp),1:length(Adj_G_comp));
 Adj_G_comp(diag_ind(:)) = 0;
 
 setMap = []
 
 
 
 %%
% ***nodes20target10_40s***
% gtsp_ILP ->sol= 1,19,8,15,14,4,17,1 =>cluse => 1, 4, 8, 9, 3, 4, 2, 10, 5, 6, 7, 1  => all fine cost 161 approz
% normal_tsp-> sol = 1,18,19,7,14,13,12,15,10,9,8,7,6,5,3,17,3,1  => clus=>  % [1,8,9,4,8,9,4,5,5,10,2,10,3,3,3,4,4,6,6,7,6] total cost: 1616/10
%  {'V20-9','V19-4','V5-6','V18-8','V14-5','V12-10','V17-7','V16-2','V8-3','V2-1'}
%  
%  
%  

%%



