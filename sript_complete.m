% clear;
% clc;
% close all;
% 
% 
% V_adj =    [0     10    0     0;
%             10     0    15     0;
%             0     15     0     20;
%             0      0     20     0];
%         
%         
%  node_name =  {'V1','V2','V3','V4'};  
%  
%   
%  G_init = digraph(V_adj, node_name);
%  
%   figure;
%  plot(G_init,'EdgeLabel', G_init.Edges.Weight);
%  
%  %%
%  
%  

 temp_fin_sol = fin_sol;
 cyclic = 1; % assumed not cyclic no split 
 cur_cycle_checked = 1;
 check_cyc_cnt = 2;
 last_node = 1;
 
 while (cyclic ==1)
     
%       repeat_check = find(~cellfun(@isempty,strfind(temp_fin_sol, temp_fin_sol{1}))); %test if there is a repeat
      conc_nod_first = strtok(strtok(temp_fin_sol{1}, '-'),'V');  %connected nod first 
      conc_nod_last = strtok(strtok(temp_fin_sol{end}, '-'),'V');  % connected last last
      
      conc_clus_first = strsplit(temp_fin_sol{1}, '-');  %connected cluster first 
      conc_clus_last = strsplit(temp_fin_sol{end}, '-');  % connected cluster last
      
%       if(length(repeat_check)>1)
%           cur_cycle_checked = 1;
%           check_cyc_cnt = 2;
%           clus_cur1 =  strsplit(temp_fin_sol{1}, '-'); % cluster of node 1 - clus_cur1{2} actually gives the cluster number
%           
%           while(cur_cycle_checked == 1)
%               clus_rep_prev =  strsplit(temp_fin_sol{check_cyc_cnt}, '-'); % check if node previous to repeated is in the same cluster otherwise cycle the array
%              
%                             
%               if (check_cyc_cnt ==(repeat_check(end)))
%                       
%                   if(isequal(clus_cur1{2}, clus_rep_prev{2}) && ~isequal(repeat_check(1), (repeat_check(end)-1)))
%                       cur_cycle_checked = 0; % current cycle is checked
%                       cyclic = 0; % it has no split in cycle
%                   else
%                       cur_cycle_checked = 0;
%                       temp_fin_sol = circshift(temp_fin_sol,1);%cycle me 
%                       cyclic = 1;
%                   end
%                   last_node = 0;
%               end 
%               
%               if(~isequal(clus_cur1{2}, clus_rep_prev{2}) && (last_node ==1))
%                   cur_cycle_checked = 0;             
%                   
%                   temp_fin_sol = circshift(temp_fin_sol,1);%cycle me 
%                   cyclic = 1;
%                          
%                         %cyclic = 0; % it is cyclic now.. cycle spit is over
%               end
%               
%               check_cyc_cnt = check_cyc_cnt+1;
%           end
          
         
     if((isequal(conc_clus_first{2},conc_clus_last{2}))|| (isequal(conc_nod_first,conc_nod_last)))
            
         temp_fin_sol = circshift(temp_fin_sol,1);%cycle me 
         cyclic = 1;
      else
           cyclic = 0; % it is cyclic already
      end
 end
 
 %%
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
         cur_clus = strsplit(fin_sol{i}, '-');
         solclus{i} = cur_clus{2};       
         solnode{i} = cur_node;
         if(~(isequal(cur_clus{2},solclus{i-1}) || isequal(cur_node,solnode{i-1})))
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

    total_cost = total_cost + distances(G_init, findnode(G_init,whole_path{end}), findnode(G_init,whole_path{1}));
    end_path = shortestpath(G_init, whole_path{end}, whole_path{1}); % path between last node and first node
    whole_path = [whole_path; end_path(2:end)'];
%[sdf';sas(2:end)']



%%
  % following we are doing to highlight the solution
    s_whole = zeros((length(whole_path)-1),1);
    t_whole = zeros((length(whole_path)-1),1);

    for i = 2:(length(whole_path))
        s_whole(i-1) = findnode(G_init,whole_path{i-1});%Out_sol(i-1);
        t_whole(i-1) = findnode(G_init,whole_path{i});%Out_sol(i);
    end

%     s_whole(G_init.numnodes) = t_whole(end-1);
%     t_whole(G_init.numnodes) = s_whole(1);
%%
    %plot the tour in the graph
    figure;
    LWidths2 = 5*(G_init.Edges.Weight+1)/max(G_init.Edges.Weight);
    p_whole = plot(G_init,'EdgeLabel', G_init.Edges.Weight, 'LineWidth',LWidths2);
    p_whole.NodeColor = 'r';
    highlight(p_whole,s_whole,t_whole, 'NodeColor','g','EdgeColor','r'); % if given edges solution is wrong this will not highlight the solution
 
%%

 figure(1);
% plot3( [guards(i,1) guards(j,1)], [guards(i,2) guards(j,2)], 0.3*[1 1], 'k', 'LineWidth', 0.5 , 'LineStyle' , '-' );

    for i = 1:length(s_whole)
        plot3( [guards(s_whole(i),1) guards(t_whole(i),1)], [guards(s_whole(i),2) guards(t_whole(i),2)], 0.3*[1 1], 'b', 'LineWidth', 3 , 'LineStyle' , '-' );
        
    end
    
    
    
    
    
    %%
    
%     %plotting
%         figure;
%     LWidths2 = 5*(G_gadget2.Edges.Weight+1)/max(G_gadget2.Edges.Weight);
%     p_gadget2 = plot(G_gadget2, 'Layout','Circle','EdgeLabel', G_gadget2.Edges.Weight, 'LineWidth',LWidths2);
%     p_gadget2.NodeColor = 'r';
  %  highlight(p_gadget2,s_sol,t_sol, 'NodeColor','g','EdgeColor','r'); % if given edges solution is wrong this will not highlight the solution
    %%
    % plot spanning tree
    G_gadget2_T = minspantree(G_gadget2);
     figure;
%    % LWidths2 = 5*(G_gadget2.Edges.Weight+1)/max(G_gadget2.Edges.Weight);
     p_gadget2_T = plot(G_gadget2_T, 'NodeLabel',G_gadget2_T.Nodes.Name,'Layout','force','EdgeLabel', G_gadget2_T.Edges.Weight); %, 'LineWidth',LWidths2);
     p_gadget2_T.NodeColor = 'r';
    %%
    
        %%
    % following we are doing to highlight the solution
    s_sol = zeros(G_gadget2.numnodes,1);
    t_sol = zeros(G_gadget2.numnodes,1);

    for i = 2:G_gadget2.numnodes
        s_sol(i-1) = Out_sol(i-1);
        t_sol(i-1) = Out_sol(i);
    end

    s_sol(G_gadget2.numnodes) = t_sol(end-1);
    t_sol(G_gadget2.numnodes) = s_sol(1);
    
    figure;
    LWidths2 = 5*(G_gadget2.Edges.Weight+1)/max(G_gadget2.Edges.Weight);
  %  p_gadget2_n = plot(G_gadget2, 'Layout','Circle','EdgeLabel', G_gadget2.Edges.Weight, 'LineWidth',LWidths2);
    p_gadget2_n  = plot(G_gadget2, 'XData', p_gadget2_T.XData, 'YData', p_gadget2_T.YData, 'EdgeLabel', G_gadget2.Edges.Weight, 'LineWidth',LWidths2);
    
     highlight(p_gadget2_n,s_sol,t_sol, 'NodeColor','g','EdgeColor','r');
