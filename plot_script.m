
%%
guards_x = guard_target_struct.guards_x;
guards_y = guard_target_struct.guards_y;

targets_x = guard_target_struct.targets_x;
targets_y = guard_target_struct.targets_y;

for i = 1:length(guard_target_struct(1).guards_x)

    node_nam{i} = sprintf('  V%d',i);

end

figure(1);

%        plot3(guard_target_struct(1).guards_x, guard_target_struct(1).guards_y, 0.4*ones(1,length(guard_target_struct(1).guards_y)), 'o','Markersize',5,'MarkerEdgeColor','k','MarkerFaceColor','k');
%        text(guard_target_struct(1).guards_x,guard_target_struct(1).guards_y, 0.4*ones(1,length(guard_target_struct(1).guards_y)), node_nam,'HorizontalAlignment','left','FontSize',8);


for i = 1:length(guard_target_struct(1).targets_x)

    target_nam{i} = sprintf(' T%d',i);

end


%%
%plot visilibility
    environment_min_x = min(environment{1}(:,1));
    environment_max_x = max(environment{1}(:,1));
    environment_min_y = min(environment{1}(:,2));
    environment_max_y = max(environment{1}(:,2));
    X_MIN = environment_min_x-0.1*(environment_max_x-environment_min_x);
    X_MAX = environment_max_x+0.1*(environment_max_x-environment_min_x);
    Y_MIN = environment_min_y-0.1*(environment_max_y-environment_min_y);
    Y_MAX = environment_max_y+0.1*(environment_max_y-environment_min_y);
    clf; set(gcf,'position',[200 500 700 600]); hold on;
    axis equal; axis on; axis([X_MIN X_MAX Y_MIN Y_MAX]);

    patch( environment{1}(:,1) , environment{1}(:,2) , 0.1*ones(1,length(environment{1}(:,1)) ) , ...
           'w' , 'linewidth' , 1.5 );
    for i = 2 : size(environment,2)
        patch( environment{i}(:,1) , environment{i}(:,2) , 0.1*ones(1,length(environment{i}(:,1)) ) , ...
               'k' , 'EdgeColor' , [0 0 0] , 'FaceColor' , [0.8 0.8 0.8] , 'linewidth' , 1.5 );
    end



    p1 = plot3(guards_x, guards_y, 0.4*ones(1,length(guards_y)), 'o','Markersize',5,'MarkerEdgeColor','k','MarkerFaceColor','k');
    text(guards_x,guards_y, 0.4*ones(1,length(guards_y)), node_nam,'HorizontalAlignment','left','FontSize',8);


    p2 = plot3(targets_x, targets_y, 0.4*ones(1,length(targets_x)), 's','Markersize',5,'MarkerEdgeColor','r','MarkerFaceColor','r');
    text(targets_x,targets_y, 0.4*ones(1,length(targets_x)), target_nam,'HorizontalAlignment','left','FontSize',13);
    
    p3 = plot3(bot_x, bot_y, 0.4*ones(1,length(bot_x)), 's','Markersize',5,'MarkerEdgeColor','b','MarkerFaceColor','b');
    text(bot_x,bot_y, 0.4*ones(1,length(bot_x)), bot_nam,'HorizontalAlignment','left','FontSize',13);

 %%
 
 %plot with path
    figure();
    environment_min_x = min(environment{1}(:,1));
    environment_max_x = max(environment{1}(:,1));
    environment_min_y = min(environment{1}(:,2));
    environment_max_y = max(environment{1}(:,2));
    X_MIN = environment_min_x-0.1*(environment_max_x-environment_min_x);
    X_MAX = environment_max_x+0.1*(environment_max_x-environment_min_x);
    Y_MIN = environment_min_y-0.1*(environment_max_y-environment_min_y);
    Y_MAX = environment_max_y+0.1*(environment_max_y-environment_min_y);
    clf; set(gcf,'position',[200 500 700 600]); hold on;
    axis equal; axis on; axis([X_MIN X_MAX Y_MIN Y_MAX]);

    patch( environment{1}(:,1) , environment{1}(:,2) , 0.1*ones(1,length(environment{1}(:,1)) ) , ...
           'w' , 'linewidth' , 1.5 );
    for i = 2 : size(environment,2)
        patch( environment{i}(:,1) , environment{i}(:,2) , 0.1*ones(1,length(environment{i}(:,1)) ) , ...
               'k' , 'EdgeColor' , [0 0 0] , 'FaceColor' , [0.8 0.8 0.8] , 'linewidth' , 1.5 );
    end



    p1 = plot3(guards_x, guards_y, 0.4*ones(1,length(guards_y)), 'o','Markersize',5,'MarkerEdgeColor','k','MarkerFaceColor','k');
    text(guards_x,guards_y, 0.4*ones(1,length(guards_y)), node_nam,'HorizontalAlignment','left','FontSize',8);


    p2 = plot3(targets_x, targets_y, 0.4*ones(1,length(targets_x)), 's','Markersize',5,'MarkerEdgeColor','r','MarkerFaceColor','r');
    text(targets_x,targets_y, 0.4*ones(1,length(targets_x)), target_nam,'HorizontalAlignment','left','FontSize',13);
    
    p3 = plot3(bot_x, bot_y, 0.4*ones(1,length(bot_x)), 's','Markersize',5,'MarkerEdgeColor','b','MarkerFaceColor','b');
    text(bot_x,bot_y, 0.4*ones(1,length(bot_x)), bot_nam,'HorizontalAlignment','left','FontSize',13);
    
    X_path = cell(2,1);
    Y_path = cell(2,1);
    for i = 1:num_bots
        for j = 1:length(whole_path_nodes{i,1})
           if(j==1)
               cur_x = bot_x(str2num(whole_path_nodes{i,1}{j}(2:end)));
               X_path{i} = [X_path{i} cur_x]; 
               
               cur_y = bot_y(str2num(whole_path_nodes{i,1}{j}(2:end)));
               Y_path{i} = [Y_path{i} cur_x]; 
           elseif(j>1 && j< length(node_num_forw_cell{i}))
               cur_x = guards_x(str2num(whole_path_nodes{i,1}{j}(2:end)));
               X_path{i} = [X_path{i} cur_x];
               
               cur_y = guards_y(str2num(whole_path_nodes{i,1}{j}(2:end)));
               Y_path{i} = [Y_path{i} cur_x]; 
           elseif (j==length(node_num_forw_cell{i}))
                
               cur_x = bot_x(str2num(whole_path_nodes{i,1}{j}(2:end)));
               X_path{i} = [X_path{i} cur_x]; 
               
               cur_y = bot_y(str2num(whole_path_nodes{i,1}{j}(2:end)));
               Y_path{i} = [Y_path{i} cur_x];
           end
        end
        
    end
    
    
    p4 = plot3(X_path{1},Y_path{1},  0.4*ones(1,length(X_path{1})), '-r');
    
    
%%
    %break just after G_comp_temp = digraph([], []); line 74
    pdata_right = [3.98084484112396,7.95533836363301;9.51039443389633,0.414757900564329;3.15039773815711,0.414757900564329;3.15722177054659,4.00932637212185;5.94136112806105,3.20478998863015;0.411811535288686,6.76618145828282; 5.94136112806098,0.414757900564329; 0.411811535288686,3.17981057099319]; 
    IN_tab_Zeroedge_forplot = IN_transform_Zeroedge(G_comp);

    G_comp_temp_forplot = addedge(G_comp_temp, IN_tab_Zeroedge_forplot);
%     G_comp_temp1 = G_comp_temp; %addnode(G_comp_temp, IN_tab_Nodes); 
    figure;
    plot(G_comp_temp_forplot, 'XData',  pdata_right(:,1), 'YData',  pdata_right(:,2));
    title('I-N transformation ducplicate nodes');
    
    %%
    %solve gtsp line breakpoint before outsol and use following code
        [Out_sol, time_concorde_struct] = TSP_tour_Dat(G_gadget2,'/home/ashishkb/softwares/concorde/concorde/TSP/concorde');
    
    Out_solName = G_gadget2.Nodes.Name(Out_sol);
    countir = 1;
    for i = 1:3:length(Out_solName)
         Outsol_gadget{countir} = Out_solName{i}(3:end);
         countir  = countir + 1;
    end
    Outsol_gadget_num = findnode(G_gadget, Outsol_gadget);
    s_solg = zeros(G_gadget.numnodes,1);
    t_solg = zeros(G_gadget.numnodes,1);

    for i = 2:G_gadget.numnodes
        s_solg(i-1) = Outsol_gadget_num(i-1);
        t_solg(i-1) = Outsol_gadget_num(i);
    end

    s_solg(G_gadget.numnodes) = t_solg(end-1);
    t_solg(G_gadget.numnodes) = s_solg(1);
    idxOut = ~(findedge(G_gadget_temp,s_solg,t_solg)~=0);
    G_gadget_temp1 = addedge(G_gadget_temp,s_solg(idxOut),t_solg(idxOut), zeros(length(s_solg(idxOut)),1));
    P_gadget2 = plot(G_gadget_temp1, 'XData',  P_gadge.XData, 'YData',  P_gadge.YData);
    highlight(P_gadget2,s_solg,t_solg, 'NodeColor','g','EdgeColor','r');

%%

format long;


%Robustness constant
epsilon = 0.000000001;


%Read environment and guards geometry from files
environment = read_vertices_from_file('./source_code/example2.environment');



%Calculate a good plot window (bounding box) based on outer polygon of environment
environment_min_x = min(environment{1}(:,1));
environment_max_x = max(environment{1}(:,1));
environment_min_y = min(environment{1}(:,2));
environment_max_y = max(environment{1}(:,2));
X_MIN = environment_min_x-0.1*(environment_max_x-environment_min_x);
X_MAX = environment_max_x+0.1*(environment_max_x-environment_min_x);
Y_MIN = environment_min_y-0.1*(environment_max_y-environment_min_y);
Y_MAX = environment_max_y+0.1*(environment_max_y-environment_min_y);


%Clear plot and form window with desired properties
clf; set(gcf,'position',[200 500 700 600]); hold on;
axis equal; axis on; axis([X_MIN X_MAX Y_MIN Y_MAX]);


%Plot Environment
patch( environment{1}(:,1) , environment{1}(:,2) , 0.1*ones(1,length(environment{1}(:,1)) ) , ...
       'w' , 'linewidth' , 1.5 );
for i = 2 : size(environment,2)
    patch( environment{i}(:,1) , environment{i}(:,2) , 0.1*ones(1,length(environment{i}(:,1)) ) , ...
           'k' , 'EdgeColor' , [0 0 0] , 'FaceColor' , [0.8 0.8 0.8] , 'linewidth' , 1.5 );
end


%%

 %%       
            plot3(guard_target_struct(1).targets_x, guard_target_struct(1).targets_y, 0.4*ones(1,length(guard_target_struct(1).targets_x)),...
            's','Markersize',5,'MarkerEdgeColor','r','MarkerFaceColor','r');
            text(guard_target_struct(1).targets_x,guard_target_struct(1).targets_y, 0.4*ones(1,length(guard_target_struct(1).targets_x)),...
            target_nam,'HorizontalAlignment','left','FontSize',13);
        
      
%%
     figure(15);
     LWidths2 = 5*(G_init.Edges.Weight+1)/max(G_init.Edges.Weight);
     p_whole = plot(G_init,'EdgeLabel', G_init.Edges.Weight, 'LineWidth',LWidths2);        
        
        
%%
        
        
      
 %       tic 
%for i= 1:10
    
node_vec = cellfun(@(x) strsplit(x(2:(regexp(x,'-','start')-1)),'-') , G_comp.Edges.EndNodes,'uni',0);
%cell2mat(cellfun(@(x) str2double(strsplit(x((regexp(x,'-','start')+1):end),'-')) , G_comp.Edges.EndNodes(:,1),'uni',0));
%end
%toc

%%
