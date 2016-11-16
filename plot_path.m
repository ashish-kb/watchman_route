

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
               Y_path{i} = [Y_path{i} cur_y]; 
           elseif(j>1 && j< length(whole_path_nodes{i,1}))
               cur_x = guards_x(str2num(whole_path_nodes{i,1}{j}(2:end)));
               X_path{i} = [X_path{i} cur_x];
               
               cur_y = guards_y(str2num(whole_path_nodes{i,1}{j}(2:end)));
               Y_path{i} = [Y_path{i} cur_y]; 
           elseif (j==length(whole_path_nodes{i,1}))
                
               cur_x = bot_x(str2num(whole_path_nodes{i,1}{j}(2:end)));
               X_path{i} = [X_path{i} cur_x]; 
               
               cur_y = bot_y(str2num(whole_path_nodes{i,1}{j}(2:end)));
               Y_path{i} = [Y_path{i} cur_y];
           end
        end
        
    end
    
    hold on;
    
    
    p4 = plot3(X_path{1},Y_path{1},  0.4*ones(1,length(X_path{1})), '-r');
    
    p4 = plot3(X_path{2},Y_path{2},  0.4*ones(1,length(X_path{2})), '-b');
    