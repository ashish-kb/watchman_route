%*****change filepath 

%Clear the desk
clear all; close all; clc;

addpath(genpath('/home/ashishkb/Dropbox/MATLAB/RAASLAB/gtsp/solved_cases/'));
addpath(genpath('/home/ashishkb/Dropbox/MATLAB/RAASLAB/gtsp/source_code/'));


format long;


%Robustness constant
epsilon = 0.000000001;


%Read environment and guards geometry from files
environment = read_vertices_from_file('./source_code/gazebo_rviz.environment');



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



% prompt = 'Input number of guards. Auto guards placed at environment corner  ';
% nm_guards = inputdlg(prompt);
% 
% prompt2 = 'Input number of targets';
% nm_targets = inputdlg(prompt2);
% field1 = 'fin_sol';  value1 = {};
% field2 = 'fin_rm_redunt';  value2 = {};
% field3 = 'G_init';  value3 = {};
% field4 = 'G_gadget';  value4 = {};
% 
% s_sol_struct = struct(field1,value1,field2,value2,field3,value3,field4,value4);

max_guards = 30; 
max_targets = 30; 
max_average = 1;
nodes_pushed_tsp = zeros(max_guards*max_targets, 1);
time_script_auto_concorde = zeros(max_guards, max_targets, max_average);%struct('total_time',[],'concorde_time',[]);
time_script_auto_total = zeros(max_guards, max_targets, max_average);

%%






environ_guards = cell2mat(environment(:));
%* guard_target_struct(1).environ_guardX = environ_guards(:,1);
%* guard_target_struct(1).environ_guardY = environ_guards(:,2);


 environ_guardX = environ_guards(:,1);
 environ_guardY = environ_guards(:,2);


%%
counter_struct = 1; % itself not a structure just counting through the structure
for rand_guard_required = 5:5:max_guards
    for rand_target_required = 5:5:max_targets
%         [x,y, key] = ginput(1);
%         if (key=='e')
%             display('Paused press some other key to continue');
%             pause;
%         
%         end     
%        msgbox()msgbox(sprintf('These Values are incorrect WARNING %f warn',k))
        
      %   uiwait(msgbox(sprintf('Target = %d, GuardsNew = %d',rand_target_required, rand_guard_required)));
        for average_generator = 1:max_average % this for loop ensures that average time for particular number of targets and guards is found
            guard_target_struct = struct('guards_x',{},'guards_y',{}, 'targets_x', {}, 'targets_y', {}, 'environ_guardX', {}, 'environ_guardY', {}, 'guards', {});
%             nodes_totsp = 0;
            
            guard_target_struct(1).environ_guardX = environ_guardX;
            guard_target_struct(1).environ_guardY = environ_guardY;
            
            got_guard_rands = 1;
            got_target_rands = 1;

            guard_randX = zeros(rand_guard_required,1);
            guard_randY = zeros(rand_guard_required,1);


            target_randX = zeros(rand_target_required,1);
            target_randY = zeros(rand_target_required,1);

            % generate random guards
            rng('shuffle');

            while (got_guard_rands <= rand_guard_required)

                X_try =    randi([ceil(environment_min_x)  floor(environment_max_x)],1 ,'double') + rand(1, 'double') ; % one gives int other gives fraction...some other way could be better
                Y_try =    randi([ceil(environment_min_y)  floor(environment_max_y)],1 ,'double') + rand(1, 'double') ;

                if(in_environment( [X_try Y_try] , environment , epsilon ))
                    guard_randX(got_guard_rands) = X_try;
                    guard_randY(got_guard_rands) = Y_try;

                    got_guard_rands  = got_guard_rands + 1;
                end

            end



            %plot3(guard_randX, guard_randY, 0.4*ones(1,length(guard_randY)), 'o','Markersize',5,'MarkerEdgeColor','k','MarkerFaceColor','k');


            % generate random targets

            while (got_target_rands <= rand_target_required)

                X_try =  randi([ceil(environment_min_x)  floor(environment_max_x)],1 ,'double') + rand(1, 'double') ; % one gives int other gives fraction...some other way could be better
                Y_try =  randi([ceil(environment_min_y)  floor(environment_max_y)],1 ,'double') + rand(1, 'double') ;

                if(in_environment([X_try Y_try] , environment , epsilon))
                    target_randX(got_target_rands) = X_try;
                    target_randY(got_target_rands) = Y_try;

                    got_target_rands  = got_target_rands + 1;
                end

            end


            %plot3(target_randX, target_randY, 0.4*ones(1,length(target_randY)), 's','Markersize',5,'MarkerEdgeColor','r','MarkerFaceColor','r');

            %guards = zeros((size(guard_randX,1) + size(environ_guardX,1) + size(target_randX,1)),2);

            guard_target_struct(1).guards = [[guard_randX, guard_randY]; [guard_target_struct(1).environ_guardX, guard_target_struct(1).environ_guardY];[target_randX, target_randY]];  % guards and targets concatenated
% 
% 
            guard_target_struct(1).targets_x = target_randX;
            guard_target_struct(1).targets_y = target_randY;
% % 
            guard_target_struct(1).guards_x = [guard_randX; guard_target_struct(1).environ_guardX];
            guard_target_struct(1).guards_y = [guard_randY; guard_target_struct(1).environ_guardY];
            
            %%


%*              visibility_adjacency_matrix = visibility_graph(guard_target_struct(1).guards, environment, epsilon);
             visibility_adjacency_matrix = visibility_graph(guard_target_struct(1).guards, environment, epsilon);




%            [fin_sol, fin_rm_redunt, G_init, G_gadget2, nodes_totsp, total_cost, whole_path, time_auto_for_struct] = auto_for_content_noon_bean(visibility_adjacency_matrix, guard_target_struct, counter_struct);
             [outfin_sol, outfin_cost,Out_solName, Out_sol, G_init, edges_totsp, nodes_totsp, time_auto_for_struct] = auto_for_content_noon_bean(visibility_adjacency_matrix, guard_target_struct, counter_struct);
             %[fin_sol, fin_rm_redunt, G_init, G_gadget2, nodes_totsp, total_cost, whole_path, time_auto_for_struct]
             filename = [num2str(rand_guard_required) '_' num2str(rand_target_required) '_' num2str(average_generator) '.mat'];
             
             time_script_auto_concorde(rand_guard_required, rand_target_required, average_generator) = time_auto_for_struct(1).concorde_time;
             time_script_auto_total(rand_guard_required, rand_target_required, average_generator) = time_auto_for_struct(1).total_time;
             nodes_pushed_tsp(rand_guard_required, rand_target_required, average_generator) = nodes_totsp;
             
             save(filename, 'outfin_sol', 'outfin_cost', 'Out_solName', 'Out_sol', 'G_init', 'edges_totsp', 'nodes_totsp', 'time_auto_for_struct', 'visibility_adjacency_matrix', 'guard_target_struct','environment');
           
             % counter_struct = counter_struct+1;
        end
    end 
    
    
end



%%
            %     for i = 1:length(guards_x)
            % 
            %         node_nam{i} = sprintf('  V%d',i);
            % 
            %     end
            %     
            %     figure(counter_struct);
            %     
            %     plot3(guards_x, guards_y, 0.4*ones(1,length(guards_y)),
            %     'o','Markersize',5,'MarkerEdgeColor','k','MarkerFaceColor','k');
            %     text(guards_x,guards_y, 0.4*ones(1,length(guards_y)),
            %     node_nam,'HorizontalAlignment','left','FontSize',8);
            % 
            % 
            %     for i = 1:length(targets_x)
            % 
            %         target_nam{i} = sprintf(' T%d',i);
            % 
            %     end
            % 
            %     plot3(targets_x, targets_y, 0.4*ones(1,length(targets_x)),
            %     's','Markersize',5,'MarkerEdgeColor','r','MarkerFaceColor','r');
            %     text(targets_x,targets_y, 0.4*ones(1,length(targets_x)),
            %     target_nam,'HorizontalAlignment','left','FontSize',13);






%%

%    figure(1);
% % plot3( [guards(i,1) guards(j,1)], [guards(i,2) guards(j,2)], 0.3*[1 1], 'k', 'LineWidth', 0.5 , 'LineStyle' , '-' );
% 
%     for i = 1:length(s_whole)
%         plot3( [guards(s_whole(i),1) guards(t_whole(i),1)], [guards(s_whole(i),2) guards(t_whole(i),2)], 0.3*[1 1], 'b', 'LineWidth', 3 , 'LineStyle' , '-' );
%         
%     end
%     
%     
%     
%     
%     
%     %%
%     
% %     %plotting
% %         figure;
% %     LWidths2 = 5*(G_gadget2.Edges.Weight+1)/max(G_gadget2.Edges.Weight);
% %     p_gadget2 = plot(G_gadget2, 'Layout','Circle','EdgeLabel', G_gadget2.Edges.Weight, 'LineWidth',LWidths2);
% %     p_gadget2.NodeColor = 'r';
%   %  highlight(p_gadget2,s_sol,t_sol, 'NodeColor','g','EdgeColor','r'); % if given edges solution is wrong this will not highlight the solution
%     %%
%     % plot spanning tree
%      G_gadget2_T = minspantree(G_gadget2);
%      figure;
% %    % LWidths2 = 5*(G_gadget2.Edges.Weight+1)/max(G_gadget2.Edges.Weight);
%      p_gadget2_T = plot(G_gadget2_T, 'NodeLabel',G_gadget2_T.Nodes.Name,'Layout','force','EdgeLabel', G_gadget2_T.Edges.Weight); %, 'LineWidth',LWidths2);
%      p_gadget2_T.NodeColor = 'r';
%     %%
%     
%         %%
%     % following we are doing to highlight the solution
%     s_sol = zeros(G_gadget2.numnodes,1);
%     t_sol = zeros(G_gadget2.numnodes,1);
% 
%     for i = 2:G_gadget2.numnodes
%         s_sol(i-1) = Out_sol(i-1);
%         t_sol(i-1) = Out_sol(i);
%     end
% 
%     s_sol(G_gadget2.numnodes) = t_sol(end-1);
%     t_sol(G_gadget2.numnodes) = s_sol(1);
%     
%     figure;
%     LWidths2 = 5*(G_gadget2.Edges.Weight+1)/max(G_gadget2.Edges.Weight);
%     p_gadget2_n = plot(G_gadget2, 'Layout','Circle','EdgeLabel', G_gadget2.Edges.Weight, 'LineWidth',LWidths2);
%     %p_gadget2_n  = plot(G_gadget2, 'XData', p_gadget2_T.XData, 'YData', p_gadget2_T.YData, 'EdgeLabel', G_gadget2.Edges.Weight, 'LineWidth',LWidths2);
%     
%      highlight(p_gadget2_n,s_sol,t_sol, 'NodeColor','g','EdgeColor','r');







