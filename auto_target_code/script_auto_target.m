clear all; close all; clc;

addpath(genpath('/home/ashishkb/Dropbox/MATLAB/RAASLAB/gtsp/solved_cases/'));
addpath(genpath('/home/ashishkb/Dropbox/MATLAB/RAASLAB/gtsp/source_code/'));





format long;


%Robustness constant
epsilon = 0.000000001;

load('for_random_target.mat');

num_guards = 63;
max_targets = 15;
min_targets = 3;
max_average = 1; % we can increase this if we want to average the data on one machine
% but we will run on multiple machines(10) so each one
% runs 1
nodes_pushed_tsp = zeros(max_targets, 1);
time_script_auto_concorde = zeros(max_targets, max_average);%struct('total_time',[],'concorde_time',[]);
time_script_auto_total = zeros(max_targets, max_average);


%%
counter_struct = 1; % not used anywhere just kept will remove later(or not)

for rand_target_required = 3:1:max_targets
    for average_generator = 1:max_average
        guard_target_struct = struct('guards_x',{},'guards_y',{}, 'targets_x',...
            {}, 'targets_y', {}, 'bot_x', {}, 'bot_y', {}, 'environ_guardX',...
            {}, 'environ_guardY', {}, 'guards', {});
        
         guard_target_struct(1).environ_guardX = environ_guardX;
         guard_target_struct(1).environ_guardY = environ_guardY;
        
        
         target_randX = zeros(rand_target_required,1);
         target_randY = zeros(rand_target_required,1);
         
         rng('shuffle');
        
         got_target_rands = 1; % checks if right number of target have been
                               % generated 

          while (got_target_rands <= rand_target_required)
                % a stupid way to generate random numbers
                X_try =  randi([ceil(environment_min_x) floor(environment_max_x)]...
                    ,1 ,'double') + rand(1, 'double') ; 
                
                Y_try =  randi([ceil(environment_min_y) floor(environment_max_y)]...
                    ,1 ,'double') + rand(1, 'double') ;

                if(in_environment([X_try Y_try] , environment , epsilon))
                    target_randX(got_target_rands) = X_try;
                    target_randY(got_target_rands) = Y_try;

                    got_target_rands  = got_target_rands + 1;
                end

          end
          
          guard_target_struct(1).guards = [[guards_x, guards_y];[target_randX,...
              target_randY]; [bot_x, bot_y]];  

          guard_target_struct(1).targets_x = target_randX;
          guard_target_struct(1).targets_y = target_randY;
          
          guard_target_struct(1).guards_x = [guards_x];
          guard_target_struct(1).guards_y = [guards_y];

          guard_target_struct(1).bot_x = [bot_x];
          guard_target_struct(1).bot_y = [bot_y];

               
          
          visibility_adjacency_matrix = ...
              visibility_graph(guard_target_struct(1).guards, environment, epsilon);
          
          filepath_targets = sprintf('target_%d.txt', rand_target_required);
          
          fid = fopen(filepath_targets,'W');
          if fid < 0
              error('Cannot create  file');
              return;
          end
          
          fprintf(fid,'//target x-y coordinates \n');
          
          for i = 1:length(target_randX)
              fprintf(fid,'%d %d\n', target_randX(i), target_randY(i));
          end

          fclose(fid);
                    
          [outfin_sol, outfin_cost, Out_solName, whole_path_nodes, G_init,...
              G_nodebot, edges_totsp, nodes_totsp, time_concorde_struct,...
              time_auto_for_struct]...              
              = auto_for_content_Multi_noon_bean(visibility_adjacency_matrix,...
              guard_target_struct, counter_struct);
          
          
          filename = [num2str(rand_target_required)...
              '_' num2str(average_generator) '.mat'];
          
          
          time_script_auto_concorde(rand_target_required, average_generator)...
              = time_auto_for_struct(1).concorde_time;
          
          time_script_auto_total(rand_target_required, average_generator)...
              = time_auto_for_struct(1).total_time;
          
          nodes_pushed_tsp(rand_target_required, average_generator)...
              = nodes_totsp;
            
          
          save(filename, 'outfin_sol', 'outfin_cost', 'Out_solName', 'whole_path_nodes',...
              'G_init', 'G_nodebot','edges_totsp', 'nodes_totsp', 'time_concorde_struct'...
              ,'time_auto_for_struct', 'visibility_adjacency_matrix',...
              'guard_target_struct','environment');
        
        
    end
end