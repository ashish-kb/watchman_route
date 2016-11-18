%Clear the desk
clc; figure;
addpath(genpath('/home/ashishkb/Dropbox/MATLAB/RAASLAB/gtsp/solved_cases/'));
addpath(genpath('/home/ashishkb/Dropbox/MATLAB/RAASLAB/gtsp/source_code/'));



filepath_guards = 'ginput2.guards'; % where u want to store that data
filepath_nodetarget = 'ginput2.nodetarget'; % just number of nodes and targets


format long;


%Robustness constant
epsilon = 0.000000001;


%Read environment and guards geometry from files
environment = read_vertices_from_file('../source_code/gazebo_rviz.environment');



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

%message = sprintf('Click nodes; press enter to exit');
%uiwait(msgbox(message));


%[guards_x, guards_y] = ginputc('Color', 'b','ShowPoints', true);%, 'ConnectPoints', true);

% guards_file = read_vertices_from_file('./source_code/guards.environment');
% guards_x = guards_file{1,1}(:,1);
% guards_y = guards_file{1,1}(:,2);

guards = [guards_x, guards_y];


for i = 1:length(guards_x)

    node_nam{i} = sprintf('  V%d',i);

end

plot3(guards(:,1), guards(:,2), 0.4*ones(1,length(guards(:,2))), 'o','Markersize',5,'MarkerEdgeColor','k','MarkerFaceColor','k');
text(guards_x,guards_y, 0.4*ones(1,length(guards_y)), node_nam,'HorizontalAlignment','left','FontSize',8);

message = sprintf('Click targets; press enter to exit');
%uiwait(msgbox(message));
%[targets_x, targets_y] = ginputc('Color', 'r','ShowPoints', true);   %, 'ConnectPoints', true);
% targets = read_vertices_from_file('./source_code/targets.environment');
% targets_x = targets{1,1}(:,1);
% targets_y = targets{1,1}(:,2);

guards = [guards; [targets_x, targets_y]]; % concatenate targets with guards
%guards = read_vertices_from_file('./source_code/example1.guards'); 
%guards = guards{1}; % converts cell to mat
target_nam  = cell(1, length(targets_x));
    
for i = 1:length(targets_x)

    target_nam{i} = sprintf(' T%d',i);

end

plot3(targets_x, targets_y, 0.4*ones(1,length(targets_x)), 's','Markersize',10,'MarkerEdgeColor','r','MarkerFaceColor','r');
text(targets_x,targets_y, 0.4*ones(1,length(targets_x)), target_nam,'HorizontalAlignment','left','FontSize',13);

%Plot Guards
% plot3(guards_x, guards_y, 0.4*ones(1,length(guards_y)), ...
%     'o','Markersize',9,'MarkerEdgeColor','y','MarkerFaceColor','k');
message = sprintf('Click bots; press enter to exit');
uiwait(msgbox(message));
[bot_x, bot_y] = ginputc('Color', 'g','ShowPoints', true);%, 'ConnectPoints', true);
for i = 1:length(bot_x)

    bot_nam{i} = sprintf('  B%d',i);

end

plot3(bot_x, bot_y, 0.4*ones(1,length(bot_x)), '*','Markersize',5,'MarkerEdgeColor','b','MarkerFaceColor','b');
text(bot_x,bot_y, 0.4*ones(1,length(bot_x)), bot_nam,'HorizontalAlignment','left','FontSize',8);

environ_guards = cell2mat(environment(:));
%* guard_target_struct(1).environ_guardX = environ_guards(:,1);
%* guard_target_struct(1).environ_guardY = environ_guards(:,2);


environ_guardX = environ_guards(:,1);
environ_guardY = environ_guards(:,2);
guard_target_struct = struct('guards_x',{},'guards_y',{}, 'targets_x', {}, 'targets_y', {}, 'bot_x', {}, 'bot_y', {}, 'environ_guardX', {}, 'environ_guardY', {}, 'guards', {});
%             nodes_totsp = 0;


guard_target_struct(1).guards = [[guards_x, guards_y];[targets_x, targets_y]; [bot_x, bot_y]];  % guards and targets concatenated
% 
% 
guard_target_struct(1).targets_x = targets_x;
guard_target_struct(1).targets_y = targets_y;
% % 
guard_target_struct(1).guards_x = [guards_x];
guard_target_struct(1).guards_y = [guards_y];

guard_target_struct(1).bot_x = [bot_x];
guard_target_struct(1).bot_y = [bot_y];

guard_target_struct(1).environ_guardX = [environ_guardX]; % ** not there in previous versions 
guard_target_struct(1).environ_guardY = [environ_guardY];




%[guards_x, guards_y] = ginputc('Color', 'b','ShowPoints', true);%, 'ConnectPoints', true);


%Compute and Plot visibility Graph
visibility_adjacency_matrix = visibility_graph(guard_target_struct(1).guards, environment, epsilon);
% for i = 1 : size( visibility_adjacency_matrix , 1 )
%     for j = 1 : size( visibility_adjacency_matrix , 2 )
%         if( visibility_adjacency_matrix(i,j) == 1 )
%             plot3( [guards(i,1) guards(j,1)], [guards(i,2) guards(j,2)], 0.3*[1 1], ...
%                 'g', 'LineWidth', 0.5 , 'LineStyle' , '-' );
%         end
%     end
% end

% for i = 1 : size( guards_x, 1 )
%     for j = 1 : size( guards_x , 1 )
%         if( visibility_adjacency_matrix(i,j) == 1 )
%             plot3( [guards(i,1) guards(j,1)], [guards(i,2) guards(j,2)], 0.3*[1 1], ...
%                 'k', 'LineWidth', 0.0 , 'LineStyle' , '--' );
%         end
%     end
% end

%writing to file

fid = fopen(filepath_guards,'W');
if fid < 0
	error('Cannot create  file');
	return;
end

fprintf(fid,'//Guard x-y coordinates \n');

for i = 1:length(guards)
    fprintf(fid,'%d %d\n', guards(i,1), guards(i,2));
end

fclose(fid);


fid = fopen(filepath_nodetarget,'W');
if fid < 0
	error('Cannot create file');
	return;
end
fprintf(fid, '%d  %d %d', length(guards_x), length(targets_x), length(bot_x));

fclose(fid);
