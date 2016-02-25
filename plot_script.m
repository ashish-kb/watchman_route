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
            for i = 1:length(guard_target_struct(1).guards_x)
        
                node_nam{i} = sprintf('  V%d',i);
        
            end
            
            figure(1);
            
            plot3(guard_target_struct(1).guards_x, guard_target_struct(1).guards_y, 0.4*ones(1,length(guard_target_struct(1).guards_y)), 'o','Markersize',5,'MarkerEdgeColor','k','MarkerFaceColor','k');
            text(guard_target_struct(1).guards_x,guard_target_struct(1).guards_y, 0.4*ones(1,length(guard_target_struct(1).guards_y)), node_nam,'HorizontalAlignment','left','FontSize',8);
        
        
            for i = 1:length(guard_target_struct(1).targets_x)
        
                target_nam{i} = sprintf(' T%d',i);
        
            end
        
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

 Cutting on node 431
OVERFLOW in CCbigguy_addmult (4)
BIGGUY errors are fatal
FATAL ERROR - received signal SIGABRT (6/6)
sleeping 1 more hours to permit debugger access