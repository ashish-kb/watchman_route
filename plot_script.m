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
        
 %       tic 
%for i= 1:10
    
node_vec = cellfun(@(x) strsplit(x(2:(regexp(x,'-','start')-1)),'-') , G_comp.Edges.EndNodes,'uni',0);
%cell2mat(cellfun(@(x) str2double(strsplit(x((regexp(x,'-','start')+1):end),'-')) , G_comp.Edges.EndNodes(:,1),'uni',0));
%end
%toc