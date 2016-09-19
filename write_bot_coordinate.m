for i = 1:length(outfin_sol)
    cur_list = outfin_sol{i};
    
    if(~isempty(cur_list))
        cur_list_unique = [cur_list(1) unique(cur_list(2:(end-1)), 'stable') cur_list(end)];
        coord_mat = zeros(length(cur_list_unique), 2);
        init_bot_loc = -1*cur_list_unique(1);        
        fin_bot_loc = -1*cur_list_unique(end);    
        coord_mat(1,:) = [bot_x(init_bot_loc) bot_y(init_bot_loc)];
        coord_mat(end,:) = [bot_x(fin_bot_loc) bot_y(fin_bot_loc)];
        coord_mat(2:(end-1),:) = [guards_x(cur_list_unique(2:(end-1))) guards_y(cur_list_unique(2:(end-1)))];
    else
        coord_mat = [];
    end
    
    in_filename = sprintf('robot%d_goals.txt', i);

    
    fid = fopen(in_filename,'W');
    if fid < 0
        error('Cannot create temp file');
        return;
    end
    
    if(~isempty(coord_mat))
        fprintf(fid,'%f\t%f\n', coord_mat(1:end-1,:)');
        fprintf(fid,'%f\t%f', coord_mat(end,:)');
    end

    fclose(fid);

    
    
    
end