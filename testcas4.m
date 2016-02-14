% Open the temporary tsplib format file

% Open the temporary tsplib format file
%
%     'EDGE_DATA_FORMAT : EDGE_LIST\n' ...'NODE_COORD_TYPE : TWOD_COORDS
%    'EDGE_WEIGHT_FORMAT :  WEIGHT_LIST\n' ...



fid = fopen('temp.tsp','W');
if fid < 0
	error('Cannot create temp file');
	return;
end



fprintf(fid, ['NAME : temp\n' ...
							'TYPE : TSP\n' ...
							'COMMENT : temp\n' ...
							'DIMENSION : ' num2str(4) '\n' ...
							'EDGE_WEIGHT_TYPE : EUC_2D\n' ...
                            'EDGE_DATA_FORMAT : EDGE_LIST\n' ... 
							'NODE_COORD_TYPE : TWOD_COORDS\n' ...
							]);
                        
                        %'EDGE_DATA_FORMAT : EDGE_LIST\n' ... 

% for i = 1 : n
% 	fprintf(fid, [num2str(i-1) ' ' num2str(scaled_X(i,1)) ' ' num2str(scaled_X(i,2)) '\n']); 
% end
fprintf(fid, ['EDGE_DATA_SECTION : \n']); 
fprintf(fid, ['0' ' ' '1' '\n']); 
fprintf(fid, ['1' ' ' '2' '\n']); 
fprintf(fid, ['2' ' ' '3' '\n']); 
fprintf(fid, ['3' ' ' '0' '\n']); 

fprintf(fid, ['-1' '\n']); 

fprintf(fid, ['NODE_COORD_SECTION : \n']); 

fprintf(fid, ['0' ' ' '0' ' ' '0' '\n']); 
fprintf(fid, ['1' ' ' '10' ' ' '0' '\n']); 
fprintf(fid, ['2' ' ' '10' ' ' '40' '\n']); 
fprintf(fid, ['3' ' ' '0' ' ' '50' '\n']); 


fprintf(fid, ['-1' '\n']); 



% fprintf(fid, ['0' ' ' '1' ' ' '10' '\n']); 
% fprintf(fid, ['1' ' ' '2' ' ' '40' '\n']); 
% fprintf(fid, ['2' ' ' '3' ' ' '70' '\n']); 
% fprintf(fid, ['3' ' ' '0' ' ' '50' '\n']);                   

%fprintf(fid, ['-1' '\n']); 

fclose(fid);




cmd = ['/home/ashishkb/softwares/concorde/concorde_build/TSP/concorde' ' -x temp.tsp'];
system(cmd);

%'EDGE_DATA_FORMAT : EDGE_LIST\n' ...
%'EDGE_DATA_SECTION : \n' ...
%test