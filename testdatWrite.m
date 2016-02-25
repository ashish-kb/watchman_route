
node_s = G_gadget2.numnodes;
edge_s = G_gadget2.numedges;
fid = fopen('tempdat.tsp','W');
if fid < 0
	error('Cannot create temp file');
	return;
end


[s_gadget2 t_gadget2] = findedge(G_gadget2);
weights_mat = G_gadget2.Edges.Weight(:);

G_gadget2_edges = [s_gadget2(:) t_gadget2(:) weights_mat(:)];

%%
fprintf(fid, [num2str(node_s) ' ' num2str(edge_s) '\n']);


%%
% fprintf(fid, [repmat('%d %d %d', 1, size(s_gadget2,1)) '\n'], s_gadget2(:)', t_gadget2(:)', weights_mat(:)');
% fprintf(fid, [repmat('%f\t', 1, size(G_gadget2_edges, 2)) '\n'], G_gadget2_edges);
 fprintf(fid,'%d %d %d\n',G_gadget2_edges(:,1:3).');

fclose(fid);

%%
[s_gadget2 t_gadget2] = findedge(G_gadget2); % s - start -> t- target this is a list of edges in terms of index


    tot_edgeweight = sum(G_gadget2.Edges.Weight);
    gadget2_compMat = (tot_edgeweight+100)*ones(G_gadget2.numnodes, G_gadget2.numnodes) - diag((tot_edgeweight+99)*ones(1,G_gadget2.numnodes));

%     gadget2_compMat1 = gadget2_compMat;

%     for i = 1:G_gadget2.numedges
% 
%          gadget2_compMat(s_gadget2(i),t_gadget2(i))  =  G_gadget2.Edges.Weight(i);
%          gadget2_compMat(t_gadget2(i), s_gadget2(i))  =  G_gadget2.Edges.Weight(i);
% 
% 
% 
%     end
    ind_gadget2 = sub2ind(size(gadget2_compMat), s_gadget2(:), t_gadget2(:));
    ind_gadget2 = [ind_gadget2;sub2ind(size(gadget2_compMat), t_gadget2(:), s_gadget2(:))];
    
    weight_feed_mat = G_gadget2.Edges.Weight(:); % will feed this below
    gadget2_compMat(ind_gadget2(:)) = [weight_feed_mat; weight_feed_mat];