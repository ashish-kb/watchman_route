num_nodes = 63;

s  = ones(num_nodes,1);
t  = ones(num_nodes,1);
weigh  = ones(num_nodes,1);


s = [1:num_nodes]';
t(1:(num_nodes-1)) =  [2:num_nodes]';

t(num_nodes) = s(1);

%%

tot_edgeweight = sum(weigh);
gadget2_compMat = (tot_edgeweight*50+100)*ones(num_nodes, num_nodes) - diag((tot_edgeweight*50+99)*ones(1,num_nodes));



for i = 1:num_nodes
    
     gadget2_compMat(s(i),t(i))  =  weigh(i);
     gadget2_compMat(t(i), s(i))  =  weigh(i);
        
     
    
end

%s[1;2;3;2;5]
%t[4;3;4;5;6]