function [ll_cand_new, conn_map_new] = DecomposeConnMap(ll_cand, conn_map)

% Decompose the input linlet let set according to their connection state
n_ll = 1;
ll_usage = false(size(ll_cand,1), 1);
ll_cand_new = cell(1);
conn_map_new = cell(1);

for i = 1:size(ll_cand,1)
    if ll_usage(i), continue; end
    
    stack_ll = i;
    ll_tmp = [];
    while ~isempty(stack_ll) 
        cur_ll = stack_ll(1);   stack_ll(1) = []; % pop 
        ll_usage(cur_ll) = true;      
    
        % Search connected linelets
        idx_nb = find(conn_map(cur_ll,:) ~= 0 & ~ll_usage' );
                 
        stack_ll = [stack_ll idx_nb];
        ll_tmp = [ll_tmp; ll_cand(cur_ll,:)];
    end
    
    if isempty(ll_tmp), continue; end
    
    [~, idx_srt] = sort(ll_tmp(:,1));
    ll_tmp = ll_tmp(idx_srt,:);
    
    ll_cand_new{n_ll,1} = ll_tmp;
    n_ll = n_ll + 1;
end

for i=1:n_ll-1
    
    conn_map_new{i,1} = Get_ConnMap(ll_cand_new{i});
end
