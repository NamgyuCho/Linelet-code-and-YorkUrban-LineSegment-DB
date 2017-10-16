function [ cons_map ] = Get_ConsistencyMap( ll_cand, conn_map )
%GET_CONSISTENCYMAP Summary of this function goes here
%   Detailed explanation goes here

num_ll = size(ll_cand, 1);

% Build a direction consistency map; how many neighbor has direction consistent neighbors
cons_map = zeros(num_ll, 2); % (num of connected neighbors of left neighbor, right neighbor)

for i = 1:num_ll
    % leftward
    idx_l = find(abs(conn_map(i,:)) == 1);    
    for j = 1:length(idx_l)
        num_const_l_nb = 1;
        curr_nb = idx_l(j);
        cur_sign = conn_map(i, curr_nb);
        next_nb = find( conn_map(curr_nb,:) == cur_sign );
        
        stack_nb = [curr_nb next_nb];
        
        while length(stack_nb) == 2
            num_const_l_nb = num_const_l_nb + 1;
            
            curr_nb = stack_nb(1);   stack_nb(1) = [];
            next_nb = stack_nb(1);   stack_nb(1) = [];
            cur_sign = conn_map(curr_nb, next_nb);
            
            idx_next = find( conn_map(next_nb,:) == cur_sign );
            if ~isempty(idx_next)
                stack_nb = [next_nb idx_next];
            end
        end
        if conn_map(i, idx_l(j)) > 0
            cons_map(i,1) = num_const_l_nb;
        else
            cons_map(i,1) = -num_const_l_nb;
        end
    end
    
    % rightward
    idx_r = find(abs(conn_map(i,:)) == 2);
    for j = 1:length(idx_r)
        num_const_r_nb = 1;
        curr_nb = idx_r(j);
        cur_sign = conn_map(i, curr_nb);
        next_nb = find( conn_map(curr_nb,:) == cur_sign );
        
        stack_nb = [curr_nb next_nb];
        
        while length(stack_nb) == 2
            num_const_r_nb = num_const_r_nb + 1;
            
            curr_nb = stack_nb(1);   stack_nb(1) = [];
            next_nb = stack_nb(1);   stack_nb(1) = [];
            cur_sign = conn_map(curr_nb, next_nb);
            
            idx_next = find( conn_map(next_nb,:) == cur_sign );
            if ~isempty(idx_next)
                stack_nb = [next_nb idx_next];
            end
        end
        if conn_map(i, idx_r(j)) > 0
            cons_map(i,2) = num_const_r_nb;
        else
            cons_map(i,2) = -num_const_r_nb;
        end
    end
end

