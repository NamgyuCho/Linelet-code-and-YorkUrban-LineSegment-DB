function conn_map = Get_ConnMap(ll_cand)

conn_map = zeros(size(ll_cand, 1)); % ll_cand = (left, right, y, length, angle)
for i = 1:size(ll_cand,1)
    % Find the nearest neighbors at both sides
    ind_nb = find( ( abs(ll_cand(i,1) - ll_cand(:,2)) == 1 |...
                     abs(ll_cand(i,2) - ll_cand(:,1)) == 1 ) &...
                     abs(ll_cand(:,3) - ll_cand(i,3)) == 1 );
    
    hor_dir = ll_cand(ind_nb,1) - ll_cand(i,1) > 0;
    conn_map(i, ind_nb) = (ll_cand(ind_nb,3) - ll_cand(i,3)) .* (hor_dir + 1);
end
