function conn_map = DisconnectLinelet(conn_map, idx1, idx2)
    conn_map(idx1, idx2) = 0;
    conn_map(idx2, idx1) = 0;
end
