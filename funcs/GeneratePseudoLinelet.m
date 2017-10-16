function [ ll_pseudo ] = GeneratePseudoLinelet( ls1, ls2, ll_Set1, ll_Set2, ll_type )
    %GENERATEPSEUDOLINELET Summary of this function goes here
    %   Detailed explanation goes here
    
    % ls1 = cur_ls; ls2 = ls_est(lsCandSet(k),:); ll_Set1 = ll_inst(ll_idxSet{ii},:); ll_Set2 = ll_inst(ll_idxSet{lsCandSet(k)},:);
    
    if strcmp( ll_type, 'hor' )
        major_idx1 = 1;
        major_idx2 = 3;
        minor_idx1 = 2;
        minor_idx2 = 2;
    elseif strcmp( ll_type, 'ver' )
        major_idx1 = 2;
        major_idx2 = 4;
        minor_idx1 = 1;
        minor_idx2 = 1;
    else % diag
        major_idx1 = 1;
        major_idx2 = 3;
        minor_idx1 = 2;
        minor_idx2 = 4;
    end
    
    ll_End1 = GetEndLinelet( ll_Set1 );
    ll_End2 = GetEndLinelet( ll_Set2 );
    
    % Generate an imagenary linelet(s) between two line segments
    if ls1(major_idx1) < ls2(major_idx1)
        % When ls1 is located at the left side
        [xx, yy] = bresenham(ll_End1(end,3), ll_End1(end,4), ll_End2(1,1), ll_End2(1,2));
        pt = [xx yy];
        pt1 = ConvertLinelet2Pts( ll_End2(1,:), ll_type );
        pt2 = ConvertLinelet2Pts( ll_End1(end,:), ll_type );
    else
        % Otherwise
        [xx, yy] = bresenham(ll_End2(end,3), ll_End2(end,4), ll_End1(1,1), ll_End1(1,2));
        pt = [xx yy];        
        pt1 = ConvertLinelet2Pts( ll_End2(end,:), ll_type );
        pt2 = ConvertLinelet2Pts( ll_End1(1,:), ll_type );
    end
    ll_pseudo = ConvertPts2Linelet( [pt; pt1; pt2], ll_type );
end

