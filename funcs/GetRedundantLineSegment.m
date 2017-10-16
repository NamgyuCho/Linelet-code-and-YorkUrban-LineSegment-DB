function [idxRedun] = GetRedundantLineSegment( cur_ls, ls_est )
    %GETREDUNDANTLINESEGMENT Summary of this function goes here
    %   Detailed explanation goes here
    tau_ang = pi/8;
    tau_perpdist = sqrt(2);
    
    num_set = size( ls_est, 1 );
    
    % ------------------------
    % Find candidates from line segemnts 
    % ------------------------
    % 1. angle tolerance
    lsCand_LsAngle = (bAngleAligned( cur_ls(3), ls_est(:,3), tau_ang ));
    
    % 2. Perpenduciular distance
    % For both horizontal and vertical line segment/linelet, distance on the minor axis is considered, for
    % diagonal case, perpendicular distance is considered
    lsPt1 = ls_est(:,1:2) + repmat( ls_est(:,4)/2, 1, 2 ).*[cos( ls_est(:,3) ) sin( ls_est(:,3) )];
    lsPt2 = ls_est(:,1:2) - repmat( ls_est(:,4)/2, 1, 2 ).*[cos( ls_est(:,3) ) sin( ls_est(:,3) )];
        
    perpDist1 = GetPerpDist( cur_ls(1:2), ls_est(:,1:2), cur_ls(3) );
    perpDist2 = GetPerpDist( cur_ls(1:2), lsPt1, cur_ls(3) );
    perpDist3 = GetPerpDist( cur_ls(1:2), lsPt2, cur_ls(3) );
    
    tmp = (perpDist1 + perpDist2 + perpDist3)'/3;
    lsCand_PerpDist = (tmp <= tau_perpdist);
    
    % 3. Central distance
    dcentLs = ls_est(:,1:2) - repmat( cur_ls(1:2), num_set, 1 );
    dcentLs = sqrt( sum( dcentLs.^2, 2 ) );
    dbtLs = dcentLs - ls_est(:,4)/2 - cur_ls(4)/2;
    lsCand_CentDist = dbtLs <= cur_ls(5) & dbtLs < 0;
        
    idxRedun = find(lsCand_LsAngle & lsCand_PerpDist & lsCand_CentDist);    
    % ------------------------
end

