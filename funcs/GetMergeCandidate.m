function [ lsFinalCandL, lsFinalCandR, llFinalCandL, llFinalCandR ] = GetMergeCandidate( ...
    cur_ls, idxSource, ls_est, ll_inst, ll_idxSet, ll_type, bAscend, level_t, bDraw )
    %GETMERGECANDIDATE This function returns a set of indices of linelets
    %   which satisfy the mering criterion wrt the current line
    %   segment. More specifically, a set within a perpendicular distance,
    %   angle difference, and central distance from the source line
    %   segment. Line segments are considered first otherwise linelets.
    %   When there is a redundancy between a line segment and linelet, the 
    %   linelet is ignored. On the otherhand, when only a linelet -- which 
    %   composes other line segment -- is detected, the linelet is excluded
    %   from the original line segment. However, the last procedure is not
    %   a final one, thus, should be considered again. !!!!!!!!!!!!!!!!!!!
    
    newSetL = [];         newSetR = [];
    lsFinalCandL = [];    lsFinalCandR = [];
    llFinalCandL = [];    llFinalCandR = [];
    lsOrgL = [];          lsOrgR = [];
    
    if strcmp(ll_type, 'ver')
        majorAxisIdx1 = 2;
        majorAxisIdx2 = 4;
        minorAxisIdx1 = 1;
        minorAxisIdx2 = 3;
    else
        majorAxisIdx1 = 1;
        majorAxisIdx2 = 3;
        minorAxisIdx1 = 2;
        minorAxisIdx2 = 4;
    end
    
    if bAscend
        minor_diff_signR = 1;
    else
        minor_diff_signR = -1;
    end
    
    tau_ang = pi/4;
    tau_perpdist = sqrt(2);
    
    num_set = size( ls_est, 1 );
    num_inst = size( ll_inst, 1 );
    instMid1 = [sum( ll_inst(:, 1:2:3), 2 ) sum( ll_inst(:, 2:2:4), 2 )]/2;
    instMid2 = [ll_inst(:, 1) ll_inst(:, 2)];
    instMid3 = [ll_inst(:, 3) ll_inst(:, 4)];
    
    
    % Linelet of current line
    llCur = ll_inst(ll_idxSet{idxSource},:);    
   
    % ------------------------
    % Find candidates from line segemnts 
    % ------------------------
    % 1. angle tolerance
    lsCand_LsAngle = (bAngleAligned( cur_ls(3), ls_est(:,3), tau_ang ));
    
    % 2. vector angle
    dxy = ls_est(:,1:2) - repmat( cur_ls(1:2), num_set, 1 );
    angtmp = atan2(dxy(:,2), dxy(:,1));
    lsCand_VecAngle = (bAngleAligned( cur_ls(3), angtmp, tau_ang ));
    
    % 3. Perpenduciular distance
    % For both horizontal and vertical line segment/linelet, distance on the minor axis is considered, for
    % diagonal case, perpendicular distance is considered
    lsPt1 = ls_est(:,1:2) + repmat( ls_est(:,4)/2, 1, 2 ).*[cos( ls_est(:,3) ) sin( ls_est(:,3) )];
    lsPt2 = ls_est(:,1:2) - repmat( ls_est(:,4)/2, 1, 2 ).*[cos( ls_est(:,3) ) sin( ls_est(:,3) )];
    
    
    perpDist1 = GetPerpDist( cur_ls(1:2), ls_est(:,1:2), cur_ls(3) );
    perpDist2 = GetPerpDist( cur_ls(1:2), lsPt1, cur_ls(3) );
    perpDist3 = GetPerpDist( cur_ls(1:2), lsPt2, cur_ls(3) );
    
    tmp = (perpDist1 + perpDist2 + perpDist3)'/3;
    lsCand_PerpDist = (tmp <= tau_perpdist);
    
    % 4. Linelet difference
    lsCand_LlLen = abs( ls_est(:,5) - cur_ls(5) ) <= max(5, ceil(cur_ls(5) * .5)); 
    
    % 5. Central distance
    dcentLs = ls_est(:,1:2) - repmat( cur_ls(1:2), num_set, 1 );
    dcentLs = sqrt( sum( dcentLs.^2, 2 ) );
    dbtLs = dcentLs - ls_est(:,4)/2 - cur_ls(4)/2;
    btw_ls_len = min([cur_ls(4)*.25, level_t+2]);
    
    lsCand_CentDist = abs(dbtLs) <= btw_ls_len & dbtLs > 0;
        
    lsCand_sofar = lsCand_LsAngle & lsCand_PerpDist & lsCand_CentDist;
    % Leftward
    lsFinalCandL = find(lsCand_sofar & (ls_est(:,1) - cur_ls(1) < 0));
    lsFinalCandL(lsFinalCandL == idxSource) = [];    
    % pick the closest one
    if ~isempty(lsFinalCandL)        
        if length( lsFinalCandL ) > 1
            [~, idxNear] = min(dbtLs(lsFinalCandL));
            lsFinalCandL = lsFinalCandL(idxNear);
        end
    end   
    
    % Rightward
    lsFinalCandR = find(lsCand_sofar & (ls_est(:,1) - cur_ls(1) > 0));
    lsFinalCandR(lsFinalCandR == idxSource) = [];    
    % pick the closest one
    if ~isempty(lsFinalCandR)        
        if length( lsFinalCandR ) > 1
            [~, idxNear] = min(dbtLs(lsFinalCandR));
            lsFinalCandR = lsFinalCandR(idxNear);
        end
    end       
    % ------------------------
    
    % ------------------------
    % Find candidates from linelets
    % ------------------------
    llCur_MajorL = min(llCur(:,majorAxisIdx1)); % minimum position on major axis
    llCur_MajorH = max(llCur(:,majorAxisIdx2)); % maximum position on major axis
    llCur_MinorL = min(llCur(:,minorAxisIdx1)); % minimum position on minor axis
    llCur_MinorH = max(llCur(:,minorAxisIdx2)); % maximum position on minor axis
    
    % Distance
    %llCand3 = abs( ll_inst(:,5) - cur_ls(5) ) < 3; 
    tmp = ll_inst(:,majorAxisIdx2) - llCur_MajorL;
    llCand_DistMajorL = tmp >= 0 & tmp <= cur_ls(5);
    
    tmp = ll_inst(:,majorAxisIdx1) - llCur_MajorH;
    llCand_DistMajorR = tmp >= 0 & tmp <= cur_ls(5);    
    
    tmp = ll_inst(:,minorAxisIdx1) - llCur_MinorL;
    llCand_DistMinorL = tmp == 0 | tmp == -minor_diff_signR;
    
    tmp = ll_inst(:,minorAxisIdx2) - llCur_MinorH;
    llCand_DistMinorR = tmp == 0 | tmp == minor_diff_signR;
    
    % 
    perpDist1 = GetPerpDist(cur_ls(1:2), instMid1, cur_ls(3));
    perpDist2 = GetPerpDist(cur_ls(1:2), instMid2, cur_ls(3));
    perpDist3 = GetPerpDist(cur_ls(1:2), instMid3, cur_ls(3));
    perpDist = (abs( perpDist1 ) + abs( perpDist2 ) + abs( perpDist3 )) / 3;
    
    llCand_PerpDist = abs( perpDist )' < tau_perpdist;
    
    % 
    tmp = instMid1 - repmat( cur_ls(1:2), num_inst, 1 );
    tmp = sqrt( sum( tmp.^2, 2 ) );
    dbtLl = tmp - ll_inst(:,5)/2 - cur_ls(4)/2;
    llCand_CentDist = dbtLl <= min( repmat( cur_ls(5), num_inst, 1 ), ll_inst(:, 5) ) & dbtLl >= 0;
    
    % Leftward linelet candidate
    llCand_Left = (ll_inst(:,3) - cur_ls(1)) < 0;
    llFinalCandL = find( llCand_DistMajorL & llCand_DistMinorL & llCand_PerpDist & llCand_CentDist & llCand_Left );
    % pick the closest one        
    if ~isempty(llFinalCandL)
        if length( llFinalCandL ) > 1
            [~, idxNear] = min(dbtLl(llFinalCandL));
            llFinalCandL = llFinalCandL(idxNear);
        end
    end
    
    % Rightward linelet candidate
    llCand_Right = (ll_inst(:,3) - cur_ls(1)) > 0;
    llFinalCandR = find( llCand_DistMajorR & llCand_DistMinorR & llCand_PerpDist & llCand_CentDist & llCand_Right );
    % pick the closest one        
    if ~isempty(llFinalCandR)
        if length( llFinalCandR ) > 1
            [~, idxNear] = min(dbtLl(llFinalCandR));
            llFinalCandR = llFinalCandR(idxNear);
        end
    end
    % ------------------------
    
    if 0
        % ------------------------
    % Find candidates, from line segemnts and linelets (then remove redundant ones)
    % 1. angle tolerance
    lsCand_LsAngle = (bAngleAligned( cur_ls(3), ls_est(:,3), tau_ang ));
    
    % 2. vector angle
    dxy = ls_est(:,1:2) - repmat( cur_ls(1:2), num_set, 1 );
    angtmp = atan2(dxy(:,2), dxy(:,1));
    lsCand_VecAngle = (bAngleAligned( cur_ls(3), angtmp, tau_ang ));
    
    % 3. Perpenduciular distance
    % For both horizontal and vertical line segment/linelet, distance on the minor axis is considered, for
    % diagonal case, perpendicular distance is considered
    % 4-1. line segment    
    n = [-sin( cur_ls(3) ), cos( cur_ls(3) )]; % unit normal vector of the line
    perpDist = GetPerpDist( cur_ls(1:2), ls_est(:,1:2), cur_ls(3) );
    lsCand_PerpDist = (perpDist' < tau_perpdist);
    % 3. Linelet length
    % !!!!!!!!!!!!!!
    % This should be done in more elegant way, not just checking the diffrence between the length of the
    % longest linelet of each, but doing a kind of imagenary linelet generation
    lsCand_LlLen = abs( ls_est(:,5) - cur_ls(5) ) <= max(3, ceil(cur_ls(5) * .3)); 
%     lsCand_LlLen = true(num_set, 1 );
    
    llCur = ll_inst(ll_idxSet{idxSource},:);    
    llCur_MajorL = min(llCur(:,majorAxisIdx1));
    llCur_MajorH = max(llCur(:,majorAxisIdx2));
    llCur_MinorL = min(llCur(:,minorAxisIdx1));
    llCur_MinorH = max(llCur(:,minorAxisIdx2));
    
    llCand3 = abs( ll_inst(:,5) - cur_ls(5) ) < 3; 
    llCand3LMajor = abs( ll_inst(:,majorAxisIdx2) - llCur_MajorL ) <= cur_ls(5);
    llCand3RMajor = abs( ll_inst(:,majorAxisIdx1) - llCur_MajorH ) <= cur_ls(5);    
    
    llCand3LMinor = ll_inst(:,minorAxisIdx1) - llCur_MinorL;
    llCand3LMinor = llCand3LMinor == 0 | llCand3LMinor == -minor_diff_signR;
    llCand3RMinor = ll_inst(:,minorAxisIdx2) - llCur_MinorH;
    llCand3RMinor = llCand3RMinor == 0 | llCand3RMinor == minor_diff_signR;
    
    
    
    % 4-2. linelet
    perpDist1 = n*(instMid1 - repmat( cur_ls(1:2), num_inst, 1 ))';
    perpDist2 = n*(instMid2 - repmat( cur_ls(1:2), num_inst, 1 ))';
    perpDist3 = n*(instMid3 - repmat( cur_ls(1:2), num_inst, 1 ))';
    perpDist = (abs( perpDist1 ) + abs( perpDist2 ) + abs( perpDist3 )) / 3;
    
    llCand4 = abs( perpDist )' < tau_perpdist;
    
    % 5. Central distance
    % 5-1. line segemnts
    dcentLs = ls_est(:,1:2) - repmat( cur_ls(1:2), num_set, 1 );
    dcentLs = sqrt( sum( dcentLs.^2, 2 ) );
    dbtLs = dcentLs - ls_est(:,4)/2 - cur_ls(4)/2;
    lsCand_CentDist = dbtLs <= ( repmat( cur_ls(5), num_set, 1 ) + ls_est(:, 5) )/2;% & dbtLs > 0;
    
    % 5-2. linelet
    dcentLl = instMid1 - repmat( cur_ls(1:2), num_inst, 1 );
    dcentLl = sqrt( sum( dcentLl.^2, 2 ) );
    dbtLl = dcentLl - ll_inst(:,5)/2 - cur_ls(4)/2;
    llCand5 = dbtLl <= min( repmat( cur_ls(5), num_inst, 1 ), ll_inst(:, 5) ) & dbtLl >= 0;
    llCand6 = ll_inst(:, 6) == cur_ls(6);
    
    % wrt the target direction
    % Leftward line segment candidate
    lsCand5L = ls_est(:,1) - cur_ls(1) <= 0;
    lsFinalCandL = find( lsCand_LsAngle & lsCand_LlLen & lsCand_PerpDist & lsCand_CentDist & lsCand5L );
    
    if ~isempty(lsFinalCandL)
        % pick the closest one
        if length( lsFinalCandL ) > 1
            [~, idxNear] = min(dbtLs(lsFinalCandL));
            lsFinalCandL = lsFinalCandL(idxNear);
        end
        lsOrgL = ll_idxSet{lsFinalCandL};
        lsFinalCandL(lsFinalCandL == idxSource) = [];
    end    
    
    % Leftward linelet candidate
    llCand5L = (ll_inst(:,3) - cur_ls(1)) <= 0;
    llFinalCandL = find( llCand3LMajor & llCand3LMinor & llCand4 & llCand5 & llCand5L );
    if ~isempty(llFinalCandL)
        % pick the closest one
        if length( llFinalCandL ) > 1
            [~, idxNear] = min(dbtLl(llFinalCandL));
            llFinalCandL = llFinalCandL(idxNear);
        end
        if ~isempty( intersect( lsOrgL, llFinalCandL ) )
%             llFinalCandL = [];
        end
    end
    newSetL = lsOrgL;
    
    %% Rightward candidate
    lsCand5R = ls_est(:,1) - cur_ls(1) >= 0;
    lsFinalCandR = find( lsCand_LsAngle & lsCand_LlLen & lsCand_PerpDist & lsCand_CentDist & lsCand5R );
        
    if ~isempty(lsFinalCandR)
        % pick the closest one
        if length( lsFinalCandR ) > 1
            [~, idxNear] = min(dbtLs(lsFinalCandR));
            lsFinalCandR = lsFinalCandR(idxNear);
        end
        lsOrgR = ll_idxSet{lsFinalCandR};
        lsFinalCandR(lsFinalCandR == idxSource) = [];
    end
    % Rightward candidate
    llCand5R = (ll_inst(:,1) - cur_ls(1)) >= 0;
    llFinalCandR = find( llCand3RMajor & llCand3RMinor & llCand4 & llCand5 & llCand5R );
    if ~isempty(llFinalCandR)
        % pick the closest one
        if length( llFinalCandR ) > 1
            [~, idxNear] = min(dbtLl(llFinalCandR));
            llFinalCandR = llFinalCandR(idxNear);
        end
        if ~isempty( intersect( lsOrgR, llFinalCandR ) )
%             llFinalCandR = [];
        end
    end    
    newSetR = lsOrgR;
    end
    
    %% Debug display
    if bDraw        
        % Original line segmnet
%         figure(21); imshow(im_gray); hold on;
        x1 = ls_est(idxSource,1:2) + ls_est(idxSource,4)/2*[cos(ls_est(idxSource,3)) sin(ls_est(idxSource,3))];
        x2 = ls_est(idxSource,1:2) - ls_est(idxSource,4)/2*[cos(ls_est(idxSource,3)) sin(ls_est(idxSource,3))];
        plot( [x1(1) x2(1)], [x1(2) x2(2)], 'color', [1 0 0], 'linewidth', 2 );
        plot( ls_est(idxSource,1), ls_est(idxSource,2), 'ro' );
      
        idxTarget = lsFinalCandL;
        for k = 1:length(idxTarget)
            x1 = ls_est(idxTarget(k),1:2) + ls_est(idxTarget(k),4)/2*[cos(ls_est(idxTarget(k),3)) sin(ls_est(idxTarget(k),3))];
            x2 = ls_est(idxTarget(k),1:2) - ls_est(idxTarget(k),4)/2*[cos(ls_est(idxTarget(k),3)) sin(ls_est(idxTarget(k),3))];
            plot( [x1(1) x2(1)], [x1(2) x2(2)], 'color', [0 1 0], 'linewidth', 1 ); % updated
            plot(ls_est(idxTarget(k),1), ls_est(idxTarget(k),2), 'go');
        end

        idxTarget = lsFinalCandR;
        for k = 1:length(idxTarget)
            x1 = ls_est(idxTarget(k),1:2) + ls_est(idxTarget(k),4)/2*[cos(ls_est(idxTarget(k),3)) sin(ls_est(idxTarget(k),3))];
            x2 = ls_est(idxTarget(k),1:2) - ls_est(idxTarget(k),4)/2*[cos(ls_est(idxTarget(k),3)) sin(ls_est(idxTarget(k),3))];
            plot( [x1(1) x2(1)], [x1(2) x2(2)], 'color', [0 0 1], 'linewidth', 1 ); % updated
            plot(ls_est(idxTarget(k),1), ls_est(idxTarget(k),2), 'bo');
        end
        title( ['Current index: ' num2str(idxSource)] );
    end
end

