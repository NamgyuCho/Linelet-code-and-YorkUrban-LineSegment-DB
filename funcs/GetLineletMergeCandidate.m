function [ llFinalCandL, llFinalCandR ] = GetLineletMergeCandidate( ls_est, llCurSet, ll_inst, ll_used, ll_type, bAscend, bDraw )
    %GETLINELETMERGECANDIDATE Summary of this function goes here
    %   Detailed explanation goes here
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
    
    tau_ang = pi/8;
    tau_perpdist = .5;%sqrt(2);
    
    num_set = size( ls_est, 1 );
    num_inst = size( ll_inst, 1 );
    instMid1 = [sum( ll_inst(:, 1:2:3), 2 ) sum( ll_inst(:, 2:2:4), 2 )]/2;
    instMid2 = [ll_inst(:, 1) ll_inst(:, 2)];
    instMid3 = [ll_inst(:, 3) ll_inst(:, 4)];
    
    lsPt1 = ls_est(:,1:2) + repmat( ls_est(:,4)/2, 1, 2 ).*[cos( ls_est(:,3) ) sin( ls_est(:,3) )];
    lsPt2 = ls_est(:,1:2) - repmat( ls_est(:,4)/2, 1, 2 ).*[cos( ls_est(:,3) ) sin( ls_est(:,3) )];
    
    
    % Find candidates, from line segemnts and linelets (then remove redundant ones)
    % 1. ll type
    if strcmp(ll_type, 'hor')
        llCand1 = ll_inst(:, 6) == 1 | ll_inst(:, 6) == 1;
    elseif strcmp(ll_type, 'ver')
        llCand1 = ll_inst(:, 6) == 3 | ll_inst(:, 6) == 3;
    else
        llCand1 = ll_inst(:, 6) == 2;
    end
    
    % 2. vector angle
    % 3. Linelet length    
    llCur_MajorL = min(llCurSet(:,majorAxisIdx1));
    llCur_MajorH = max(llCurSet(:,majorAxisIdx2));
    llCur_MinorL = min(llCurSet(:,minorAxisIdx1));
    llCur_MinorH = max(llCurSet(:,minorAxisIdx2));
    
    llCand3 = abs( ll_inst(:,5) - ls_est(5) ) < 1;
    llCand3LMajor = abs( ll_inst(:,majorAxisIdx2) - llCur_MajorL ) <= ls_est(5);
    llCand3RMajor = abs( ll_inst(:,majorAxisIdx1) - llCur_MajorH ) <= ls_est(5);    
    
    llCand3LMinor = ll_inst(:,minorAxisIdx1) - llCur_MinorL;
    llCand3LMinor = llCand3LMinor == 0 | llCand3LMinor == -minor_diff_signR;
    llCand3RMinor = ll_inst(:,minorAxisIdx2) - llCur_MinorH;
    llCand3RMinor = llCand3RMinor == 0 | llCand3RMinor == minor_diff_signR;
    
    % 4. Perpenduciular distance
    % For both horizontal and vertical line segment/linelet, distance on the minor axis is considered, for
    % diagonal case, perpendicular distance is considered
    % 4-1. line segment
    n = [-sin( ls_est(3) ), cos( ls_est(3) )]; % unit normal vector of the line
    % Perpendicular distance -- line normal * distance vector between one of original line point and estimated one's
    % 4-2. linelet
    perpDist1 = n*(instMid1 - repmat( ls_est(1:2), num_inst, 1 ))';
    perpDist2 = n*(instMid2 - repmat( ls_est(1:2), num_inst, 1 ))';
    perpDist3 = n*(instMid3 - repmat( ls_est(1:2), num_inst, 1 ))';
    perpDist = (abs( perpDist1 ));% + abs( perpDist2 ) + abs( perpDist3 )) / 3;
    
    idx_PerpDist = abs( perpDist )' <= tau_perpdist;
    
    % 5. Central distance
    
    % 5-2. linelet
    dcentLl = instMid1 - repmat( ls_est(1:2), num_inst, 1 );
    dcentLl = sqrt( sum( dcentLl.^2, 2 ) );
    dbtLl = dcentLl - ll_inst(:,5)/2 - ls_est(4)/2;
    llCand5 = dbtLl <= min( repmat( ls_est(5), num_inst, 1 ), ll_inst(:, 5) ) & dbtLl >= 0;
    
    % wrt the target direction
    % Leftward linelet candidate
    idx_Left = (ll_inst(:,majorAxisIdx2) - ls_est(majorAxisIdx1)) <= 0;
%     llFinalCandL = find( llCand1 & llCand3LMajor & llCand3LMinor & llCand4 & llCand5 & llCand5L & ~ll_used );
    llFinalCandL = find( llCand3LMajor & llCand3LMinor & llCand5 & idx_PerpDist & idx_Left & ~ll_used );
    if ~isempty(llFinalCandL)
        % pick the closest one
        if length( llFinalCandL ) > 1
            [~, idxNear] = min(dbtLl(llFinalCandL));
            llFinalCandL = llFinalCandL(idxNear);
        end    
    end
    
    %% Rightward candidate
    idx_Right = (ll_inst(:,1) - ls_est(1)) >= 0;
%     llFinalCandR = find( llCand1 & llCand3RMajor & llCand3RMinor & idx_PerpDist & llCand5 & llCand5R & ~ll_used );
    llFinalCandR = find( llCand3RMajor & llCand3RMinor & llCand5 & idx_PerpDist & idx_Right & ~ll_used );
    if ~isempty(llFinalCandR)
        % pick the closest one
        if length( llFinalCandR ) > 1
            [~, idxNear] = min(dbtLl(llFinalCandR));
            llFinalCandR = llFinalCandR(idxNear);
        end    
    end
    
    %%
    if bDraw        
        % Original line segmnet
%         figure; imshow(im_gray); hold on;
        x1 = ls_est(idxSource,1:2) + ls_est(idxSource,4)/2*[cos(ls_est(idxSource,3)) sin(ls_est(idxSource,3))];
        x2 = ls_est(idxSource,1:2) - ls_est(idxSource,4)/2*[cos(ls_est(idxSource,3)) sin(ls_est(idxSource,3))];
        plot( [x1(1) x2(1)], [x1(2) x2(2)], 'color', [1 0 0], 'linewidth', 2 );
        % Merged line segment
        lsFinalCandL = find( lsCand1 & lsCand2 &  lsCand4 & lsCand5 & lsCand5L );
      
        idxTarget = lsFinalCandL;
        for k = 1:length(idxTarget)
            x1 = ls_est(idxTarget(k),1:2) + ls_est(idxTarget(k),4)/2*[cos(ls_est(idxTarget(k),3)) sin(ls_est(idxTarget(k),3))];
            x2 = ls_est(idxTarget(k),1:2) - ls_est(idxTarget(k),4)/2*[cos(ls_est(idxTarget(k),3)) sin(ls_est(idxTarget(k),3))];
            plot( [x1(1) x2(1)], [x1(2) x2(2)], 'color', [0 1 0], 'linewidth', 1 ); % updated
        end
        
        idxTarget = lsFinalCandR;
        for k = 1:length(idxTarget)
            x1 = ls_est(idxTarget(k),1:2) + ls_est(idxTarget(k),4)/2*[cos(ls_est(idxTarget(k),3)) sin(ls_est(idxTarget(k),3))];
            x2 = ls_est(idxTarget(k),1:2) - ls_est(idxTarget(k),4)/2*[cos(ls_est(idxTarget(k),3)) sin(ls_est(idxTarget(k),3))];
            plot( [x1(1) x2(1)], [x1(2) x2(2)], 'color', [0 1 0], 'linewidth', 2 ); % updated
        end
        title( ['Current index: ' num2str(idxSource)] );
    end
    
end

