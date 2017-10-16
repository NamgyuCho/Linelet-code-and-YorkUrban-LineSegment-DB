function [ret_dist, ls_est, curInstSet] = LineletInference( curllSet, curUpot1, curUpot2, ll_type, bAscend, theta_space0 )
    %LINELETINFERENCE Summary of this function goes here
    %   Detailed explanation goes here
    
    if strcmp(ll_type, 'ver')
        minorAxisIdx = 1;
        majorAxisIdx = 2;
    else
        minorAxisIdx = 2;
        majorAxisIdx = 1;
    end
    size_theta = length( theta_space0 );
    
    if bAscend
        if strcmp(ll_type, 'ver')            
            theta_space0 = pi/2 - theta_space0;
        elseif strcmp(ll_type, 'diag')
            theta_space0 = theta_space0(size_theta:-1:1);
        end
    else
        if strcmp(ll_type, 'ver')
            theta_space0 = theta_space0 - pi/2;
        elseif strcmp(ll_type, 'diag')
            theta_space0 = -theta_space0(size_theta:-1:1);
        else
            theta_space0 = -theta_space0;        
        end
    end
    
    %ret_dist = zeros(1, size_theta);
    ls_est = [];
    
    % Sort instance and potential wrt the left position of instance
    % Sort instance
    [~, srt_idx] = sort(curllSet(:, majorAxisIdx));
    curllSet = curllSet(srt_idx,:);
    
    % Sort potential
    curUpot1 = curUpot1(srt_idx,:);
    curUpot2 = curUpot2(srt_idx,:);
    
    % Sort index
    num_curInstSet0 = size(curllSet, 1);
    curInstSet0 = 1:num_curInstSet0;
    %curCandIdx0 = curCandIdx0(srt_idx);
    
    curPpot1 = 0;
    curPpot2 = zeros(1, size_theta);
    
    %------------------------------------------------------------------
    % Select major linelets
    %------------------------------------------------------------------
    switch num_curInstSet0
        case 1
            curInstSet = curInstSet0;
        case 2
            [valLong, idxLong] = max(curllSet(curInstSet0,5));
            idxValid = idxLong;
            curInstSet = curInstSet0(idxValid);
        otherwise
            curInstSet = curInstSet0(2:end-1);
    end
    num_CurInstSet = length(curInstSet);
    %------------------------------------------------------------------
    
    if num_CurInstSet == 1% Single linelet
        if num_curInstSet0 == 1
            curPpot1 = 0;
            curPpot2 = curPpot2 + curUpot2(curInstSet,:) + curUpot2(curInstSet,:);
            curPpot2(1) = max(curPpot2) + eps;
        else
            curPpot1 = 0;
            curPpot2 = curPpot2 + curUpot2(curInstSet,:);
        end
    else % More than two linelets
        for ij = 1:num_CurInstSet-1
            if abs(curllSet(curInstSet(ij), 5) - curllSet(curInstSet(ij+1), 5)) <= 3 &&... % length difference
                    abs(curllSet(curInstSet(ij), minorAxisIdx) - curllSet(curInstSet(ij+1), minorAxisIdx)) == 1 &&... % minor position difference
                    1 % slope difference
                curPpot1 = curPpot1 + 1;
                curPpot2 = curPpot2 + curUpot2(curInstSet(ij),:) + curUpot2(curInstSet(ij+1),:);
            end            
        end
    end
    
    % Score so far
    ret_dist = sum(curUpot2(curInstSet,:),1) + curPpot1;% + curPpot2;   %  figure; bar(ret_dist)
    ret_dist = ret_dist / sum( ret_dist );
    
    if (num_CurInstSet * num_curInstSet0) == 1 
        ret_dist = (size_theta:-1:1);
        ret_dist = ret_dist / sum( ret_dist );
    end
    
    % Global fitting
    [ ptSet ] = Linelet2PtSet( curllSet, ll_type );
    cent_est = [curllSet(curInstSet(1), 1) + curllSet(curInstSet(end), 3),...
        curllSet(curInstSet(1), 2) + curllSet(curInstSet(end), 4)]/2;
    
    [~, idxTopCand] = sort(ret_dist, 'descend');
    idxTopCand = 1:length( theta_space0 );%idxTopCand(1:300);
    perpDist = zeros( 1, length( idxTopCand ) );
    for i = 1:length( idxTopCand )
        tmp = GetPerpDist( cent_est, ptSet, theta_space0(idxTopCand(i)) );
        perpDist(i) = sum(tmp) / size(ptSet, 1);
        
    end
    fit_score = 1./(perpDist.^2);
    
    % Final score
    ret_dist = ret_dist + fit_score;
    ret_dist = ret_dist / sum(ret_dist);        
    
    [~, maxIdx] = max(ret_dist);
    
    cent_est = [curllSet(curInstSet(1), 1) + curllSet(curInstSet(end), 3),...
        curllSet(curInstSet(1), 2) + curllSet(curInstSet(end), 4)]/2;
    ang_est = theta_space0(maxIdx);
    
    len_est = sqrt( (max( [curllSet(:,1); curllSet(:,3)] ) - min( [curllSet(:,1); curllSet(:,3)] ))^2 +...
                    (max( [curllSet(:,2); curllSet(:,4)] ) - min( [curllSet(:,2); curllSet(:,4)]) )^2 );
    ls_est = [cent_est, ang_est, len_est, max(curllSet(:, 5)), perpDist(maxIdx) ];
    
    if bAscend
        len_ldiff = sqrt( ( min([curllSet(:,1); curllSet(:,3)] ) - cent_est(1) )^2 +...
            ( min([curllSet(:,2); curllSet(:,4)] ) - cent_est(2) )^2 );
        len_rdiff = sqrt( ( max([curllSet(:,1); curllSet(:,3)] ) - cent_est(1) )^2 +...
            ( max([curllSet(:,2); curllSet(:,4)] ) - cent_est(2) )^2 );
    else
        len_ldiff = sqrt( ( min([curllSet(:,1); curllSet(:,3)] ) - cent_est(1) )^2 +...
            ( max([curllSet(:,2); curllSet(:,4)] ) - cent_est(2) )^2 );
        len_rdiff = sqrt( ( max([curllSet(:,1); curllSet(:,3)] ) - cent_est(1) )^2 +...
            ( min([curllSet(:,2); curllSet(:,4)] ) - cent_est(2) )^2 );
    end
    adj_diff = len_ldiff - len_rdiff;
    ls_est(1:2) = ls_est(1:2) - adj_diff/2*[cos(ls_est(3)) sin(ls_est(3))];    
end

