function [ ll_idxSet, ll_Assigned, lsEst, pAng, bAscend ] = LineletCleanUp( ll_idxSet, ll_Assigned, ll_inst, lsEst, pAng, bAscend, bDraw )
    %LINELETCLEANUP Summary of this function goes here
    %   Detailed explanation goes here
    
    % ll_idxSet = ll_idxSet_Hor; ll_Assigned = ll_Assigned_Hor; ll_inst = ll_Hor; lsEst = ls_est_Hor; pAng = pAng_Hor; bAscend = bAscend_Hor; bDraw = false;
    % ll_idxSet = ll_idxSet_Ver; ll_Assigned = ll_Assigned_Ver; ll_inst = ll_Ver; lsEst = ls_est_Ver; pAng = pAng_Ver; bAscend = bAscend_Ver; bDraw = false;
    
    % ------------------------------------------------------------
    % Find line segments having redundant linelet
    tmp = cellfun(@length, ll_Assigned);
    idxTarget = find(tmp > 1);
    lsTarget = cell(length(idxTarget),1);
    [lsTarget{:}] = ll_Assigned{idxTarget,1};
    
    % remove redundant instances -- simply drop shorter one?
    idxDrop = [];
    for ii = 1:length(idxTarget)
        % Usually, a linelet can be assigned twice at most, but should be careful for outliers
        %for k = 1:length(ll_Assigned{idxTarget(ii)}),         end -- use this loop for the outliers
        idx1 = ll_Assigned{idxTarget(ii)}(1);
        idx2 = ll_Assigned{idxTarget(ii)}(2);
        
        ls1 = lsEst(idx1,:);
        ls2 = lsEst(idx2,:);
        
        % This is mostly a case where ascending and descending line segments conflict each other
        if length(ll_idxSet{idx1}) <= 2 || length(ll_idxSet{idx2}) <= 2
            if ls1(4) < ls2(4)
                %                 if length(ll_idxSet{idx1}) < 3
                % Drop instance
                idxDrop = [idxDrop; idx1];
                %                 else
                % Drop linelet only
                %                 end
            else
                idxDrop = [idxDrop; idx2];
            end
        else
            continue;
        end
        
        
        
        if bAngleAligned( ls1(3), ls2(3), pi/16 )
%             idxNew = unique( [ll_idxSet{idx1}; ll_idxSet{idx2}] );
%             
%             curllSet = ll_inst(idxNew,:);
%             curUpot1 = upot1(idxNew, :);
%             curUpot2 = upot2(idxNew, :);
%             
%             % Calculate following properties: Kurtosis and fitness
%             [pAng_tmp, ls_est_tmp, ~] = LineletInference( curllSet, curUpot1, curUpot2, ll_type, bAscend(idx1), theta_space );
%             
%             if kurtosis( pAng_tmp ) >= kurtosis( pAngNew(idx1,:) ) && ls_est_tmp(6) >= lsEst(ii, 6)
%                 ll_idxSet{idx1,1} = idxNew;
%                 pAngNew(idx1,:) = pAng_tmp;
%                 lsEst(idx1,:) = ls_est_tmp;
%                 bMerged(idx1,1) = true;
%                 bMerged(idx2,1) = true;
%                 
%                 bLsMerged(idx1) = true;
%                 bLsMerged(idx2) = true;
%                 idxDrop = [idxDrop; idx2];
%             else
%                 if ls1(4) < ls2(4), idxDrop = [idxDrop; idx1];
%                 else                idxDrop = [idxDrop; idx2];
%                 end
%             end
%         else
%             if ls1(4) < ls2(4)
%                 %                 if length(ll_idxSet{idx1}) < 3
%                 % Drop instance
%                 idxDrop = [idxDrop; idx1];
%                 %                 else
%                 % Drop linelet only
%                 %                 end
%             else
%                 idxDrop = [idxDrop; idx2];
%             end
        end
        
        if bDraw
            figure(13); imagesc(im_gray); axis image; hold on; colormap gray
            % Original line segmnet
            x1 = ls1(1:2) + ls1(4)/2*[cos(ls1(3)) sin(ls1(3))];
            x2 = ls1(1:2) - ls1(4)/2*[cos(ls1(3)) sin(ls1(3))];
            plot( [x1(1) x2(1)], [x1(2) x2(2)], 'g-', 'linewidth', 1 );
            x1 = ls2(1:2) + ls2(4)/2*[cos(ls2(3)) sin(ls2(3))];
            x2 = ls2(1:2) - ls2(4)/2*[cos(ls2(3)) sin(ls2(3))];
            plot( [x1(1) x2(1)], [x1(2) x2(2)], 'r-.', 'linewidth', 1 ); % updated
        end
    end
    
    % Clean up unwanted instances
    pAng(idxDrop,:) = [];
    lsEst(idxDrop,:) = [];
    bAscend(idxDrop) = [];
    idxSetValid = true(length(ll_idxSet),1);
    idxSetValid(idxDrop) = false;
    ll_idxSetTmp = cell(length(find(idxSetValid)),1);
    [ll_idxSetTmp{:}] = ll_idxSet{idxSetValid};
    ll_idxSet = ll_idxSetTmp;
    % ------------------------------------------------------------
    
    if bDraw
        figure; imagesc(im_gray); axis image; axis off; hold on; colormap gray
        line_own = lsEst;  lcol = [1 0 0];
        for k = 1:size(line_own,1)
            %         if line_own(k,4) <= 4, continue; end
            x1 = line_own(k,1:2) + line_own(k,4)/2*[cos(line_own(k,3)) sin(line_own(k,3))];
            x2 = line_own(k,1:2) - line_own(k,4)/2*[cos(line_own(k,3)) sin(line_own(k,3))];
            plot( [x1(1) x2(1)], [x1(2) x2(2)], '-', 'linewidth', 1, 'color', lcol ); % updated
            %         text(line_own(k,1), line_own(k,2)-1, num2str(k), 'color', [1 0 1], 'fontsize', 12);
            n0 = n0+1;
        end
    end
end

