function [ ll_idxSetNew, pAngNew, lsEstNew, bAscendNew, bMerged, bLlMerged, ll_Assigned, ll_ValidNew, bUpdated ] = LineletAggregationHi(...
    ll_idxSet, pAng, ls_est, ll_inst, upot1, upot2,...
    bAscend, ll_type, level_t, theta_space, im_gray, im_grad, im_dir, bMerged, bLlMerged, ll_Assigned, ll_Valid, LineletNum, param, bDraw )
    %LINELETAGGREGATION Summary of this function goes here
   
    bUpdated = false;
        
    im_dird = im_dir * 180 / pi;
    size_im = size( im_grad );
    num_set = size( ll_idxSet, 1 );
    num_inst = size( ll_inst, 1 );
    bLsMerged = false( num_set, 1 );
    ll_Assigned = cell(num_inst,1);
    
    bLsVisited = zeros( num_set, 2 );
    ll_idxSetNew = [];
    pAngNew = [];
    lsEstNew = [];
    bAscendNew = [];
    ll_ValidNew =[];
    
    strMerge = '';
    
    if level_t == 1
        bMerged = [];
    end
    
    num_newSet = 1;
    idx_valid = true(size(ls_est,1),1);
    idx_merged_later = [];
    
    %%
    for ii = 1:num_set 
        try            
            if bLsVisited(ii,1) == 1
                continue;
            else
                bLsVisited(ii,1) = 1;
            end            
            
            if ~idx_valid(ii) || (level_t ~= 1 && ~bMerged(ii))
                ll_idxSetNew{num_newSet,1} = ll_idxSet{ii};
                pAngNew(num_newSet,:) = pAng(ii,:);
                lsEstNew(num_newSet,:) = ls_est(ii,:);
                bAscendNew(num_newSet,1) = bAscend(ii);
                bMerged(num_newSet,1) = false;
                ll_ValidNew(num_newSet,:) = ll_Valid(ii,:);
                bLsVisited(ii,2) = num_newSet;
                                
                % Record linelet assignments
                ll_cand_idx = ll_idxSetNew{num_newSet,1};
                for k = 1:length(ll_cand_idx)
                    ll_Assigned{ll_cand_idx(k)} = unique([ll_Assigned{ll_cand_idx(k)}; num_newSet]);
                end
                num_newSet = num_newSet + 1;
                continue; 
            end
            
            if bLsMerged(ii)
                continue; 
            end
            
            cur_ls = ls_est(ii,:);
            bUnmerged = true;
            bCandLs = false;
            bCandLL = false;
            [ lsCandL, lsCandR, llCandL, llCandR ] = GetMergeCandidate( cur_ls, ii, ls_est, ll_inst, ll_idxSet, ll_type, bAscend(ii), level_t, bDraw );
            
            lsCandSet = [];
            
            if ~isempty( lsCandL ) && ~bLsMerged(lsCandL)
                lsCandSet = [lsCandSet, lsCandL];
                bCandLs = true;
            end
            
            if ~isempty( lsCandR ) && ~bLsMerged(lsCandR)
                lsCandSet = [lsCandSet, lsCandR];
                bCandLs = true;
            end
            
            
            if bCandLs
                p_fin = zeros(length(lsCandSet), 2);
                p_f = zeros(length(lsCandSet), 2);
                
                p_fin = [];
                ll_Set = [];
                
                for k = 1:length(lsCandSet)
                    ll_NewCand = unique( [ll_idxSet{ii}; ll_idxSet{lsCandSet(k)}] );
                    ll_pseudo = GeneratePseudoLinelet( cur_ls, ls_est(lsCandSet(k),:), ll_inst(ll_idxSet{ii},:), ll_inst(ll_idxSet{lsCandSet(k)},:), ll_type );
                    curllSet = LineletRearrange(ll_inst(ll_NewCand,:), ll_type);
                    
                    [ curUpot1, curUpot2, ~, ~ ] = UnaryPot( curllSet, ll_type, LineletNum, im_gray, im_grad, im_dir, theta_space );
                    bAscendCur = CheckAscend( curllSet, ll_type );
                    
                    [pAng_tmp, ls_est_tmp, ~] = LineletInference( curllSet, curUpot1, curUpot2, ll_type, bAscendCur, theta_space );
                    p_fin(k,:) = FgBgProbability( ls_est_tmp, curllSet, im_gray, im_grad, im_dir, ll_type, param );
                    ll_Set{k,1} = curllSet;
                end
                
                [~, idxEst] = max(p_fin(:,1));
                
                if p_fin(idxEst,1) > p_fin(idxEst,2)
                    ll_NewCand = unique( [ll_idxSet{ii}; ll_idxSet{lsCandSet(idxEst)}] );
                    ls_NewCand = lsCandSet(idxEst);
                    ll_ValidNew(num_newSet,:) = p_fin(idxEst,:);
                    
                    curllSet = ll_Set{idxEst};
                    [ curUpot1, curUpot2, ~, ~ ] = UnaryPot( curllSet, ll_type, LineletNum, im_gray, im_grad, im_dir, theta_space );
                    bAscendCur = CheckAscend( curllSet, ll_type );
                    
                    [pAng_tmp, ls_est_tmp, ~] = LineletInference( curllSet, curUpot1, curUpot2, ll_type, bAscendCur, theta_space );
                    
                    % Calculate following properties: Kurtosis and fitness
                    if kurtosis( pAng_tmp ) >= kurtosis( pAng(ii,:) ) && ls_est_tmp(6) <= sqrt(2)
                        bUnmerged = false;
                    else
                        bUnmerged = true;
                    end
                end                
            end
               
            if bUnmerged == false
                ll_idxSetNew{num_newSet,1} = ll_NewCand;
                pAngNew(num_newSet,:) = pAng_tmp;
                lsEstNew(num_newSet,:) = ls_est_tmp;
                bAscendNew(num_newSet,1) = bAscendCur;
                bMerged(num_newSet,1) = true;
                
                bLsMerged(ii) = true;
                
                if ~bLsMerged(ls_NewCand,1)
                    idx_merged_later = [idx_merged_later; ls_NewCand];
                end
                
                bLsMerged(ls_NewCand) = true;
                bLlMerged(ll_NewCand) = true;
                bUpdated = true;
                bLsVisited(ii,2) = num_newSet;
                num_newSet = num_newSet + 1;
            else
                ll_idxSetNew{num_newSet,1} = ll_idxSet{ii};
                pAngNew(num_newSet,:) = pAng(ii,:);
                lsEstNew(num_newSet,:) = ls_est(ii,:);
                bAscendNew(num_newSet,1) = bAscend(ii);
                bMerged(num_newSet,1) = false;
                bLsVisited(ii, 2) = num_newSet;
                
                ll_ValidNew(num_newSet,:) = ll_Valid(ii,:);
                num_newSet = num_newSet + 1;
            end
            
            % Record linelet assignments
            ll_cand_idx = ll_idxSetNew{num_newSet-1,1};
            for k = 1:length(ll_cand_idx)                
                ll_Assigned{ll_cand_idx(k)} = unique([ll_Assigned{ll_cand_idx(k)}; num_newSet-1]);
            end
        catch err
            fprintf('Error at LineletAggregationHi of ii = %d, ll_type = %s.\n', ii, ll_type);
            rethrow(err);
        end
    end
    
    try
        if ~isempty(idx_merged_later)
            idx_remove = bLsVisited(idx_merged_later,2);
            idx_remove(idx_remove == 0) = [];
            ll_idxSetNew(idx_remove,:) = [];
            pAngNew(idx_remove,:) = [];
            lsEstNew(idx_remove,:) = [];
            bAscendNew(idx_remove,:) = [];
            bMerged(idx_remove,:) = [];
            ll_ValidNew(idx_remove,:) = [];
        end
    catch err
        fprintf('Error at LineletAggregationHi during dropping redundant ones --  ll_type = %s.\n', ii, ll_type);
        rethrow(err);
    end
end