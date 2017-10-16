function [ ll_idxSet, bAscend, im_ascend, im_descend, ll_usage_post, ll_Assigned, ll_Valid ] = ...
    LineletGrouping( ll_inst, upot, ll_type, level_t, im_gray, im_grad, im_dir, bDraw, param )
   
    tau_len = param.thres_tau_len; 
    tau_dist = 1;
        
    size_im = size(im_grad);
    num_ll_inst = size(ll_inst, 1);
    bAscend = [];
    
    im_ascend = zeros(size_im);
    im_descend = zeros(size_im);
    ll_Assigned = cell(num_ll_inst,1);
    ll_usage_post = [];
    ll_Valid = [];
    
    
    if strcmp(ll_type, 'hor')
        major_idx1 = 1;
        major_idx2 = 3;
        minor_idx1 = 2;
        minor_idx2 = 2;
    elseif strcmp(ll_type, 'ver')
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
    
    % ---------------------------------------------------
    % Ascending direction (1st octan)
    % ---------------------------------------------------
    LL_inst_new = [];
    ll_idxSet = {};
    num_inst_new = 1;
    ll_usage = zeros(num_ll_inst, 6); 
    idx_asds = 1;
    dist_minor = 1;
    for ii = 1:num_ll_inst
        try
            if ll_usage(ii, idx_asds) ~= 0 || ll_inst(ii, 6) == 0, continue; end
            
            stack_ll = ii;
            ll_cand = [];
            ll_cand_idx = [];
            
            
            cur_lsc_id = -1;
            while ~isempty(stack_ll)
                cur_ll = stack_ll(1);   stack_ll(1) = []; % pop
                
                if isempty(ll_cand)
                    base_len = ll_inst(cur_ll, 5);
                end
                    
                % Search possible connection with other linelets. At this step,
                % also check possible connection over one empty neighbor. This is
                % related to the minor axis update, but it is not certain yet
                
                % candidate region id
                idx_cr = ll_inst(cur_ll, 6) == ll_inst(:, 6);
                
                % Leftward distance
                idx_left = (ll_inst(cur_ll,major_idx1) - ll_inst(:,major_idx2) <= tau_dist) &...
                    (ll_inst(cur_ll,major_idx1) - ll_inst(:,major_idx2) >= 1) &...
                    (ll_inst(:,minor_idx1) - ll_inst(cur_ll,minor_idx2) == -dist_minor);
                
                % Rightward distance
                idx_right = (ll_inst(:,major_idx1) - ll_inst(cur_ll,major_idx2) <= tau_dist) &...
                    (ll_inst(:,major_idx1) - ll_inst(cur_ll,major_idx2) >= 1) &...
                    (ll_inst(:,minor_idx1) - ll_inst(cur_ll,minor_idx2) == dist_minor);
                
                % Length difference
                idx_len = abs(ll_inst(:,5) - ll_inst(cur_ll,5)) <= tau_len;                
                
                idx_nb = find((idx_right | idx_left) & idx_len & ll_usage(:,idx_asds) == 0 & idx_cr);
                idx_nb(idx_nb == cur_ll) = [];
                
                stack_ll = [stack_ll; idx_nb];
                
                ll_usage(idx_nb,idx_asds) = num_inst_new;
                ll_usage(cur_ll,idx_asds) = num_inst_new;
                
                ll_cand = [ll_cand; ll_inst(cur_ll,:)];
                ll_cand_idx = [ll_cand_idx; cur_ll];
                
                if isempty(stack_ll)
                    % ----------------------------
                    % Attach short linelet at both sides, this step is not considered for a moment
                    [~, idx1] = sort(ll_cand(:,major_idx1));
                    
                    % Leftward distance
                    idx_left = (ll_cand(idx1(1),major_idx1) - ll_inst(:,major_idx2) == 1) &...
                        (ll_inst(:,minor_idx1) - ll_cand(idx1(1),minor_idx2) == -dist_minor);
                    
                    % Rightward distance
                    idx_right = (ll_inst(:,major_idx1) - ll_cand(idx1(end),major_idx2) == 1) &...
                        (ll_inst(:,minor_idx1) - ll_cand(idx1(end),minor_idx2) == dist_minor);
                
                    idx_nb = find((idx_right | idx_left) & ll_usage(:,idx_asds) == 0  & idx_cr);
                    
                    ll_cand = [ll_cand; ll_inst(idx_nb,:)];
                    ll_cand_idx = [ll_cand_idx; idx_nb];
                    % ----------------------------
                end                
            end
            
            % ----------------------------
            % Line segment candidatie division
            % ----------------------------
            if length(ll_cand_idx) > 1
                [~, idx1] = sort(ll_cand(:,major_idx1));
                ll_cand = ll_cand(idx1,:);
                ll_cand_idx = ll_cand_idx(idx1);
                
                ptCent = [min(ll_cand(:,1)) + max(ll_cand(:,3)), min(ll_cand(:,2)) + max(ll_cand(:,4))]/2;
                angle = atan2(ll_cand(end, 4) - ll_cand(1, 2), ll_cand(end, 3) - ll_cand(1, 1));
                [ ptSet ] = Linelet2PtSet( ll_cand, ll_type );                
                
                perpDist = GetPerpDist(ptCent, ptSet, angle);
                [vMax, iMax] = max(perpDist);
                
                if vMax >= 2
                    idx = find((ll_cand(:,1) <= ptSet(iMax,1) & ll_cand(:,3) >= ptSet(iMax,1)) &...
                        (ll_cand(:,2) <= ptSet(iMax,2) & ll_cand(:,4) >= ptSet(iMax,2)));
                    ll_cand = ll_cand(1:idx, :);
                    ll_usage(ll_cand_idx(idx+1:end), idx_asds:idx_asds+2) = 0;
                    ll_cand_idx = ll_cand_idx(1:idx);
                end
            end
            % ----------------------------
            
            % ----------------------------
            % Validation step
            % ----------------------------
            % Measure candidate validity
            if sum(ll_cand(:, 5)) <= 4
                % this operation is only necessary for reducing computation. 
                ll_usage(ll_cand_idx, idx_asds:idx_asds+2) = 0;
                continue; 
            end
            
            tmpValid = LineletForegroundProbability( ll_cand, im_gray, im_grad, im_dir, ll_type );                        
            % ----------------------------            
            
            LL_inst_new{num_inst_new,1} = ll_cand;
            ll_idxSet{num_inst_new,1} = ll_cand_idx;
            ll_Valid(num_inst_new,:) = tmpValid;
            
            ll_usage(ll_cand_idx, idx_asds) = num_inst_new;
            ll_usage(ll_cand_idx, idx_asds+1) = sum(ll_cand(:,5));
            ll_usage(ll_cand_idx, idx_asds+2) = size(ll_cand, 1);
            
            for k = 1: length(ll_cand_idx)
                ll_Assigned{ll_cand_idx(k)} = unique([ll_Assigned{ll_cand_idx(k)}; num_inst_new]);
            end
            
            im_ascend = im_ascend + LineletDraw(ll_cand, ll_type, size_im) * num_inst_new;
            
            num_inst_new = num_inst_new + 1;            
        catch err
            fprintf('Error occured at %d of ascending.\n', ii);
            rethrow(err);
        end
    end
    if num_inst_new ~= 1
        bAscend(1:num_inst_new-1) = 1;
    end
    if bDraw
        figure; imagesc(im_ascend); axis off; axis image;
        title('Ascending sets');
    end
    % ---------------------------------------------------
    
    % ---------------------------------------------------
    % Descending direction (1st octan)
    % ---------------------------------------------------    
    dist_minor = -1;
    idx_asds=4;
    
    for ii = 1:num_ll_inst%\
        try
            if ll_usage(ii,idx_asds) ~= 0|| ll_inst(ii, 6) == 0, continue; end
            
            stack_ll = ii;
            ll_cand = [];
            ll_cand_idx = [];
            
            
            while ~isempty(stack_ll)
                cur_ll = stack_ll(1);   stack_ll(1) = []; % pop
                
                if isempty(ll_cand)
                    base_len = ll_inst(cur_ll, 5);
                end
                
                % Search possible connection with other linelets. At this step,
                % also check possible connection over one empty neighbor. This is
                % related to the minor axis update, but it is not certain yet
                
                % candidate region id
                idx_cr = ll_inst(cur_ll, 6) == ll_inst(:, 6);
                                
                % Leftward distance
                idx_left = (ll_inst(cur_ll,major_idx1) - ll_inst(:,major_idx2) <= tau_dist) &...
                    (ll_inst(cur_ll,major_idx1) - ll_inst(:,major_idx2) >= 1) &...
                    (ll_inst(:,minor_idx1) - ll_inst(cur_ll,minor_idx2) == -dist_minor);
                
                % Rightward distance
                idx_right = (ll_inst(:,major_idx1) - ll_inst(cur_ll,major_idx2) <= tau_dist) &...
                    (ll_inst(:,major_idx1) - ll_inst(cur_ll,major_idx2) >= 1) &...
                    (ll_inst(:,minor_idx1) - ll_inst(cur_ll,minor_idx2) == dist_minor);
                
                % Length difference
                idx_len = abs(ll_inst(:,5) - ll_inst(cur_ll,5)) <= tau_len;
                
                idx_nb = find((idx_right | idx_left) & idx_len & ll_usage(:,idx_asds) == 0 & idx_cr);
                idx_nb(idx_nb == cur_ll) = [];
                
                stack_ll = [stack_ll; idx_nb];
                
                ll_usage(idx_nb,idx_asds) = num_inst_new;
                ll_usage(cur_ll,idx_asds) = num_inst_new;
                
                ll_cand = [ll_cand; ll_inst(cur_ll,:)];
                ll_cand_idx = [ll_cand_idx; cur_ll];
                
                if isempty(stack_ll)
                    % ----------------------------
                    % Attach short linelet at both sides, this step is not considered for a moment
                    [~, idx1] = sort(ll_cand(:,major_idx1));
                    
                    % Leftward distance
                    idx_left = (ll_cand(idx1(1),major_idx1) - ll_inst(:,major_idx2) == 1) &...
                        (ll_inst(:,minor_idx1) - ll_cand(idx1(1),minor_idx2) == -dist_minor);
                    
                    % Rightward distance
                    idx_right = (ll_inst(:,major_idx1) - ll_cand(idx1(end),major_idx2) == 1) &...
                        (ll_inst(:,minor_idx1) - ll_cand(idx1(end),minor_idx2) == dist_minor);
                    
                    idx_nb = find((idx_right | idx_left) & ll_usage(:,idx_asds) == 0 & idx_cr);
                    
                    ll_cand = [ll_cand; ll_inst(idx_nb,:)];
                    ll_cand_idx = [ll_cand_idx; idx_nb];
                    % ----------------------------
                end
            end
           
            % ----------------------------
            % Line segment candidatie division
            % ----------------------------
            if length(ll_cand_idx) > 1
                [~, idx1] = sort(ll_cand(:,major_idx1));
                ll_cand = ll_cand(idx1,:);
                ll_cand_idx = ll_cand_idx(idx1);
                
                ptCent = [min(ll_cand(:,1)) + max(ll_cand(:,3)), min(ll_cand(:,2)) + max(ll_cand(:,4))]/2;
                angle = atan2(ll_cand(end, 4) - ll_cand(1, 2), ll_cand(end, 3) - ll_cand(1, 1));
                [ ptSet ] = Linelet2PtSet( ll_cand, ll_type );                
                
                perpDist = GetPerpDist(ptCent, ptSet, angle);
                [vMax, iMax] =max(perpDist);
                
                if vMax >= 2
                    idx = find((ll_cand(:,1) <= ptSet(iMax,1) & ll_cand(:,3) >= ptSet(iMax,1)) &...
                        (ll_cand(:,2) <= ptSet(iMax,2) & ll_cand(:,4) >= ptSet(iMax,2)));
                    ll_cand = ll_cand(1:idx, :);
                    ll_usage(ll_cand_idx(idx+1:end), idx_asds:idx_asds+2) = 0;
                    ll_cand_idx = ll_cand_idx(1:idx);
                end
            end
            % ----------------------------
            
            % ----------------------------
            % Validation step
            % ----------------------------
            % Measure candidate validity
            if sum(ll_cand(:, 5)) <= 4
                % this operation is only necessary for reducing computation. 
                ll_usage(ll_cand_idx, idx_asds:idx_asds+2) = 0;
                continue; 
            end
            
            tmpValid = LineletForegroundProbability( ll_cand, im_gray, im_grad, im_dir, ll_type );                        
            % ----------------------------            
            
            LL_inst_new{num_inst_new,1} = ll_cand;
            ll_idxSet{num_inst_new,1} = ll_cand_idx;
            ll_Valid(num_inst_new,:) = tmpValid;
            
            ll_usage(ll_cand_idx, idx_asds) = num_inst_new;
            ll_usage(ll_cand_idx, idx_asds+1) = sum(ll_cand(:,5));
            ll_usage(ll_cand_idx, idx_asds+2) = size(ll_cand, 1);
                        
            for k = 1: length(ll_cand_idx)
                ll_Assigned{ll_cand_idx(k)} = unique([ll_Assigned{ll_cand_idx(k)}; num_inst_new]);
            end
            
            im_descend = im_descend + LineletDraw(ll_cand, ll_type, size_im) * num_inst_new;
            
            num_inst_new = num_inst_new + 1;
        catch err
            fprintf('Error occured at %d of descending.\n', ii);
            rethrow(err);
        end
    end
    if bDraw
        figure; imagesc(im_descend); axis off; axis image;
        title('Descending sets');
    end
    
    if num_inst_new ~= 1
        bAscend(end+1:num_inst_new-1) = 0;
    else
        bAscend = 0;
    end    
    % ---------------------------------------------------
   
    
    
    % ---------------------------------------------------    
    % Prune redundant assignment    
    if level_t == 1
        ll_usage_post = ll_usage;
          im_ascend = zeros(size_im);
         im_descend = zeros(size_im);
   
        % Remove descending one that was assigned only once and whlie corresponding ascending one was more
        idxValid = ll_usage_post(:,1) ~= 0 & ll_usage_post(:,4) ~= 0;
        idx_comp = find(ll_usage_post(:,2) >= ll_usage_post(:,5) & ll_usage_post(:,3) > 1 & ll_usage_post(:,6) == 1 & idxValid);
        if ~isempty(idx_comp),  ll_usage_post(idx_comp, 4:end) = 0; end
        
        % Remove ascending one that was assigned only once and whlie corresponding ascending one was more
        idxValid = ll_usage_post(:,1) ~= 0 & ll_usage_post(:,4) ~= 0;
        idx_comp = find(ll_usage_post(:,2) <= ll_usage_post(:,5) & ll_usage_post(:,3) == 1 & ll_usage_post(:,6) > 1 & idxValid);
        if ~isempty(idx_comp),  ll_usage_post(idx_comp, 1:3) = 0; end
        
        % Remove either ascending one or descending one where both were assigned only once
        idxValid = ll_usage_post(:,1) ~= 0 & ll_usage_post(:,4) ~= 0;
        idx_comp = find(ll_usage_post(:,2) == ll_usage_post(:,5) & ll_usage_post(:,3) == 1 & ll_usage_post(:,6) == 1 &idxValid);
        if ~isempty(idx_comp),  ll_usage_post(idx_comp, 1:3) = 0; end
        
        % Any one having common linelet instance -- principle is drop the linelet within less relevant one -- use p(f(+)|O)
        idxValid = ll_usage_post(:,1) ~= 0 & ll_usage_post(:,4) ~= 0;
        idx_comp = find(ll_usage_post(:,3) > 1 & ll_usage_post(:,6) > 1 & idxValid);
        if ~isempty(idx_comp)
            p1 = ll_Valid(ll_usage_post(idx_comp, 1), 1);
            p2 = ll_Valid(ll_usage_post(idx_comp, 4), 1);
            idx = p1 >= p2;
            ll_usage_post(idx_comp(idx), 4:end) = 0;
            ll_usage_post(idx_comp(~idx), 1:3) = 0;
        end
           
        % Do re-assesment
        LL_set = []; num_set = 1;
        bAscend = [];
        
        valUnq1 = unique(ll_usage_post(:,1));
        for i = 1:length(valUnq1)
            if valUnq1(i) == 0, continue; end
            
            LL_set{num_set,1} = find(ll_usage_post(:,1) == valUnq1(i));
            num_set = num_set + 1;
            
            im_ascend = im_ascend + LineletDraw(ll_inst(LL_set{num_set-1,1},:), ll_type, size_im) * (num_set-1);
            
        end
        if num_set ~= 1
            num_set = num_set - 1;
            bAscend(1:num_set) = 1;
        end
        
        valUnq2 = unique(ll_usage_post(:,4));
        for i = 1:length(valUnq2)
            if valUnq2(i) == 0, continue; end
            
            LL_set{num_set,1} = find(ll_usage_post(:,4) == valUnq2(i));
            num_set = num_set + 1;
            
            im_descend = im_descend + LineletDraw(ll_inst(LL_set{num_set-1,1},:), ll_type, size_im) * (num_set-1);
        end
        if num_set ~= 1
        num_set = num_set - 1;
        bAscend(end+1:num_set) = 0;
        end
        ll_idxSet = LL_set;
    end
    
    if bDraw
        figure; imagesc(im_ascend); axis off; axis image;
        title('Ascending sets after clean up');
        figure; imagesc(im_descend); axis off; axis image;
        title('Descending sets after clean up');
    end   
end

