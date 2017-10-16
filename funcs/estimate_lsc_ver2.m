function [line_segment, linelet_set] = estimate_lsc_ver2(ll_seeds, lkMap, gradMap, fgMap, im_dird, param, bHorMajor)

    global im_ll_map

    linelet_set = [];

    % gradMap = cur_lsc;
    bDebugDisplay = false;
    size_im = size(gradMap);
    num_seeds = size(ll_seeds, 1);
    line_segment = [];
    
    [rr,cc] = find(fgMap);
    pt_fg = [cc'; rr'];
    
    
    if strcmp(param.est_method, 'PROPOSED_CONSAC')
        % Get the initial angle and tracking offset
        %---------------------------------------------------------------------------------------
        % Estimate line segments via CONSAC
        pts = [ll_seeds(:,1)'; ll_seeds(:,2)'];% [xx'; yy'];
        
        num_sample = 2;
        num_iter = size(pts,2);
        num_pts = size(pts,2);
        
        thres_dist = .75;
        thres_inlier_ratio = .1;
        thres_inlier_num = round(thres_inlier_ratio*num_pts);
        
        num_inlier = zeros(1,num_iter);
        num_aligned = zeros(1, num_iter);
        
        idx_inlier = cell(1,num_iter);
        
        theta1 = zeros(1,num_iter);
        rho1 = zeros(1,num_iter);
        cts = zeros(num_iter, 2);
        lts = zeros(1, num_iter);
        
        for p = 1:num_iter
            % 1. fit using 2 random points
            sampleIdx = randIndex(num_pts,num_sample);
            ptSample = pts(:,sampleIdx);
            
            
            d = ptSample(:,2)-ptSample(:,1);
            d = d/norm(d); % direction vector of the line
            
            % 2. count the inliers, if more than thres_inlier_num, refit; else iterate
            n = [-d(2),d(1)]; % unit normal vector of the line
            dist1 = n*(pts-repmat(ptSample(:,1),1,num_pts));
            inlier1 = find(abs(dist1) < thres_dist);
            
            num_inlier(p) = length(inlier1);
            idx_inlier{p} = inlier1;
            
            tmp = lkMap(inlier1, inlier1);
            num_aligned(p) = sum(tmp(:)) / length(inlier1)^2;
            
            if length(inlier1) < thres_inlier_num, continue; end
            [ev, scr] = princomp(pts(:,inlier1)');
            
            d1 = ev(:,1);
            theta1(p) = -atan2(d1(2),d1(1)); % save the coefs
            rho1(p) = [-d1(2),d1(1)]*mean(pts(:,inlier1),2);
            cts(p,:) = mean(pts(:,inlier1)');
            lts(p) = sqrt(sum(d.^2));
        end
        
        % 3. choose the coef with the most inliers
        [~,idx] = max(num_inlier);
        
        %         if inlrDirVar(idx) > 10, continue; end
        theta = theta1(idx);
        %rho = rho1(idx);
        
        opt_pts = [ll_seeds(idx_inlier{idx}, 3) ll_seeds(idx_inlier{idx}, 5);
            ll_seeds(idx_inlier{idx}, 4) ll_seeds(idx_inlier{idx}, 6)];
        
        line_len = sqrt( sum((max(opt_pts) - min(opt_pts)).^2) );
        center = (max(opt_pts) + min(opt_pts))/2;
        
        angle = -theta;
        
        line_segment = [line_segment; center angle line_len];
    end
    
    if strcmp(param.est_method, 'PROPOSED_CONSAC1')
        ll_usage = false(num_seeds,1);
        
        [~, idx_srt] = sort(ll_seeds(:,4), 'descend');
        ll_seeds = ll_seeds(idx_srt,:); %or "LL_rep";
        idx = ll_seeds(:,5) < 0;   ll_seeds(idx,5) = ll_seeds(idx,5) + pi;
       
        for ii = 1:num_seeds
            if ll_usage(ii), continue; end
            
            stack_ll = ii;
            ll_cand = [];
            
            len_thres_major = max(max(ll_seeds(:,4))/2, 4); %4
            len_thres_minor = 1;
            len_thres_diff = 3;
            
            while ~isempty(stack_ll)
                cur_ll = stack_ll(1);   stack_ll(1) = []; % pop
                
                % Search possible connection with other linelets. At this step, also check possible connection over one
                % empty neighbor. This is related to the minor axis update, but it is not certain yet
                idx_left = (abs(ll_seeds(cur_ll,1) - ll_seeds(:,2)) <= len_thres_major) & (abs(ll_seeds(:,3) - ll_seeds(cur_ll,3)) <= len_thres_minor); % Leftward distance
                idx_right = (abs(ll_seeds(cur_ll,2) - ll_seeds(:,1)) <= len_thres_major) & (abs(ll_seeds(:,3) - ll_seeds(cur_ll,3)) <= len_thres_minor); % Rightward distance
                idx_ang = bAngleAligned( ll_seeds(cur_ll,5), ll_seeds(:,5), param.thres_angle_diff ); % Gradient angle            
                idx_len = abs(ll_seeds(:,4) - ll_seeds(cur_ll,4)) <= len_thres_diff;
                
                %idx_nb = find((idx_left | idx_right) & idx_len & ~ll_usage);
                %idx_nb = find((idx_left | idx_right) & idx_ang & ~ll_usage);
                idx_nb = find((idx_left | idx_right) & ~ll_usage);
                idx_nb(idx_nb == cur_ll) = [];
                
                stack_ll = [stack_ll; idx_nb];
                
                ll_usage(idx_nb) = true;
                ll_usage(cur_ll) = true;
                
                ll_cand = [ll_cand; ll_seeds(cur_ll,:)];
            end
%             figure; imagesc(fgMap'); hold on; plot(ll_cand(:,1), ll_cand(:,3), 'go'); hold on
            
            
            if isempty(ll_cand), continue; end
            %         im_tmp = DrawLL(ll_cand, size_im); fig91 = figure(91); imshowpair(lbMap_LL_Set1>0, im_tmp>0, 'blend'); ii
            
            if size(ll_cand,1) == 1
                line_segment = [line_segment; ll_cand(1) + ll_cand(4)/2 ll_cand(3) 0 ll_cand(4)];
            else
                pts = [ll_cand(:,1)' + ll_cand(:,4)'/2; ll_cand(:,3)'];% [xx'; yy'];
                pts_weight = ll_cand(:,4) ./ sum(ll_cand(:,4));
                        
                num_sample = 2;
                num_iter = size(pts,2);
                num_pts = size(pts,2);
                
                thres_dist = .5;
                thres_inlier_ratio = .1;
                thres_inlier_num = round(thres_inlier_ratio*num_pts);
                
                % scores to optimize
                num_inlier = zeros(1,num_iter);
                num_aligned = zeros(1, num_iter);                
                std_linelet_len = zeros(1, num_iter);
                score_linelet = zeros(1, num_iter);
                
                idx_inlier = cell(1,num_iter);
                
                theta1 = zeros(1,num_iter);
                rho1 = zeros(1,num_iter);
                cts = zeros(num_iter, 2);
                lts = zeros(1, num_iter);
                
                for p = 1:num_iter
                    % 1. fit using 2 random points
                    sampleIdx = randIndex(num_pts,num_sample);
                    ptSample = pts(:,sampleIdx);
                    
                    
                    d = ptSample(:,2)-ptSample(:,1);
                    d = d/norm(d); % direction vector of the line
                    
                    % 2. count the inliers, if more than thres_inlier_num, refit; else iterate
                    n = [-d(2),d(1)]; % unit normal vector of the line
                    dist1 = n*(pts-repmat(ptSample(:,1),1,num_pts));
                    inlier1 = find(abs(dist1) < thres_dist);
                    
                    num_inlier(p) = length(inlier1);
                    idx_inlier{p} = inlier1;
                    
                    tmp = lkMap(inlier1, inlier1);
                    num_aligned(p) = sum(tmp(:)) / length(inlier1)^2;
                    std_linelet_len(p) = 1/std(ll_cand(inlier1,4)); % std of linelet lengths
                    score_linelet(p) = sum(pts_weight(inlier1));
                    
                    if length(inlier1) < thres_inlier_num, continue; end
                    [ev, scr] = princomp(pts(:,inlier1)');
                    
                    d1 = ev(:,1);
                    theta1(p) = -atan2(d1(2),d1(1)); % save the coefs
                    rho1(p) = [-d1(2),d1(1)]*mean(pts(:,inlier1),2);
                    cts(p,:) = mean(pts(:,inlier1)');
                    lts(p) = sqrt(sum(d.^2));
                end
                
                % 3. choose the coef with the most inliers
                %[~,idx] = max(num_inlier + num_aligned + std_linelet_len + score_linelet);
                [~,idx] = max(num_inlier + num_aligned + score_linelet);
                
                %         if inlrDirVar(idx) > 10, continue; end
                theta = theta1(idx);
                %rho = rho1(idx);
                
                opt_pts = [ll_cand(idx_inlier{idx}, 1) ll_cand(idx_inlier{idx}, 3);
                    ll_cand(idx_inlier{idx}, 2) ll_cand(idx_inlier{idx}, 3)];
%                 figure; imagesc(fgMap); hold on; plot(opt_pts(:,1), opt_pts(:,2), 'ro');
            
                
                line_len = sqrt( sum((max(opt_pts) - min(opt_pts)).^2) );
                center = (max(opt_pts) + min(opt_pts))/2;
                
                angle = -theta;
                %angle = atan2( opt_pts(end,2) - opt_pts(1,2), opt_pts(end,1) - opt_pts(1,1) );
                
                ls = [center angle line_len];
                line_segment = [line_segment; ls];
                
%                  for k=1:size(ls,1)
%                     x1 = ls(k,1:2) + ls(k,4)/2*[cos(ls(k,3)) sin(ls(k,3))];
%                     x2 = ls(k,1:2) - ls(k,4)/2*[cos(ls(k,3)) sin(ls(k,3))];
%                     plot([x1(1) x2(1)], [x1(2) x2(2)], 'r-', 'linewidth', 1); % updated
%                 end
            end
            
            %         [ ll_cand_new, conn_map_new ] = CompleteLineletSet2(ll_cand, param); % strictly, "CompleteLineletSet2" doesn't act exactly as expected as CompleteSet scheme in the paper
            %
            %         % estimate
            %         for i = 1:size(ll_cand_new, 1)
            %             ll_cand = ll_cand_new{i};
            %
            %             %[BN_cand, Pot] = LineletPot(gradMap, ll_cand, BN_Model);
            %             %[~, iMax] = max(Pot);
            %             %ang = theta_BN(BN_cand(iMax));                                               % Estimate by the energy model
            %             ang = atan2( ll_cand(end,3) - ll_cand(1,3), ll_cand(end,2) - ll_cand(1,1) ); % Estimate by rectangular approaximation
            %
            %             if (ang > pi/4 && ang < pi*3/4) || (ang < -pi/4 && ang > -pi*3/4), continue; end
            %
            %             len = sqrt( (ll_cand(end,2)-ll_cand(1,1))^2 + (ll_cand(end,3) - ll_cand(1,3))^2 );
            %
            %             %%%%%%%%%%%%%%%%%%%%%%%%%%% THIS PART SHOULD BE MODIFIED %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             if len <= param.thres_min_linelet_len, continue; end  %if len <= param.thres_linelet_len, continue; end
            %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %             line_segment = [line_segment; (ll_cand(1,1) + ll_cand(end,2))/2, (ll_cand(1,3) + ll_cand(end,3))/2, ang, len];
            %
            %             if bDebugDisplay
            %                 x1 = line_segment(end,1:2) + line_segment(end,4)/2*[cos(line_segment(end,3)) sin(line_segment(end,3))];
            %                 x2 = line_segment(end,1:2) - line_segment(end,4)/2*[cos(line_segment(end,3)) sin(line_segment(end,3))];
            %                 plot([x1(1) x2(1)], [x1(2) x2(2)], 'r-', 'linewidth', 2); % updated
            %                 waitforbuttonpress;
            %             end
            %         end
        end
    else
        
        
        %         % Generate seed points
        %         BW = im_lsc_track > 0; % im_lsc or im_lsc_track
        %         BWL = bwlabel(im_lsc, 8);
        %         CC = bwconncomp(BW, 8);
        %         im_tmp = zeros(size_im);
        %         for j = 1:CC.NumObjects
        %             ind = CC.PixelIdxList{j};
        %             if length(ind) < 2, continue; end
        %
        %             cur_grad = im_lsc(ind);
        %             [~, max_ind] = max(cur_grad);
        %             [yy, xx] = ind2sub(size_im, ind);
        %             seed_pts = [xx(max_ind) yy(max_ind); round(median([xx yy]))];
        %
        %             [coef, scr] = princomp([xx yy]);
        %             init_angle = atan2(coef(2,1), coef(1,1));
        %             init_angle_op = get_opposite_angle( init_angle, 'posneg');
        %             [offsets, sub_offset, bUseGX] = get_octan_rot_angle_and_offset(init_angle);
        %             [offsets_op, sub_offset_op, bUseGX_op] = get_octan_rot_angle_and_offset(init_angle_op);
        %
        %
        %             line_score = zeros(size(seed_pts,1),1);
        %             pixel_set = cell(size(seed_pts,1),1);
        %             end_points = zeros(size(seed_pts,1)*2,2);
        %             % Track from each seed point
        %             for i = 1:size(seed_pts,1)
        %                 pixel_visited = false(size_im);
        %
        %                 % Track two directions
        %                 [pixel_visited, end_points(i*2-1,:)] = track_local_maximum(im_lsc_track, pixel_visited, size_im, seed_pts(i,:), offsets);
        %                 [pixel_visited, end_points(i*2,:)] = track_local_maximum(im_lsc_track, pixel_visited, size_im, seed_pts(i,:), offsets_op);
        %
        %                 vind = find(pixel_visited ~= 0);
        %                 [rr, cc] = ind2sub(size_im, vind);
        %                 pts = [cc rr];
        %                 num_tracked = size(pts, 1);
        %
        %                 % Prune invalid pixel
        %                 [pts, ~] = return_valid_loc( pts, size_im, sub_offset );
        %                 sub_pts = pts + repmat(sub_offset, size(pts,1), 1);
        %
        %                 if size(pts,1) == 1, continue; end
        %
        %                 ind_main = sub2ind(size_im, pts(:,2), pts(:,1));
        %                 ind_sub  = sub2ind(size_im, sub_pts(:,2), sub_pts(:,1));
        %
        %                 grad_main = gradMap(ind_main);
        %                 grad_sub  = gradMap(ind_sub);
        %
        %                 % Calculate each lsc score (= energy)
        %                 [pro_x, pro_y] = bresenham(end_points(i*2-1, 1), end_points(i*2-1, 2), end_points(i*2, 1), end_points(i*2, 2));
        %                 num_model = size(pro_x,1);
        %
        %                 intersec_idx = intersect([pro_x, pro_y], pts, 'rows');
        %                 num_inter = size(intersec_idx, 1);
        %
        %                 if bUseGX
        %                     hist_pt = [pts(:,2) pts(:,1)];
        %                 else
        %                     hist_pt = pts;
        %                 end
        %
        %                 [unq_val, num_unq] = count_unique(pts(:,2));
        %                 lsc_hist = zeros(1, length(unq_val));
        %                 for k = 1:length(unq_val)
        %                     valid_idx = find(pts(:,2) == unq_val(k));
        %                     valid_idx = sub2ind(size_im, pts(valid_idx,2), pts(valid_idx,1));
        %                     lsc_hist(k) = sum(im_lsc_track(valid_idx)) / length(valid_idx);
        %                 end
        %
        %                 %         num_inter / num_model
        %                 %         figure(555);
        %                 %         tmp = zeros(size_im);
        %                 % %         tmp(sind) = 1;
        %                 %         tmp(ind_main) = 2;
        %                 %         tmp(sub2ind(size_im, pro_y, pro_x)) = 3;
        %                 %         imagesc(tmp);%imshowpair(im_gray, pixel_visited, 'blend');
        %                 %         waitforbuttonpress;
        %
        %                 %             if (num_inter / num_model) < .5
        %                 %                 line_score(i) = 0;
        %                 % %                 0
        %                 %             else
        %                 line_score(i) = num_inter^2 / (num_tracked * num_model) + 1 / std(grad_main + grad_sub);
        %                 %             end
        %                 pixel_set{i} = pts;
        %             end
        %             % Measure line candidate score
        %             [opt_val, opt_ind] = max(line_score);
        %             pts = pixel_set{opt_ind};
        %             pts_ext = end_points(opt_ind*2-1:opt_ind*2,:);
        %             line_len = pts_ext(1,:) - pts_ext(end,:);
        %             ang_val = atan2(line_len(2), line_len(1));
        %             line_len = sqrt(sum(line_len.^2));
        %
        %             if line_len < param.thres_min_ls_length || opt_val == 0 || opt_val == Inf || size(pts,1) == 1 || isempty(pts), continue; end
        %
        %             %         [yy, xx] = ind2sub(
        %             %         pts = [xx yy];
        %             a = tan(ang_val);
        %             b = -1;
        %             c = -a*pts(1,1) + pts(1,2);
        %             dist1 = abs(a * pts(:,1) + b * pts(:,2) + c) / sqrt(a^2 + b^2);
        %
        %             [val, ind] = max(dist1);
        %
        %
        %             if val >= 1.5
        %                 % split if the maximum distance exceeds the threshold
        %                 pts1 = [pts(1,:); pts(ind,:)];
        %                 line_len1 = pts1(1,:) - pts1(end,:);
        %                 ang_val1 = atan2(line_len1(2), line_len1(1));
        %                 line_len1 = sqrt(sum(line_len1.^2));
        %
        %                 line_segment(num_lsc,:) = [mean(pts1), line_len1, ang_val1, line_score(opt_ind)];
        %                 pixel_cluster{num_lsc,1} = sub2ind(size_im, pts(1:ind,2), pts(1:ind,1));
        %                 im_tmp(pixel_cluster{num_lsc,1}) = im_tmp(pixel_cluster{num_lsc,1}) + 1;
        %                 num_lsc = num_lsc+1;
        %
        %                 pts1 = [pts(ind,:); pts(end,:)];
        %                 line_len1 = pts1(1,:) - pts1(end,:);
        %                 ang_val1 = atan2(line_len1(2), line_len1(1));
        %                 line_len1 = sqrt(sum(line_len1.^2));
        %
        %                 line_segment(num_lsc,:) = [mean(pts1), line_len1, ang_val1, line_score(opt_ind)];
        %                 pixel_cluster{num_lsc,1} = sub2ind(size_im, pts(ind:end,2), pts(ind:end,1));
        %                 im_tmp(pixel_cluster{num_lsc,1}) = im_tmp(pixel_cluster{num_lsc,1}) + 1;
        %                 num_lsc = num_lsc+1;
        %             else
        %
        %                 line_segment(num_lsc,:) = [mean(pts), line_len, ang_val, line_score(opt_ind)];
        %                 pixel_cluster{num_lsc,1} = sub2ind(size_im, pts(:,2), pts(:,1));
        %                 im_tmp(pixel_cluster{num_lsc,1}) = im_tmp(pixel_cluster{num_lsc,1}) + 1;
        %                 num_lsc = num_lsc+1;
        %             end
        %         end
        
    end
    
    %%
    if strcmp(param.est_method, 'PROPOSED_CONSAC2')
        ll_usage = false(num_seeds,1);
        
        [~, idx_srt] = sort(ll_seeds(:,1));
        ll_seeds = ll_seeds(idx_srt,:); %or "LL_rep";
        idx = ll_seeds(:,5) < 0;   ll_seeds(idx,5) = ll_seeds(idx,5) + pi;
       
        n_est = 1;
        for ii = 1:num_seeds
            %%
            if ll_usage(ii), continue; end
            
            stack_ll = ii;
            ll_cand = [];
            
            len_thres_major = max(max(ll_seeds(:,4)), 4); %4
            len_thres_minor = 2;
            len_thres_diff = 3;
            
            while ~isempty(stack_ll)
                cur_ll = stack_ll(1);   stack_ll(1) = []; % pop
                
%                 len_thres_major = max(len_thres_major * 0.3 + ll_seeds(cur_ll,4) * 0.7, 1); %4
%                 len_thres_minor = 1;
            
                
                % Search possible connection with other linelets. At this step, also check possible connection over one
                % empty neighbor. This is related to the minor axis update, but it is not certain yet
                idx_left = (abs(ll_seeds(cur_ll,1) - ll_seeds(:,2)) <= len_thres_major) & (abs(ll_seeds(:,3) - ll_seeds(cur_ll,3)) <= len_thres_minor); % Leftward distance
                idx_right = (ll_seeds(:,1) - ll_seeds(cur_ll,2) <= len_thres_major) &...
                            (ll_seeds(:,1) - ll_seeds(cur_ll,2) > 0) &...
                            (abs(ll_seeds(:,3) - ll_seeds(cur_ll,3)) <= len_thres_minor); % Rightward distance
                idx_ang = bAngleAligned( ll_seeds(cur_ll,5), ll_seeds(:,5), param.thres_angle_diff ); % Gradient angle            
                idx_len = abs(ll_seeds(:,4) - ll_seeds(cur_ll,4)) <= len_thres_diff;
                
                %idx_nb = find((idx_left | idx_right) & idx_len & ~ll_usage);
                %idx_nb = find((idx_left | idx_right) & idx_ang & ~ll_usage);
                idx_nb = find((idx_right) & ~ll_usage);
                idx_nb(idx_nb == cur_ll) = [];
                
                % pick the closest one
                if length(idx_nb) > 1
                    [~, idx] = min(ll_seeds(idx_nb,1) - ll_seeds(cur_ll,2) + abs(ll_seeds(idx_nb,3) - ll_seeds(cur_ll,3)));
                    idx_nb = idx_nb(idx);
                end
                
                stack_ll = [stack_ll; idx_nb];
                
                ll_usage(idx_nb) = true;
                ll_usage(cur_ll) = true;
                
                ll_cand = [ll_cand; ll_seeds(cur_ll,:)];
            end
%             figure; imagesc(fgMap); hold on; plot(ll_cand(:,1), ll_cand(:,3), 'go'); hold on
            
            %%
            if isempty(ll_cand), continue; end
            %         im_tmp = DrawLL(ll_cand, size_im); fig91 = figure(91); imshowpair(lbMap_LL_Set1>0, im_tmp>0, 'blend'); ii
            
            %
            % i) Build the connection map and conistency map
            %
            % Build linelet connection map: non-zero value of conn_map(i,j) means that linelet i and j are connected such that
            %                               positive/negative for i is below/above j, and 1/2 for i is left/right to j
            [ ll_cand_new, conn_map_new ] = CompleteLineletSet3(ll_cand, param); % strictly, "CompleteLineletSet2" doesn't act exactly as expected as CompleteSet scheme in the paper
            
            % estimate
            for i = 1:size(ll_cand_new, 1)
                ll_cand = ll_cand_new{i};
                
%                 if bHorMajor
%                     im_ll_map = im_ll_map + DrawLL(ll_cand, size_im);
%                 else
%                     im_ll_map = im_ll_map + DrawLL(ll_cand, [size_im(2) size_im(1)])';
%                 end

%                 % validation
%                 if NFA_linelet(ll_cand, im_dird) <= param.thres_log_eps, continue; end
                
                %if ~HelmholtzValidation(ll_cand, size_im), continue; end
                
                ls = estimate_lineseg(ll_cand, param.est_line_seg);
                
                if (ls(3) > pi/4 && ls(3) < pi*3/4) || (ls(3) < -pi/4 && ls(3) > -pi*3/4), continue; end
                
%                 [BN_cand, Pot] = LineletPot(gradMap, ll_cand, BN_Model);
%                 [~, iMax] = max(Pot);
%                 ang = theta_BN(BN_cand(iMax));                                               % Estimate by the energy model
%                 ang = atan2( ll_cand(end,3) - ll_cand(1,3), ll_cand(end,2) - ll_cand(1,1) ); % Estimate by rectangular approaximation
                
                
                
%                 len = sqrt( (ll_cand(end,2)-ll_cand(1,1))^2 + (ll_cand(end,3) - ll_cand(1,3))^2 );
%                 
                line_segment = [line_segment; ls];
                linelet_set{n_est,1} = ll_cand;
                
                n_est = n_est + 1;
                
                if bDebugDisplay
                    x1 = line_segment(end,1:2) + line_segment(end,4)/2*[cos(line_segment(end,3)) sin(line_segment(end,3))];
                    x2 = line_segment(end,1:2) - line_segment(end,4)/2*[cos(line_segment(end,3)) sin(line_segment(end,3))];
                    plot([x1(1) x2(1)], [x1(2) x2(2)], 'r-', 'linewidth', 2); % updated
%                     waitforbuttonpress;
                end
            end            
        end    
    end