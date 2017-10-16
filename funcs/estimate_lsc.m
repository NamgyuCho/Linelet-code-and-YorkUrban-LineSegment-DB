function [line_segment, pixel_cluster] = estimate_lsc(lc, weight, im_grad, im_gx, im_gy, im_dird, im_nms, param)

size_im = size(im_nms);
center = [];
angle = [];
lengths = [];
pixel_cluster = {};
line_segment = [];

if strcmp(param.est_method, 'PCA') % PCA
    % Get main pc's angle
    [coef, scr] = princomp(lc);            
    pc1_angle = atan2(coef(2,1), coef(1,1));
    pc2_angle = atan2(coef(2,2), coef(1,2));
    lengths = max(scr);
    
    % Get histogram along minor axis
    bin_center_size = 1;
    bin_width = ceil(abs([lengths(1), lengths(2)]));
    
    gradient_bin = zeros(1, sum(bin_width)+1);
    gradient_bin_loc = zeros(1, sum(bin_width)+1);
    gradient_bin_cnt = zeros(1, sum(bin_width)+1);
    acc_grad = 0;
    
    nms_bin = zeros(1, sum(bin_width)+1);
    acc_nms = 0;
    for k = 1:size(lc, 1)
        bin_pos = ceil(scr(k,2) + bin_width(1))+1;
        
        gradient_bin(bin_pos) = gradient_bin(bin_pos) + im_grad(lc(k,2), lc(k,1));
        gradient_bin_loc(bin_pos) = gradient_bin_loc(bin_pos) + scr(k,2) ;
        gradient_bin_cnt(bin_pos) = gradient_bin_cnt(bin_pos) + 1;
        acc_grad = acc_grad + im_grad(lc(k,2), lc(k,1));
        
        nms_bin(bin_pos) = nms_bin(bin_pos) + im_nms(lc(k,2), lc(k,1));
        acc_nms = acc_nms + im_nms(lc(k,2), lc(k,1));
    end
%     lc_grad_hist{usage_count,1} = scr(:,2);
    %     gradient_bin = smooth(gradient_bin);
    
    % Find local peaks
    if length(gradient_bin) >= 3
        
        [pks,locs] = findpeaks(double(gradient_bin) / double(acc_grad), 'MINPEAKHEIGHT', .1, 'NPEAKS', 2);
        [pks1,locs1] = findpeaks(double(nms_bin) / double(acc_nms), 'MINPEAKHEIGHT', .1, 'NPEAKS', 2);
        if length(locs) > 1
            cn = [];
            for k = 1:length(locs)
                %                 cn = [cn; center + gradient_bin_loc(k) / gradient_bin_cnt(k) *[cos(pc2_angle) sin(pc2_angle)]];
                if locs(k) > bin_width(1)
                    m_off = locs(k) - bin_width(1) - 1;
                    cn = [cn; center + m_off*[cos(pc2_angle) sin(pc2_angle)]];
                elseif locs(k) < bin_width(1)
                    m_off = bin_width(1) - locs(k) + 1;
                    cn = [cn; center - m_off*[cos(pc2_angle) sin(pc2_angle)]];
                else
                    cn = [cn; center];
                end
            end
        end
    end
    
elseif strcmp(param.est_method, 'RANSAC') % RANSAC
    sind = sub2ind(size_im, lc(:,2), lc(:,1));    
    sind = sub2ind(size_im, lc(:,2), lc(:,1));
    im_lsc = zeros(size_im);    
    im_lsc(sind) = im_grad(sind);% .* double(im_nms(sind) > 0);
%     im_lsc_track = im_lsc;
%     im_lsc_track(sind) = im_grad(sind);
%     
%     
% %     tmp = false(size_im);
% %     tmp(sind) = true;
% %     tmp1 = im_nms .* tmp;
%     tmp = im_grad; tmp(sind) = -1; ind = tmp == -1; tmp = im_grad; tmp(~ind) = 0;
%     tmp1 = tmp;%nonmaxsup(tmp, im_dird, 1.2);
    BW = im_lsc>0;%tmp1 > 0;
    CC = bwconncomp(BW);
    for k = 1:CC.NumObjects
        if length(CC.PixelIdxList{k}) < param.thres_min_ls_length, continue; end
        %---------------------------------------------------------------------------------------
        % Get line segments via RANSAC
        ind = CC.PixelIdxList{k};
        [yy, xx] = ind2sub(size_im, ind);
        pts = [xx'; yy'];
        iterNum = size(pts,2);
        thDist = 1.2;
        thInlrRatio = .1;
        thInAngVar = 3;        
        sampleNum = 2;
        ptNum = size(pts,2);
        thInlr = round(thInlrRatio*ptNum);
        inlrNum = zeros(1,iterNum);
        inlrDirVar = zeros(1,iterNum);
        theta1 = zeros(1,iterNum);
        rho1 = zeros(1,iterNum);
        cts = zeros(iterNum, 2);
        lts = zeros(1, iterNum);
        
        for p = 1:iterNum
            % 1. fit using 2 random points
            sampleIdx = randIndex(ptNum,sampleNum);
            ptSample = pts(:,sampleIdx);
            d = ptSample(:,2)-ptSample(:,1);
            d = d/norm(d); % direction vector of the line
            
            % 2. count the inliers, if more than thInlr, refit; else iterate
            n = [-d(2),d(1)]; % unit normal vector of the line
            dist1 = n*(pts-repmat(ptSample(:,1),1,ptNum));
            inlier1 = find(abs(dist1) < thDist);
            inlrNum(p) = length(inlier1);
            inlrDirVar(p) = std(im_dird(ind(inlier1)));
            if length(inlier1) < thInlr, continue; end
            [ev, scr] = princomp(pts(:,inlier1)');
            
            d1 = ev(:,1);
            theta1(p) = -atan2(d1(2),d1(1)); % save the coefs
            rho1(p) = [-d1(2),d1(1)]*mean(pts(:,inlier1),2);
            cts(p,:) = mean(pts(:,inlier1)');
            lts(p) = sqrt(sum(d.^2));
        end
        
        % 3. choose the coef with the most inliers
        [~,idx] = max(inlrNum);
        
%         if inlrDirVar(idx) > 10, continue; end
        theta = theta1(idx);
        rho = rho1(idx);
        center = cts(idx,:);
        
        angle = -theta;
        
        [ev, scr] = princomp(pts');
        line_len = max(abs(scr(:)));
        
        line_segment(k,:) = [center line_len*2 angle];
%         line_segment(num_lsc,:) = [mean(pts), line_len, ang_val, line_score(opt_ind)];
            
%         x1 = center + lsc_len*[cos(-theta) sin(-theta)];
%         x2 = center - lsc_len*[cos(-theta) sin(-theta)];
%         plot([x1(1) x2(1)], [x1(2) x2(2)], 'r-'); % updated
    end    
elseif strcmp(param.est_method, 'RANSAC_LARGE_GRAD_MAG')
    sind = sub2ind(size_im, lc(:,2), lc(:,1));    
    tmp = im_grad; tmp(sind) = -1; ind = tmp == -1; tmp = im_grad; tmp(~ind) = 0;
    tmp1 = nonmaxsup(tmp, im_dird, 1.2);
    BW = tmp1 > 0;
    CC = bwconncomp(BW);
    
    for k = 1:CC.NumObjects
        if length(CC.PixelIdxList{k}) < param.thres_min_ls_length, continue; end
        %---------------------------------------------------------------------------------------
        % Get line segments via RANSAC
        ind = CC.PixelIdxList{k};
        grad_vec = im_grad(ind);
        [val, idx_valid_pts] = sort(grad_vec, 'descend');
        
        if length(idx_valid_pts) > 3
            idx_valid_pts = idx_valid_pts(1:ceil(length(idx_valid_pts)/3));
        end
        
        [yy, xx] = ind2sub(size_im, ind);
        pts = [xx'; yy'];
        iterNum = length(idx_valid_pts) * 2;
        thDist = 1.2;
        thInlrRatio = .1;
        thInAngVar = 3;        
        sampleNum = 2;
        ptNum = size(pts,2);
        thInlr = round(thInlrRatio*ptNum);
        inlrNum = zeros(1,iterNum);
        inlrDirVar = zeros(1,iterNum);
        theta1 = zeros(1,iterNum);
        rho1 = zeros(1,iterNum);
        cts = zeros(iterNum, 2);
        lts = zeros(1, iterNum);
        lens = zeros(1, iterNum);
        
        for p = 1:iterNum
            % 1. fit using 2 random points
            sampleIdx = randIndex(size(idx_valid_pts),sampleNum);
            ptSample = pts(:,idx_valid_pts(sampleIdx));
            d = ptSample(:,2)-ptSample(:,1);
            d = d/norm(d); % direction vector of the line
            
            % 2. count the inliers, if more than thInlr, refit; else iterate
            n = [-d(2),d(1)]; % unit normal vector of the line
            dist1 = n*(pts-repmat(ptSample(:,1),1,ptNum));
            inlier1 = find(abs(dist1) < thDist);
            inlrNum(p) = length(inlier1);
            inlrDirVar(p) = std(im_dird(ind(inlier1)));
            if length(inlier1) < thInlr, continue; end
            [ev, scr] = princomp(pts(:,inlier1)');
            
            d1 = ev(:,1);
            theta1(p) = -atan2(d1(2),d1(1)); % save the coefs
            rho1(p) = [-d1(2),d1(1)]*mean(pts(:,inlier1),2);
            cts(p,:) = mean(pts(:,inlier1)');
            lts(p) = sqrt(sum(d.^2));
            lens(p) = max(scr(:));
        end
        
        % 3. choose the coef with the most inliers
        [~,idx] = max(inlrNum);
        
%         if inlrDirVar(idx) > 10, continue; end
        theta = theta1(idx);
        rho = rho1(idx);
        center = cts(idx,:);
        
        angle = -theta;
        
        
%         lps = [lps; find(pixel_visited ~= 0)];
        line_segment = [line_segment; cts(idx,:) lens(idx) lens(idx) angle];
    end    
elseif strcmp(param.est_method, 'PROPOSED')     
    sind = sub2ind(size_im, lc(:,2), lc(:,1));
    tmp = zeros(size_im);
    tmp(sind) = im_grad(sind);
    tmp1 = tmp;
%     tmp1 = im_nms .* tmp; 
%     tmp1 = nonmaxsup(tmp, im_dird, 1.2);
        
    BW = tmp1 > 0;
    CC = bwconncomp(BW, 8);
    num_pixel_cluster = 1;
    for k = 1:CC.NumObjects
        if length(CC.PixelIdxList{k}) <= param.thres_min_ls_length, continue; end
        
        ind = CC.PixelIdxList{k};
        [yy, xx] = ind2sub(size_im, ind);
        lc = [xx yy];
        
        % Get the initial angle
        [coef, scr] = princomp(lc);
        init_angle = atan2(coef(2,1), coef(1,1));
        
        % Opposite direction
        if init_angle > 0 && init_angle < pi
            init_angle_op = init_angle - pi;
        elseif init_angle < 0 && init_angle > -pi
            init_angle_op = init_angle + pi;
        elseif init_angle == 0
            init_angle_op = pi;
        elseif init_angle == -pi || init_angle == pi
            init_angle_op = 0;
        end
        [offsets, sub_offset, bUseGX] = get_octan_rot_angle_and_offset(init_angle);%get_move_offset(init_angle);
        [offsets_op, sub_offset_op, bUseGX_op] = get_octan_rot_angle_and_offset(init_angle_op);%get_move_offset(init_angle_op);
        
        % Find the next pixel location
%         sind = sub2ind(size_im, lc(:,2), lc(:,1));
%         tmp = false(size_im);
%         tmp(sind) = true;
%         tmp1 = im_nms .* tmp;
        
        pixel_visited = false(size_im);
        
        % Start from the pixel with the maximum gradient magnitude
        if bUseGX, tmp1 = abs(im_gx .* double(tmp > 0));            
        else       tmp1 = abs(im_gy .* double(tmp > 0));
        end
        
        grad_mags = tmp1(sub2ind(size_im, yy, xx));
        [~, ind] = max(grad_mags);
        
        pixel_stack = lc(ind,:);
        
        % one direction
        while ~isempty(pixel_stack)
            cur_pixel = pixel_stack(1,:);
            pixel_stack(1,:) = [];
            pixel_visited(cur_pixel(2), cur_pixel(1)) = true;
            
            next_loc = repmat(cur_pixel, size(offsets,1), 1) + offsets;
            idx1 = find(next_loc(:,1) > 0 & next_loc(:,1) <= size_im(2) &...
                next_loc(:,2) > 0 & next_loc(:,2) <= size_im(1));
            
            
            next_loc = next_loc(idx1,:);
            
            
            next_grad = tmp1(sub2ind(size_im, next_loc(:,2), next_loc(:,1)));
            [~, idx] = max(next_grad);
            
            if ~isempty(idx) &&...
                    (next_loc(idx,1) > 0 && next_loc(idx,1) <= size_im(2)) && ...
                    (next_loc(idx,2) > 0 && next_loc(idx,2) <= size_im(1)) && ...
                    ~pixel_visited(next_loc(idx,2), next_loc(idx,1)) && tmp1(next_loc(idx,2), next_loc(idx,1))
                pixel_stack = [pixel_stack; next_loc(idx,:)];
                %             tmp1(sub2ind(size_im, next_loc(:,2), next_loc(:,1)) = 0;
            end
        end
        
        % Opposite direction
        pixel_stack = lc(ind,:);
        % one direction
        while ~isempty(pixel_stack)
            cur_pixel = pixel_stack(1,:);
            pixel_stack(1,:) = [];
            pixel_visited(cur_pixel(2), cur_pixel(1)) = true;
            
            next_loc = repmat(cur_pixel, size(offsets_op,1), 1) + offsets_op;
            idx1 = find(next_loc(:,1) > 0 & next_loc(:,1) <= size_im(2) &...
                next_loc(:,2) > 0 & next_loc(:,2) <= size_im(1));
            
            next_loc = next_loc(idx1,:);
            
            
            next_grad = tmp1(sub2ind(size_im, next_loc(:,2), next_loc(:,1)));
            [~, idx] = max(next_grad);
            
            if ~isempty(idx) &&...
                    (next_loc(idx,1) > 0 && next_loc(idx,1) <= size_im(2)) && ...
                    (next_loc(idx,2) > 0 && next_loc(idx,2) <= size_im(1)) && ...
                    ~pixel_visited(next_loc(idx,2), next_loc(idx,1)) && tmp1(next_loc(idx,2), next_loc(idx,1))
                pixel_stack = [pixel_stack; next_loc(idx,:)];
                %             tmp1(sub2ind(size_im, next_loc(:,2), next_loc(:,1)) = 0;
            end
        end
        [rr, cc] = find(pixel_visited ~= 0);
        if length(rr) <= param.thres_min_ls_length, continue; end
        
        [coef, scr] = princomp([cc rr]);
        ang_est = atan2(coef(2,1), coef(1,1));
        line_len = max(scr);
        ct = mean([cc rr]);
        
        line_segment = [line_segment; ct line_len ang_est];
        
        pixel_cluster{num_pixel_cluster,1} = find(pixel_visited ~= 0);
        num_pixel_cluster = num_pixel_cluster + 1;
    end
elseif strcmp(param.est_method, 'PROPOSED_AA_PROPOSAL')
    sind = sub2ind(size_im, lc(:,2), lc(:,1));
    tmp = false(size_im);
    tmp(sind) = true;
    im_lsc_grad = im_grad .* double(tmp);
    grad_mags = im_grad(sind);
    
    % Get the initial angle
    [coef, scr] = princomp(lc);
    init_angle = atan2(coef(2,1), coef(1,1)); 
    init_angle_op = get_opposite_angle( init_angle, 'posneg');
    line_len = abs(max(scr(:,1)) - min(scr(:,1)));
    line_segment = zeros(1, 5);
    
    % Get gradient magnitudes
    grad_vec = im_grad(sind);
    [grad_vec, srt_ind] = sort(grad_vec, 'descend');
    pt_pool = lc(srt_ind,:);
    
    % Find the first end point
    [val1, ind1] = min(scr(:,1)); pt1 = lc(ind1,:);
    p_dist = pdist2(pt1, pt_pool);
    n_ind = find(p_dist <= 2);
    pts(1,:) = pt_pool(min(n_ind),:);
    angs(1,:) = [init_angle -100];
        
    % Find the second end point
    [val2, ind2] = max(scr(:,1)); pt2 = lc(ind2,:);
    p_dist = pdist2(pt2, pt_pool);
    n_ind = find(p_dist <= 2);    
    pts(2,:) = pt_pool(min(n_ind),:);
    angs(2,:) = [init_angle_op -100]; 
    
    % Find center point
    pt3 = mean(lc);
    p_dist = pdist2(pt3, pt_pool);
    n_ind = find(p_dist <= 4);
    pts(3,:) = pt_pool(min(n_ind),:);
    angs(3,:) = [init_angle init_angle_op];
    
    line_score = zeros(3,1);
    pixel_set = cell(3,1);
    for i = 1:3
        pixel_visited = false(size_im);        
        
        % One direction
        [offsets, sub_offset, bUseGX] = get_octan_rot_angle_and_offset(angs(i,1));                
        pixel_visited = track_local_maximum(im_lsc_grad, pixel_visited, size_im, pts(i,:), offsets);        
        
        % Opposite direction
        if i == 3
            [offsets_op, sub_offset_op, bUseGX_op] = get_octan_rot_angle_and_offset(angs(i,2)); 
            pixel_visited = track_local_maximum(im_lsc_grad, pixel_visited, size_im, pts(i,:), offsets_op);            
        end        
        
        pixel_set{i} = find(pixel_visited ~= 0);
        [rr, cc] = ind2sub(size_im, pixel_set{i});
        sub_pts = [cc rr] + repmat(sub_offset, size(cc,1), 1);
        
        valid_ind = find(sub_pts(:,1) > 0 & sub_pts(:,1) <= size_im(2) & cc > 0 & cc <= size_im(2) &...
                   sub_pts(:,2) > 0 & sub_pts(:,2) <= size_im(1) & rr > 0 & rr <= size_im(1));
        sub_pts = sub_pts(valid_ind,:);
        rr = rr(valid_ind); cc = cc(valid_ind);
        
        ind_main = sub2ind(size_im, rr, cc);
        ind_sub  = sub2ind(size_im, sub_pts(:,2), sub_pts(:,1));   
        
        grad_main = im_grad(ind_main);
        grad_sub  = im_grad(ind_sub);
        line_score(i) = sum(grad_main) + sum(grad_sub) + 1 / std(grad_main + grad_sub);
        pixel_set{i} = [cc, rr];        
    end
    
    % Measure line candidate score
    [val, opt_ind] = max(line_score);        
    pts = pixel_set{opt_ind};
    line_len = pts(1,:) - pts(end,:);
    ang_val = atan2(line_len(2), line_len(1));
    line_len = sqrt(sum(line_len.^2))/2;
    
    line_segment(1:2) = mean(pts);
    line_segment(3:4) = [line_len 0];
    line_segment(5) = ang_val;
    pixel_cluster{1} = sub2ind(size_im, pts(:,2), pts(:,1));
       
    if 0
        num_pt = size(pt_pool, 1);
        num_iter = ceil(size(pt_pool, 1) * 1.5);
        line_score = zeros(num_iter, 1);
        pixel_set = cell(num_iter, 1);
        angle_set = zeros(num_iter,1);
        
        for kk = 1:num_iter
            rnd_ind = randperm(num_pt, 2);
            pt_sam1 = pt_pool(rnd_ind(1),:);
            pt_sam2 = pt_pool(rnd_ind(2),:);
            
            pt_dist = pt_sam2 - pt_sam1;
            angle_set(kk) = atan2(pt_dist(2), pt_dist(1));
            
            if pt_sam1(1) == pt_sam2(1) && pt_sam1(2) == pt_sam2(2), continue; end
            %         if abs(pt_diff(1)) <= 3 || abs(pt_diff(2)) <= 3, continue; end
            
            
            
            % Generate Antialis-Bresenham proposal
            [x, y, c] = Generate_AAB_Proposal(pt_sam1(1), pt_sam2(1), pt_sam1(2), pt_sam2(2), pt1, pt2);
            
            figure(1311); clf;
            iitmp = double(tmp);
            iitmp(sub2ind(size_im, x(:,1), y(:,1))) = 2; imagesc(iitmp);
            waitforbuttonpress;
            
            if 0
                if isempty(x), continue; end
                
                if any(isnan(x)) | any(isnan(y))
                    disp('hey')
                end
                %         figure; imagesc(tmp); axis image; hold on;
                %         plot(x(:,1), y(:,1), 'yx'); plot(x(:,2), y(:,2), 'gx');
                
                if swapped
                    [x, y] = swap_vec(x, y);
                end
                
                ind = find(x(:,1) > 0 & x(:,1) <= size_im(2) & x(:,2) > 0 & x(:,2) <= size_im(2) &...
                    y(:,1) > 0 & y(:,1) <= size_im(1) & y(:,2) > 0 & y(:,2) <= size_im(1));
                x = x(ind,:); y = y(ind, :); c = c(ind, :);
                
                ind_main = sub2ind(size_im, y(:,1), x(:,1));
                ind_sub  = sub2ind(size_im, y(:,2), x(:,2));
                
                grad_main = im_grad(ind_main);
                grad_sub  = im_grad(ind_sub);
                line_score(kk) = sum(grad_main) + sum(grad_sub) + 1 / std(grad_main + grad_sub);
                pixel_set{kk} = [x, y, c];
            end
            
        end
        [val, opt_ind] = max(line_score);
        opt_pixel_set = pixel_set{opt_ind};
        pts = [opt_pixel_set(:,1) opt_pixel_set(:,3)];
        line_len = pts(1,:) - pts(end,:);
        line_len = sqrt(sum(line_len.^2));
        
        line_segment(1:2) = mean(pts);
        line_segment(3:4) = [line_len 0];
        line_segment(5) = angle_set(opt_ind);
        pixel_cluster{1} = sub2ind(size_im, pts(:,2), pts(:,1));
    end
elseif strcmp(param.est_method, 'PROPOSED_AA_PROPOSAL_GRAD_SUM_GUIDE')
    sind = sub2ind(size_im, lc(:,2), lc(:,1));
    tmp = false(size_im);
    tmp(sind) = true;
    tmp1 = im_grad .* double(tmp);
    grad_mags = im_grad(sind);
    
    % Get the initial angle
    [coef, scr] = princomp(lc);
    init_angle = atan2(coef(2,1), coef(1,1)); 
    init_angle_op = get_opposite_angle( init_angle, 'posneg');
    line_len = abs(max(scr(:,1)) - min(scr(:,1)));
    line_segment = zeros(1, 5);
    
    % Get gradient magnitudes
    grad_vec = im_grad(sind);
    [grad_vec, srt_ind] = sort(grad_vec, 'descend');
    pt_pool = lc(srt_ind,:);
    
    % Find the first end point
    [val1, ind1] = min(scr(:,1)); pt1 = lc(ind1,:);
    p_dist = pdist2(pt1, pt_pool);
    n_ind = find(p_dist <= 2);
    pts(1,:) = pt_pool(min(n_ind),:);
    angs(1,:) = [init_angle -100];
        
    % Find the second end point
    [val2, ind2] = max(scr(:,1)); pt2 = lc(ind2,:);
    p_dist = pdist2(pt2, pt_pool);
    n_ind = find(p_dist <= 2);    
    pts(2,:) = pt_pool(min(n_ind),:);
    angs(2,:) = [init_angle_op -100]; 
    
    % Find center point
    pt3 = mean(lc);
    p_dist = pdist2(pt3, pt_pool);
    n_ind = find(p_dist <= 4);
    pts(3,:) = pt_pool(min(n_ind),:);
    angs(3,:) = [init_angle init_angle_op];
    
    line_score = zeros(3,1);
    pixel_set = cell(3,1);
    for i = 1:3
        pixel_visited = false(size_im);        
        [offsets, sub_offset, bUseGX] = get_octan_rot_angle_and_offset(angs(i,1));
        
        cal_pts = return_valid_loc(pts(i,:), size_im, sub_offset);
        if isempty(cal_pts), continue; end
        ref_grad_sum = im_grad(pts(i,2), pts(i,1)) + im_grad(pts(i,2) + sub_offset(2), pts(i,1) + sub_offset(1));
        
%         [~, ind] = max(grad_mags);        
        pixel_stack = pts(i,:);
        
        % one direction
        while ~isempty(pixel_stack)
            cur_pixel = pixel_stack(1,:);
            pixel_stack(1,:) = [];
            pixel_visited(cur_pixel(2), cur_pixel(1)) = true;
            
            % Get only relevant next position
            next_loc = repmat(cur_pixel, size(offsets,1), 1) + offsets;
            next_loc = return_valid_loc( next_loc, size_im, [] );
            
%             % Get gradient magnitude 
%             next_grad = tmp1(sub2ind(size_im, next_loc(:,2), next_loc(:,1)));
%             [vals, idx] = sort(next_grad, 'descend');
%             idx = idx(1);
            
            % Get gradient sum consistency info. motivated by antialised bresenham line
            nl = return_valid_loc( next_loc, size_im, sub_offset );
            nl_sub = nl + repmat(sub_offset, size(nl,1), 1);
            grad_sum = im_grad(sub2ind(size_im, nl(:,2), nl(:,1))) +...
                       im_grad(sub2ind(size_im, nl_sub(:,2), nl_sub(:,1)));
            
            grad_diff = abs(grad_sum - repmat(ref_grad_sum, length(grad_sum), 1));
            [~, idx] = min(grad_diff);
            
%             bCheckIntensityConsistency = false;
%             if strcmp(param.est_line_tracking_method, 'BRESENHAM_CONST_INTENSITY')
%                 bCheckIntensityConsistency = true;
%             end
            
            
            
            if ~isempty(idx) &&...
                    (next_loc(idx,1) > 0 && next_loc(idx,1) <= size_im(2)) && ...
                    (next_loc(idx,2) > 0 && next_loc(idx,2) <= size_im(1)) && ...
                    ~pixel_visited(next_loc(idx,2), next_loc(idx,1)) && tmp1(next_loc(idx,2), next_loc(idx,1))
                pixel_stack = [pixel_stack; next_loc(idx,:)];
                %             tmp1(sub2ind(size_im, next_loc(:,2), next_loc(:,1)) = 0;
            end
        end
        
        % Opposite direction
        if i == 3
            [offsets_op, sub_offset_op, bUseGX_op] = get_octan_rot_angle_and_offset(angs(i,2)); 
            pixel_stack = pts(i,:);
            
            while ~isempty(pixel_stack)
                cur_pixel = pixel_stack(1,:);
                pixel_stack(1,:) = [];
                pixel_visited(cur_pixel(2), cur_pixel(1)) = true;
                
                next_loc = repmat(cur_pixel, size(offsets_op,1), 1) + offsets_op;
                idx1 = find(next_loc(:,1) > 0 & next_loc(:,1) <= size_im(2) &...
                    next_loc(:,2) > 0 & next_loc(:,2) <= size_im(1));
                
                next_loc = next_loc(idx1,:);
                
                
                next_grad = tmp1(sub2ind(size_im, next_loc(:,2), next_loc(:,1)));
                [~, idx] = max(next_grad);
                
                if ~isempty(idx) &&...
                        (next_loc(idx,1) > 0 && next_loc(idx,1) <= size_im(2)) && ...
                        (next_loc(idx,2) > 0 && next_loc(idx,2) <= size_im(1)) && ...
                        ~pixel_visited(next_loc(idx,2), next_loc(idx,1)) && tmp1(next_loc(idx,2), next_loc(idx,1))
                    pixel_stack = [pixel_stack; next_loc(idx,:)];
                    %             tmp1(sub2ind(size_im, next_loc(:,2), next_loc(:,1)) = 0;
                end
            end
            
        end        
        pixel_set{i} = find(pixel_visited ~= 0);
        [rr, cc] = ind2sub(size_im, pixel_set{i});
        sub_pts = [cc rr] + repmat(sub_offset, size(cc,1), 1);
        
        valid_ind = find(sub_pts(:,1) > 0 & sub_pts(:,1) <= size_im(2) & cc > 0 & cc <= size_im(2) &...
                   sub_pts(:,2) > 0 & sub_pts(:,2) <= size_im(1) & rr > 0 & rr <= size_im(1));
        sub_pts = sub_pts(valid_ind,:);
        rr = rr(valid_ind); cc = cc(valid_ind);
        
        ind_main = sub2ind(size_im, rr, cc);
        ind_sub  = sub2ind(size_im, sub_pts(:,2), sub_pts(:,1));   
        
        grad_main = im_grad(ind_main);
        grad_sub  = im_grad(ind_sub);
        line_score(i) = sum(grad_main) + sum(grad_sub) + 1 / std(grad_main + grad_sub);
        pixel_set{i} = [cc, rr];
    end
    
    % Measure line candidate score
    [val, opt_ind] = max(line_score);        
    pts = pixel_set{opt_ind};
    line_len = pts(1,:) - pts(end,:);
    ang_val = atan2(line_len(2), line_len(1));
    line_len = sqrt(sum(line_len.^2))/2;
    
    line_segment(1:2) = mean(pts);
    line_segment(3:4) = [line_len 0];
    line_segment(5) = ang_val;
    pixel_cluster{1} = sub2ind(size_im, pts(:,2), pts(:,1));       
elseif strcmp(param.est_method, 'PROPOSED_AA_PROPOSAL_LSC_SEPA')
    sind = sub2ind(size_im, lc(:,2), lc(:,1));
    im_lsc = zeros(size_im);    
    im_lsc(sind) = im_grad(sind) .* double(im_nms(sind) > 0);
    im_lsc_track = im_lsc;
    im_lsc_track(sind) = im_grad(sind);
    
    % Get the initial angle and tracking offset
    [coef, scr] = princomp(lc);
    init_angle = atan2(coef(2,1), coef(1,1)); 
    init_angle_op = get_opposite_angle( init_angle, 'posneg');
    [offsets, sub_offset, bUseGX] = get_octan_rot_angle_and_offset(init_angle);
    [offsets_op, sub_offset_op, bUseGX_op] = get_octan_rot_angle_and_offset(init_angle_op);
    
    line_segment = [];
    num_lsc = 1;
    
    % Generate seed points
    BW = im_lsc > 0;
    CC = bwconncomp(BW, 8);
    im_tmp = zeros(size_im);
    for j = 1:CC.NumObjects
        ind = CC.PixelIdxList{j};
        if length(ind) < 2, continue; end
        
        cur_grad = im_lsc(ind);
        [~, max_ind] = max(cur_grad);
        [yy, xx] = ind2sub(size_im, ind);
        seed_pts = [xx(max_ind) yy(max_ind); round(median([xx yy]))];
        
        line_score = zeros(size(seed_pts,1),1);
        pixel_set = cell(size(seed_pts,1),1);
        end_points = zeros(size(seed_pts,1)*2,2);
        % Track from each seed point
        for i = 1:size(seed_pts,1)
            pixel_visited = false(size_im);
            
            % Track two directions
            [pixel_visited, end_points(i*2-1,:)] = track_local_maximum(im_lsc_track, pixel_visited, size_im, seed_pts(i,:), offsets);
            [pixel_visited, end_points(i*2,:)] = track_local_maximum(im_lsc_track, pixel_visited, size_im, seed_pts(i,:), offsets_op);
            
            vind = find(pixel_visited ~= 0);
            [rr, cc] = ind2sub(size_im, vind);
            pts = [cc rr];
            
            % Prune invalid pixel
            [pts, ~] = return_valid_loc( pts, size_im, sub_offset );
            sub_pts = pts + repmat(sub_offset, size(pts,1), 1);
            
            ind_main = sub2ind(size_im, pts(:,2), pts(:,1));
            ind_sub  = sub2ind(size_im, sub_pts(:,2), sub_pts(:,1));
            
            grad_main = im_grad(ind_main);
            grad_sub  = im_grad(ind_sub);
            
            % Calculate each lsc score (= energy)
            line_score(i) = mean(grad_main + grad_sub) + 1 / std(grad_main + grad_sub);
            pixel_set{i} = pts;
        end
        % Measure line candidate score
        [~, opt_ind] = max(line_score);
        pts = pixel_set{opt_ind};
        pts_ext = end_points(opt_ind*2-1:opt_ind*2,:);
        line_len = pts_ext(1,:) - pts_ext(end,:);
        ang_val = atan2(line_len(2), line_len(1));
        line_len = sqrt(sum(line_len.^2));
        
        if line_len < param.thres_min_ls_length, continue; end
        
        line_segment(num_lsc,:) = [mean(pts), line_len, ang_val, line_score(opt_ind)];
        pixel_cluster{num_lsc,1} = sub2ind(size_im, pts(:,2), pts(:,1));
        im_tmp(pixel_cluster{num_lsc,1}) = im_tmp(pixel_cluster{num_lsc,1}) + 1;
        num_lsc = num_lsc+1;
    end
    
    % Prune redundant lsc
    % - Detect lscs overlapping each other
    valid_lsc_idx = false(size(line_segment, 1), 1);
    lsc_visited = false(size(line_segment, 1), 1);
    
    for i = 1:length(valid_lsc_idx)
        if lsc_visited(i), continue; end
        
        overlap_idx = false(size(line_segment, 1), 1);
        
        for j = 1:length(valid_lsc_idx)
            inter_idx = intersect(pixel_cluster{i}, pixel_cluster{j});
            if isempty(inter_idx), continue; end
            overlap_idx(j) = true;
        end
        % Select the best lsc
        [~,opt_ind] = max(line_segment(overlap_idx, 5));
        tm_idx = find(overlap_idx == 1);
        valid_lsc_idx(tm_idx(opt_ind)) = true;
        
        lsc_visited(overlap_idx) = true;
    end
    
%     whos pixel_cluster 
%     if any(valid_lsc_idx == true)
%         pixel_cluster = pixel_cluster{valid_lsc_idx};
%         line_segment = line_segment(valid_lsc_idx,:);
%     end
elseif strcmp(param.est_method, 'PROPOSED_AA_SCORE')
    sind = sub2ind(size_im, lc(:,2), lc(:,1));
    im_lsc = zeros(size_im);    
    im_lsc(sind) = im_grad(sind) .* double(im_nms(sind) > 0);
    im_lsc_track = im_lsc;
    im_lsc_track(sind) = im_grad(sind);
    
    % Get the initial angle and tracking offset
    [coef, scr] = princomp(lc);
    init_angle = atan2(coef(2,1), coef(1,1)); 
    init_angle_op = get_opposite_angle( init_angle, 'posneg');
    [offsets, sub_offset, bUseGX] = get_octan_rot_angle_and_offset(init_angle);
    [offsets_op, sub_offset_op, bUseGX_op] = get_octan_rot_angle_and_offset(init_angle_op);
    
    line_segment = [];
    num_lsc = 1;    
    
    % Generate seed points
    BW = im_lsc_track > 0;
    BWL = bwlabel(im_lsc, 8);
    CC = bwconncomp(BW, 8);
    im_tmp = zeros(size_im);
    for j = 1:CC.NumObjects
        ind = CC.PixelIdxList{j};
        if length(ind) < 2, continue; end
        
        cur_grad = im_lsc(ind);
        [~, max_ind] = max(cur_grad);
        [yy, xx] = ind2sub(size_im, ind);
        seed_pts = [xx(max_ind) yy(max_ind); round(median([xx yy]))];
        
        line_score = zeros(size(seed_pts,1),1);
        pixel_set = cell(size(seed_pts,1),1);
        end_points = zeros(size(seed_pts,1)*2,2);
        % Track from each seed point
        for i = 1:size(seed_pts,1)
            pixel_visited = false(size_im);
            
            % Track two directions
            [pixel_visited, end_points(i*2-1,:)] = track_local_maximum(im_lsc_track, pixel_visited, size_im, seed_pts(i,:), offsets);
            [pixel_visited, end_points(i*2,:)] = track_local_maximum(im_lsc_track, pixel_visited, size_im, seed_pts(i,:), offsets_op);
            
            vind = find(pixel_visited ~= 0);
            [rr, cc] = ind2sub(size_im, vind);
            pts = [cc rr];
            num_tracked = size(pts, 1);
            
            % Prune invalid pixel
            [pts, ~] = return_valid_loc( pts, size_im, sub_offset );
            sub_pts = pts + repmat(sub_offset, size(pts,1), 1);
            
            ind_main = sub2ind(size_im, pts(:,2), pts(:,1));
            ind_sub  = sub2ind(size_im, sub_pts(:,2), sub_pts(:,1));
            
            grad_main = im_grad(ind_main);
            grad_sub  = im_grad(ind_sub);
            
            % Calculate each lsc score (= energy)
            [pro_x, pro_y] = bresenham(end_points(i*2-1, 1), end_points(i*2-1, 2), end_points(i*2, 1), end_points(i*2, 2));
            num_model = size(pro_x,1);
            
            intersec_idx = intersect([pro_x, pro_y], pts, 'rows');
            num_inter = size(intersec_idx, 1);
            
            if bUseGX
                hist_pt = [pts(:,2) pts(:,1)];                
            else
                hist_pt = pts;
            end
            
            [unq_val, num_unq] = count_unique(pts(:,2));
            lsc_hist = zeros(1, length(unq_val));
            for k = 1:length(unq_val)
                valid_idx = find(pts(:,2) == unq_val(k));
                valid_idx = sub2ind(size_im, pts(valid_idx,2), pts(valid_idx,1));
                lsc_hist(k) = sum(im_lsc_track(valid_idx)) / length(valid_idx);
            end
            
%         num_inter / num_model
%         figure(555);
%         tmp = zeros(size_im);
% %         tmp(sind) = 1;
%         tmp(ind_main) = 2;
%         tmp(sub2ind(size_im, pro_y, pro_x)) = 3;
%         imagesc(tmp);%imshowpair(im_gray, pixel_visited, 'blend');
%         waitforbuttonpress;
            
%             if (num_inter / num_model) < .5
%                 line_score(i) = 0;
% %                 0
%             else
                line_score(i) = num_inter^2 / (num_tracked * num_model) + 1 / std(grad_main + grad_sub);
%             end
            pixel_set{i} = pts;
        end
        % Measure line candidate score
        [opt_val, opt_ind] = max(line_score);
        pts = pixel_set{opt_ind};
        pts_ext = end_points(opt_ind*2-1:opt_ind*2,:);
        line_len = pts_ext(1,:) - pts_ext(end,:);
        ang_val = atan2(line_len(2), line_len(1));
        line_len = sqrt(sum(line_len.^2));
        
        if line_len < param.thres_min_ls_length || opt_val == 0, continue; end
        
        line_segment(num_lsc,:) = [mean(pts), line_len, ang_val, line_score(opt_ind)];
        pixel_cluster{num_lsc,1} = sub2ind(size_im, pts(:,2), pts(:,1));
        im_tmp(pixel_cluster{num_lsc,1}) = im_tmp(pixel_cluster{num_lsc,1}) + 1;
        num_lsc = num_lsc+1;
    end
    
    % Prune redundant lsc
    % - Detect lscs overlapping each other
    valid_lsc_idx = false(size(line_segment, 1), 1);
    lsc_visited = false(size(line_segment, 1), 1);
elseif strcmp(param.est_method, 'PROPOSED_NMS_LSC_SEPA')
    sind = sub2ind(size_im, lc(:,2), lc(:,1));
    im_lsc = zeros(size_im);    
    im_lsc(sind) = im_grad(sind) .* double(im_nms(sind) > 0);
    im_lsc_track = im_lsc;
    im_lsc_track(sind) = im_grad(sind);
    
    % Get the initial angle and tracking offset
    line_segment = [];
    num_lsc = 1;    
    
    % Generate seed points
    BW = im_lsc_track > 0; % im_lsc or im_lsc_track
    BWL = bwlabel(im_lsc, 8);
    CC = bwconncomp(BW, 8);
    im_tmp = zeros(size_im);
    for j = 1:CC.NumObjects
        ind = CC.PixelIdxList{j};
        if length(ind) < 2, continue; end
        
        cur_grad = im_lsc(ind);
        [~, max_ind] = max(cur_grad);
        [yy, xx] = ind2sub(size_im, ind);
        seed_pts = [xx(max_ind) yy(max_ind); round(median([xx yy]))];
        
        [coef, scr] = princomp([xx yy]);
        init_angle = atan2(coef(2,1), coef(1,1));
        init_angle_op = get_opposite_angle( init_angle, 'posneg');
        [offsets, sub_offset, bUseGX] = get_octan_rot_angle_and_offset(init_angle);
        [offsets_op, sub_offset_op, bUseGX_op] = get_octan_rot_angle_and_offset(init_angle_op);
    
        
        line_score = zeros(size(seed_pts,1),1);
        pixel_set = cell(size(seed_pts,1),1);
        end_points = zeros(size(seed_pts,1)*2,2);
        % Track from each seed point
        for i = 1:size(seed_pts,1)
            pixel_visited = false(size_im);
            
            % Track two directions
            [pixel_visited, end_points(i*2-1,:)] = track_local_maximum(im_lsc_track, pixel_visited, size_im, seed_pts(i,:), offsets);
            [pixel_visited, end_points(i*2,:)] = track_local_maximum(im_lsc_track, pixel_visited, size_im, seed_pts(i,:), offsets_op);
            
            vind = find(pixel_visited ~= 0);
            [rr, cc] = ind2sub(size_im, vind);
            pts = [cc rr];
            num_tracked = size(pts, 1);
            
            % Prune invalid pixel
            [pts, ~] = return_valid_loc( pts, size_im, sub_offset );
            sub_pts = pts + repmat(sub_offset, size(pts,1), 1);
            
            if size(pts,1) == 1, continue; end
            
            ind_main = sub2ind(size_im, pts(:,2), pts(:,1));
            ind_sub  = sub2ind(size_im, sub_pts(:,2), sub_pts(:,1));
            
            grad_main = im_grad(ind_main);
            grad_sub  = im_grad(ind_sub);
            
            % Calculate each lsc score (= energy)
            [pro_x, pro_y] = bresenham(end_points(i*2-1, 1), end_points(i*2-1, 2), end_points(i*2, 1), end_points(i*2, 2));
            num_model = size(pro_x,1);
            
            intersec_idx = intersect([pro_x, pro_y], pts, 'rows');
            num_inter = size(intersec_idx, 1);
            
            if bUseGX
                hist_pt = [pts(:,2) pts(:,1)];                
            else
                hist_pt = pts;
            end
            
            [unq_val, num_unq] = count_unique(pts(:,2));
            lsc_hist = zeros(1, length(unq_val));
            for k = 1:length(unq_val)
                valid_idx = find(pts(:,2) == unq_val(k));
                valid_idx = sub2ind(size_im, pts(valid_idx,2), pts(valid_idx,1));
                lsc_hist(k) = sum(im_lsc_track(valid_idx)) / length(valid_idx);
            end
            
%         num_inter / num_model
%         figure(555);
%         tmp = zeros(size_im);
% %         tmp(sind) = 1;
%         tmp(ind_main) = 2;
%         tmp(sub2ind(size_im, pro_y, pro_x)) = 3;
%         imagesc(tmp);%imshowpair(im_gray, pixel_visited, 'blend');
%         waitforbuttonpress;
            
%             if (num_inter / num_model) < .5
%                 line_score(i) = 0;
% %                 0
%             else
                line_score(i) = num_inter^2 / (num_tracked * num_model) + 1 / std(grad_main + grad_sub);
%             end
            pixel_set{i} = pts;
        end
        % Measure line candidate score
        [opt_val, opt_ind] = max(line_score);
        pts = pixel_set{opt_ind};
        pts_ext = end_points(opt_ind*2-1:opt_ind*2,:);
        line_len = pts_ext(1,:) - pts_ext(end,:);
        ang_val = atan2(line_len(2), line_len(1));
        line_len = sqrt(sum(line_len.^2));
        
        if line_len < param.thres_min_ls_length || opt_val == 0 || opt_val == Inf || size(pts,1) == 1 || isempty(pts), continue; end
        
%         [yy, xx] = ind2sub(
%         pts = [xx yy];
        a = tan(ang_val);
        b = -1;
        c = -a*pts(1,1) + pts(1,2);
        dist1 = abs(a * pts(:,1) + b * pts(:,2) + c) / sqrt(a^2 + b^2);
        
        [val, ind] = max(dist1);
        
        
        if val >= 1.5
            % split if the maximum distance exceeds the threshold
            pts1 = [pts(1,:); pts(ind,:)];
            line_len1 = pts1(1,:) - pts1(end,:);
            ang_val1 = atan2(line_len1(2), line_len1(1));
            line_len1 = sqrt(sum(line_len1.^2));
                        
            line_segment(num_lsc,:) = [mean(pts1), line_len1, ang_val1, line_score(opt_ind)];
            pixel_cluster{num_lsc,1} = sub2ind(size_im, pts(1:ind,2), pts(1:ind,1));
            im_tmp(pixel_cluster{num_lsc,1}) = im_tmp(pixel_cluster{num_lsc,1}) + 1;            
            num_lsc = num_lsc+1;
            
            pts1 = [pts(ind,:); pts(end,:)];
            line_len1 = pts1(1,:) - pts1(end,:);
            ang_val1 = atan2(line_len1(2), line_len1(1));
            line_len1 = sqrt(sum(line_len1.^2));
            
            line_segment(num_lsc,:) = [mean(pts1), line_len1, ang_val1, line_score(opt_ind)];   
            pixel_cluster{num_lsc,1} = sub2ind(size_im, pts(ind:end,2), pts(ind:end,1));
            im_tmp(pixel_cluster{num_lsc,1}) = im_tmp(pixel_cluster{num_lsc,1}) + 1;            
            num_lsc = num_lsc+1;
        else
        
            line_segment(num_lsc,:) = [mean(pts), line_len, ang_val, line_score(opt_ind)];
            pixel_cluster{num_lsc,1} = sub2ind(size_im, pts(:,2), pts(:,1));
            im_tmp(pixel_cluster{num_lsc,1}) = im_tmp(pixel_cluster{num_lsc,1}) + 1;            
            num_lsc = num_lsc+1;
        end
    end
elseif strcmp(param.est_method, 'PROPOSED_NMS_SEEDING')    
    sind = sub2ind(size_im, lc(:,2), lc(:,1));
    im_lsc = zeros(size_im);    
    im_lsc(sind) = im_grad(sind) .* double(im_nms(sind) > 0);
    im_lsc_track = im_lsc;
    im_lsc_track(sind) = im_grad(sind);
    
    % Get the initial angle and tracking offset
    line_segment = [];
    num_lsc = 1;    
    
    % Generate seed points
    BW = im_lsc_track > 0; % im_lsc or im_lsc_track
    BWL = bwlabel(im_lsc, 8);
    CC = bwconncomp(BW, 8);
    im_tmp = zeros(size_im);
    for j = 1:CC.NumObjects
        ind = CC.PixelIdxList{j};
        if length(ind) < 2, continue; end
        
        cur_grad = im_lsc(ind);
        [~, max_ind] = max(cur_grad);
        [yy, xx] = ind2sub(size_im, ind);
        seed_pts = [xx(max_ind) yy(max_ind); round(median([xx yy]))];
        
        [coef, scr] = princomp([xx yy]);
        init_angle = atan2(coef(2,1), coef(1,1));
        init_angle_op = get_opposite_angle( init_angle, 'posneg');
        [offsets, sub_offset, bUseGX] = get_octan_rot_angle_and_offset(init_angle);
        [offsets_op, sub_offset_op, bUseGX_op] = get_octan_rot_angle_and_offset(init_angle_op);
    
        
        line_score = zeros(size(seed_pts,1),1);
        pixel_set = cell(size(seed_pts,1),1);
        end_points = zeros(size(seed_pts,1)*2,2);
        % Track from each seed point
        for i = 1:size(seed_pts,1)
            pixel_visited = false(size_im);
            
            % Track two directions
            [pixel_visited, end_points(i*2-1,:)] = track_local_maximum(im_lsc_track, pixel_visited, size_im, seed_pts(i,:), offsets);
            [pixel_visited, end_points(i*2,:)] = track_local_maximum(im_lsc_track, pixel_visited, size_im, seed_pts(i,:), offsets_op);
            
            vind = find(pixel_visited ~= 0);
            [rr, cc] = ind2sub(size_im, vind);
            pts = [cc rr];
            num_tracked = size(pts, 1);
            
            % Prune invalid pixel
            [pts, ~] = return_valid_loc( pts, size_im, sub_offset );
            sub_pts = pts + repmat(sub_offset, size(pts,1), 1);
            
            if size(pts,1) == 1, continue; end
            
            ind_main = sub2ind(size_im, pts(:,2), pts(:,1));
            ind_sub  = sub2ind(size_im, sub_pts(:,2), sub_pts(:,1));
            
            grad_main = im_grad(ind_main);
            grad_sub  = im_grad(ind_sub);
            
            % Calculate each lsc score (= energy)
            [pro_x, pro_y] = bresenham(end_points(i*2-1, 1), end_points(i*2-1, 2), end_points(i*2, 1), end_points(i*2, 2));
            num_model = size(pro_x,1);
            
            intersec_idx = intersect([pro_x, pro_y], pts, 'rows');
            num_inter = size(intersec_idx, 1);
            
            if bUseGX
                hist_pt = [pts(:,2) pts(:,1)];                
            else
                hist_pt = pts;
            end
            
            [unq_val, num_unq] = count_unique(pts(:,2));
            lsc_hist = zeros(1, length(unq_val));
            for k = 1:length(unq_val)
                valid_idx = find(pts(:,2) == unq_val(k));
                valid_idx = sub2ind(size_im, pts(valid_idx,2), pts(valid_idx,1));
                lsc_hist(k) = sum(im_lsc_track(valid_idx)) / length(valid_idx);
            end
            
%         num_inter / num_model
%         figure(555);
%         tmp = zeros(size_im);
% %         tmp(sind) = 1;
%         tmp(ind_main) = 2;
%         tmp(sub2ind(size_im, pro_y, pro_x)) = 3;
%         imagesc(tmp);%imshowpair(im_gray, pixel_visited, 'blend');
%         waitforbuttonpress;
            
%             if (num_inter / num_model) < .5
%                 line_score(i) = 0;
% %                 0
%             else
                line_score(i) = num_inter^2 / (num_tracked * num_model) + 1 / std(grad_main + grad_sub);
%             end
            pixel_set{i} = pts;
        end
        % Measure line candidate score
        [opt_val, opt_ind] = max(line_score);
        pts = pixel_set{opt_ind};
        pts_ext = end_points(opt_ind*2-1:opt_ind*2,:);
        line_len = pts_ext(1,:) - pts_ext(end,:);
        ang_val = atan2(line_len(2), line_len(1));
        line_len = sqrt(sum(line_len.^2));
        
        if line_len < param.thres_min_ls_length || opt_val == 0 || opt_val == Inf || size(pts,1) == 1 || isempty(pts), continue; end
        
%         [yy, xx] = ind2sub(
%         pts = [xx yy];
        a = tan(ang_val);
        b = -1;
        c = -a*pts(1,1) + pts(1,2);
        dist1 = abs(a * pts(:,1) + b * pts(:,2) + c) / sqrt(a^2 + b^2);
        
        [val, ind] = max(dist1);
        
        
        if val >= 1.5
            % split if the maximum distance exceeds the threshold
            pts1 = [pts(1,:); pts(ind,:)];
            line_len1 = pts1(1,:) - pts1(end,:);
            ang_val1 = atan2(line_len1(2), line_len1(1));
            line_len1 = sqrt(sum(line_len1.^2));
                        
            line_segment(num_lsc,:) = [mean(pts1), line_len1, ang_val1, line_score(opt_ind)];
            pixel_cluster{num_lsc,1} = sub2ind(size_im, pts(1:ind,2), pts(1:ind,1));
            im_tmp(pixel_cluster{num_lsc,1}) = im_tmp(pixel_cluster{num_lsc,1}) + 1;            
            num_lsc = num_lsc+1;
            
            pts1 = [pts(ind,:); pts(end,:)];
            line_len1 = pts1(1,:) - pts1(end,:);
            ang_val1 = atan2(line_len1(2), line_len1(1));
            line_len1 = sqrt(sum(line_len1.^2));
            
            line_segment(num_lsc,:) = [mean(pts1), line_len1, ang_val1, line_score(opt_ind)];   
            pixel_cluster{num_lsc,1} = sub2ind(size_im, pts(ind:end,2), pts(ind:end,1));
            im_tmp(pixel_cluster{num_lsc,1}) = im_tmp(pixel_cluster{num_lsc,1}) + 1;            
            num_lsc = num_lsc+1;
        else
        
            line_segment(num_lsc,:) = [mean(pts), line_len, ang_val, line_score(opt_ind)];
            pixel_cluster{num_lsc,1} = sub2ind(size_im, pts(:,2), pts(:,1));
            im_tmp(pixel_cluster{num_lsc,1}) = im_tmp(pixel_cluster{num_lsc,1}) + 1;            
            num_lsc = num_lsc+1;
        end
    end
elseif strcmp(param.est_method, 'PROPOSED_LSC_SEPA')
    sind = sub2ind(size_im, lc(:,2), lc(:,1));
    im_lsc = zeros(size_im);    
    im_lsc(sind) = im_grad(sind);
    im_lsc_track = im_lsc;
    
    % Get the initial angle and tracking offset
    line_segment = [];
    num_lsc = 1;    
    
    lc_grad = im_grad(sind);
    
    
    % Generate seed points
    BW = im_lsc_track > 0; % im_lsc or im_lsc_track
    BWL = bwlabel(im_lsc, 8);
    CC = bwconncomp(BW, 8);
    im_tmp = zeros(size_im);
    for j = 1:CC.NumObjects
        ind = CC.PixelIdxList{j};
        if length(ind) < 2, continue; end
        
        cur_grad = im_lsc(ind);
        [~, max_ind] = max(cur_grad);
        [yy, xx] = ind2sub(size_im, ind);
        seed_pts = [xx(max_ind) yy(max_ind); round(median([xx yy]))];
        
        [coef, scr] = princomp([xx yy]);
        init_angle = atan2(coef(2,1), coef(1,1));
        init_angle_op = get_opposite_angle( init_angle, 'posneg');
        [offsets, sub_offset, bUseGX] = get_octan_rot_angle_and_offset(init_angle);
        [offsets_op, sub_offset_op, bUseGX_op] = get_octan_rot_angle_and_offset(init_angle_op);
    
        
        line_score = zeros(size(seed_pts,1),1);
        pixel_set = cell(size(seed_pts,1),1);
        end_points = zeros(size(seed_pts,1)*2,2);
        % Track from each seed point
        for i = 1:size(seed_pts,1)
            pixel_visited = false(size_im);
            
            % Track two directions
            [pixel_visited, end_points(i*2-1,:)] = track_local_maximum(im_lsc_track, pixel_visited, size_im, seed_pts(i,:), offsets);
            [pixel_visited, end_points(i*2,:)] = track_local_maximum(im_lsc_track, pixel_visited, size_im, seed_pts(i,:), offsets_op);
            
            vind = find(pixel_visited ~= 0);
            [rr, cc] = ind2sub(size_im, vind);
            pts = [cc rr];
            num_tracked = size(pts, 1);
            
            % Prune invalid pixel
            [pts, ~] = return_valid_loc( pts, size_im, sub_offset );
            sub_pts = pts + repmat(sub_offset, size(pts,1), 1);
            
            if size(pts,1) == 1, continue; end
            
            ind_main = sub2ind(size_im, pts(:,2), pts(:,1));
            ind_sub  = sub2ind(size_im, sub_pts(:,2), sub_pts(:,1));
            
            grad_main = im_grad(ind_main);
            grad_sub  = im_grad(ind_sub);
            
            % Calculate each lsc score (= energy)
            [pro_x, pro_y] = bresenham(end_points(i*2-1, 1), end_points(i*2-1, 2), end_points(i*2, 1), end_points(i*2, 2));
            num_model = size(pro_x,1);
            
            intersec_idx = intersect([pro_x, pro_y], pts, 'rows');
            num_inter = size(intersec_idx, 1);
            
            if bUseGX
                hist_pt = [pts(:,2) pts(:,1)];                
            else
                hist_pt = pts;
            end
            
            [unq_val, num_unq] = count_unique(pts(:,2));
            lsc_hist = zeros(1, length(unq_val));
            for k = 1:length(unq_val)
                valid_idx = find(pts(:,2) == unq_val(k));
                valid_idx = sub2ind(size_im, pts(valid_idx,2), pts(valid_idx,1));
                lsc_hist(k) = sum(im_lsc_track(valid_idx)) / length(valid_idx);
            end
            
%         num_inter / num_model
%         figure(555);
%         tmp = zeros(size_im);
% %         tmp(sind) = 1;
%         tmp(ind_main) = 2;
%         tmp(sub2ind(size_im, pro_y, pro_x)) = 3;
%         imagesc(tmp);%imshowpair(im_gray, pixel_visited, 'blend');
%         waitforbuttonpress;
            
%             if (num_inter / num_model) < .5
%                 line_score(i) = 0;
% %                 0
%             else
                line_score(i) = num_inter^2 / (num_tracked * num_model) + 1 / std(grad_main + grad_sub);
%             end
            pixel_set{i} = pts;
        end
        % Measure line candidate score
        [opt_val, opt_ind] = max(line_score);
        pts = pixel_set{opt_ind};
        pts_ext = end_points(opt_ind*2-1:opt_ind*2,:);
        line_len = pts_ext(1,:) - pts_ext(end,:);
        ang_val = atan2(line_len(2), line_len(1));
        line_len = sqrt(sum(line_len.^2));
        
        if line_len < param.thres_min_ls_length || opt_val == 0 || opt_val == Inf || size(pts,1) == 1 || isempty(pts), continue; end
        
%         [yy, xx] = ind2sub(
%         pts = [xx yy];
        a = tan(ang_val);
        b = -1;
        c = -a*pts(1,1) + pts(1,2);
        dist1 = abs(a * pts(:,1) + b * pts(:,2) + c) / sqrt(a^2 + b^2);
        
        [val, ind] = max(dist1);
        
        
        if val >= 1.5
            % split if the maximum distance exceeds the threshold
            pts1 = [pts(1,:); pts(ind,:)];
            line_len1 = pts1(1,:) - pts1(end,:);
            ang_val1 = atan2(line_len1(2), line_len1(1));
            line_len1 = sqrt(sum(line_len1.^2));
                        
            line_segment(num_lsc,:) = [mean(pts1), line_len1, ang_val1, line_score(opt_ind)];
            pixel_cluster{num_lsc,1} = sub2ind(size_im, pts(1:ind,2), pts(1:ind,1));
            im_tmp(pixel_cluster{num_lsc,1}) = im_tmp(pixel_cluster{num_lsc,1}) + 1;            
            num_lsc = num_lsc+1;
            
            pts1 = [pts(ind,:); pts(end,:)];
            line_len1 = pts1(1,:) - pts1(end,:);
            ang_val1 = atan2(line_len1(2), line_len1(1));
            line_len1 = sqrt(sum(line_len1.^2));
            
            line_segment(num_lsc,:) = [mean(pts1), line_len1, ang_val1, line_score(opt_ind)];   
            pixel_cluster{num_lsc,1} = sub2ind(size_im, pts(ind:end,2), pts(ind:end,1));
            im_tmp(pixel_cluster{num_lsc,1}) = im_tmp(pixel_cluster{num_lsc,1}) + 1;            
            num_lsc = num_lsc+1;
        else
        
            line_segment(num_lsc,:) = [mean(pts), line_len, ang_val, line_score(opt_ind)];
            pixel_cluster{num_lsc,1} = sub2ind(size_im, pts(:,2), pts(:,1));
            im_tmp(pixel_cluster{num_lsc,1}) = im_tmp(pixel_cluster{num_lsc,1}) + 1;            
            num_lsc = num_lsc+1;
        end
    end
elseif strcmp(param.est_method, 'PROPOSED_TRACK_LM')
    sind = sub2ind(size_im, lc(:,2), lc(:,1));
    im_lsc = zeros(size_im);    
    im_lsc(sind) = im_grad(sind);
    im_lsc_track = im_lsc;
    
    % Get the initial angle and tracking offset
    line_segment = [];
    num_lsc = 1;    
    
    lc_grad = im_grad(sind);
    
    
    % Generate seed points
    BW = im_lsc_track > 0; % im_lsc or im_lsc_track
    BWL = bwlabel(im_lsc, 8);
    CC = bwconncomp(BW, 8);
    im_tmp = zeros(size_im);
    for j = 1:CC.NumObjects
        ind = CC.PixelIdxList{j};
        if length(ind) < 2, continue; end
        
        cur_grad = im_lsc(ind);
        [~, max_ind] = max(cur_grad);
        [yy, xx] = ind2sub(size_im, ind);
        seed_pts = [xx(max_ind) yy(max_ind); round(median([xx yy]))];
        
        [coef, scr] = princomp([xx yy]);
        init_angle = atan2(coef(2,1), coef(1,1));
        init_angle_op = get_opposite_angle( init_angle, 'posneg');
        [offsets, sub_offset, bUseGX] = get_octan_rot_angle_and_offset(init_angle);
        [offsets_op, sub_offset_op, bUseGX_op] = get_octan_rot_angle_and_offset(init_angle_op);
    
        
        line_score = zeros(size(seed_pts,1),1);
        pixel_set = cell(size(seed_pts,1),1);
        end_points = zeros(size(seed_pts,1)*2,2);
        % Track from each seed point
        for i = 1:size(seed_pts,1)
            pixel_visited = false(size_im);
            
            % Track two directions
            [pixel_visited, end_points(i*2-1,:)] = track_local_maximum(im_lsc_track, pixel_visited, size_im, seed_pts(i,:), offsets);
            [pixel_visited, end_points(i*2,:)] = track_local_maximum(im_lsc_track, pixel_visited, size_im, seed_pts(i,:), offsets_op);
            
            vind = find(pixel_visited ~= 0);
            [rr, cc] = ind2sub(size_im, vind);
            pts = [cc rr];
            num_tracked = size(pts, 1);
            
            % Prune invalid pixel
            [pts, ~] = return_valid_loc( pts, size_im, sub_offset );
            sub_pts = pts + repmat(sub_offset, size(pts,1), 1);
            
            if size(pts,1) == 1, continue; end
            
            ind_main = sub2ind(size_im, pts(:,2), pts(:,1));
            ind_sub  = sub2ind(size_im, sub_pts(:,2), sub_pts(:,1));
            
            grad_main = im_grad(ind_main);
            grad_sub  = im_grad(ind_sub);
            
            % Calculate each lsc score (= energy)
            [pro_x, pro_y] = bresenham(end_points(i*2-1, 1), end_points(i*2-1, 2), end_points(i*2, 1), end_points(i*2, 2));
            num_model = size(pro_x,1);
            
            intersec_idx = intersect([pro_x, pro_y], pts, 'rows');
            num_inter = size(intersec_idx, 1);
            
            if bUseGX
                hist_pt = [pts(:,2) pts(:,1)];                
            else
                hist_pt = pts;
            end
            
            [unq_val, num_unq] = count_unique(pts(:,2));
            lsc_hist = zeros(1, length(unq_val));
            for k = 1:length(unq_val)
                valid_idx = find(pts(:,2) == unq_val(k));
                valid_idx = sub2ind(size_im, pts(valid_idx,2), pts(valid_idx,1));
                lsc_hist(k) = sum(im_lsc_track(valid_idx)) / length(valid_idx);
            end
            
%         num_inter / num_model
%         figure(555);
%         tmp = zeros(size_im);
% %         tmp(sind) = 1;
%         tmp(ind_main) = 2;
%         tmp(sub2ind(size_im, pro_y, pro_x)) = 3;
%         imagesc(tmp);%imshowpair(im_gray, pixel_visited, 'blend');
%         waitforbuttonpress;
            
%             if (num_inter / num_model) < .5
%                 line_score(i) = 0;
% %                 0
%             else
                line_score(i) = num_inter^2 / (num_tracked * num_model) + 1 / std(grad_main + grad_sub);
%             end
            pixel_set{i} = pts;
        end
        % Measure line candidate score
        [opt_val, opt_ind] = max(line_score);
        pts = pixel_set{opt_ind};
        pts_ext = end_points(opt_ind*2-1:opt_ind*2,:);
        line_len = pts_ext(1,:) - pts_ext(end,:);
        ang_val = atan2(line_len(2), line_len(1));
        line_len = sqrt(sum(line_len.^2));
        
        if line_len < param.thres_min_ls_length || opt_val == 0 || opt_val == Inf || size(pts,1) == 1 || isempty(pts), continue; end
        
%         [yy, xx] = ind2sub(
%         pts = [xx yy];
        a = tan(ang_val);
        b = -1;
        c = -a*pts(1,1) + pts(1,2);
        dist1 = abs(a * pts(:,1) + b * pts(:,2) + c) / sqrt(a^2 + b^2);
        
        [val, ind] = max(dist1);
        
        
        if val >= 1.5
            % split if the maximum distance exceeds the threshold
            pts1 = [pts(1,:); pts(ind,:)];
            line_len1 = pts1(1,:) - pts1(end,:);
            ang_val1 = atan2(line_len1(2), line_len1(1));
            line_len1 = sqrt(sum(line_len1.^2));
                        
            line_segment(num_lsc,:) = [mean(pts1), line_len1, ang_val1, line_score(opt_ind)];
            pixel_cluster{num_lsc,1} = sub2ind(size_im, pts(1:ind,2), pts(1:ind,1));
            im_tmp(pixel_cluster{num_lsc,1}) = im_tmp(pixel_cluster{num_lsc,1}) + 1;            
            num_lsc = num_lsc+1;
            
            pts1 = [pts(ind,:); pts(end,:)];
            line_len1 = pts1(1,:) - pts1(end,:);
            ang_val1 = atan2(line_len1(2), line_len1(1));
            line_len1 = sqrt(sum(line_len1.^2));
            
            line_segment(num_lsc,:) = [mean(pts1), line_len1, ang_val1, line_score(opt_ind)];   
            pixel_cluster{num_lsc,1} = sub2ind(size_im, pts(ind:end,2), pts(ind:end,1));
            im_tmp(pixel_cluster{num_lsc,1}) = im_tmp(pixel_cluster{num_lsc,1}) + 1;            
            num_lsc = num_lsc+1;
        else
        
            line_segment(num_lsc,:) = [mean(pts), line_len, ang_val, line_score(opt_ind)];
            pixel_cluster{num_lsc,1} = sub2ind(size_im, pts(:,2), pts(:,1));
            im_tmp(pixel_cluster{num_lsc,1}) = im_tmp(pixel_cluster{num_lsc,1}) + 1;            
            num_lsc = num_lsc+1;
        end
    end    
end