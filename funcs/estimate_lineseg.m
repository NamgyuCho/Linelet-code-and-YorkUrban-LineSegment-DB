function [ line_seg ] = estimate_lineseg( ll_cand, str_opt )
    %ESTIMATE_LINESEG Summary of this function goes here
    %   Detailed explanation goes here

    % str_opt \in {'major_linelet', 'sub_major_linelet', 'rect'}
    if strcmp(str_opt, 'major_signle_linelet')
        % Pick the longest linelet
        [val, idx] = max(ll_cand(:,4));

        % -------------- Get the line equation centered and angled wrt the longest one ----------------------- %    
        % Angle
        if      ll_cand(end,3) - ll_cand(1,3) > 0,   ang = atan2(1, val);
        elseif  ll_cand(end,3) - ll_cand(1,3) < 0,   ang = atan2(-1, val);        
        else                                         ang = 0;
        end

        % Center
        ct = [mean(ll_cand(idx,1:2)), ll_cand(idx,3)];

        % Length
        pair_dist = pdist2(ct, [ll_cand(1,1) ll_cand(1,3); ll_cand(end,2) ll_cand(end,3)]);

        % Move the center to the real one
        ct = [ll_cand(1,1) ll_cand(1,3)] + diff(pair_dist)/2*[cos(ang) sin(ang)];

        line_seg = [ct, ang, pair_dist(1) + pair_dist(end)];
    end
    
    % When try to pick the most influential (more than two) linelets
    if strcmp(str_opt, 'major_group_linelet')
        
        if size(ll_cand, 1) > 3
            if ll_cand(1,4) >= ll_cand(2,4),    idx_1 = 1;
            else                                idx_1 = 2;
            end

            if ll_cand(end,4) >= ll_cand(end-1,4),  idx_2 = size(ll_cand,1);
            else                                    idx_2 = size(ll_cand,1) - 1;
            end
        else
            idx_1 = 1;
            idx_2 = size(ll_cand, 1);
        end
        
        % Angle
        ang = atan2(ll_cand(idx_2,3) - ll_cand(idx_1,3), ll_cand(idx_2,2) - ll_cand(idx_1,1));
        
        % Center
        ct = [ll_cand(idx_1,1) + ll_cand(idx_2,2) ll_cand(idx_1,3) + ll_cand(idx_2,3)] / 2;
        
        % Length
        pair_dist = pdist2(ct, [ll_cand(1,1) ll_cand(1,3); ll_cand(end,2) ll_cand(end,3)]);
        
        % Move the center to the real one
        ct = ct + diff(pair_dist)/2*[cos(ang) sin(ang)];
        
        line_seg = [ct, ang, pair_dist(1) + pair_dist(end)];   
    end
    
    if strcmp(str_opt, 'sub_major_linelet')
        if size(ll_cand, 1) > 3
            if ll_cand(1,4) >= ll_cand(2,4),    idx_1 = 1;
            else                                idx_1 = 2;
            end

            if ll_cand(end,4) >= ll_cand(end-1,4),  idx_2 = size(ll_cand,1);
            else                                    idx_2 = size(ll_cand,1) - 1;
            end
        else
            idx_1 = 1;
            idx_2 = size(ll_cand, 1);
        end
        
        % Angle
        ang = atan2(ll_cand(idx_2,3) - ll_cand(idx_1,3), ll_cand(idx_2,2) - ll_cand(idx_1,1));
        
        % Center
        ct = [ll_cand(idx_1,1) + ll_cand(idx_2,2) ll_cand(idx_1,3) + ll_cand(idx_2,3)] / 2;
        
        % Length
        pair_dist = pdist2(ct, [ll_cand(1,1) ll_cand(1,3); ll_cand(end,2) ll_cand(end,3)]);
        
        % Move the center to the real one
        ct = ct + diff(pair_dist)/2*[cos(ang) sin(ang)];
        
        line_seg = [ct, ang, pair_dist(1) + pair_dist(end)];     
    end
          
    if strcmp(str_opt, 'rect')
        idx_1 = 1;
        idx_2 = size(ll_cand, 1);
        
        % Angle
        ang = atan2(ll_cand(idx_2,3) - ll_cand(idx_1,3), ll_cand(idx_2,2) - ll_cand(idx_1,1));
        
        % Center
        ct = [ll_cand(idx_1,1) + ll_cand(idx_2,2) ll_cand(idx_1,3) + ll_cand(idx_2,3)] / 2;
        
        % Length
        pair_dist = pdist2(ct, [ll_cand(1,1) ll_cand(1,3); ll_cand(end,2) ll_cand(end,3)]);
        
        % Move the center to the real one
        ct = ct + diff(pair_dist)/2*[cos(ang) sin(ang)];
        
        line_seg = [ct, ang, pair_dist(1) + pair_dist(end)];      
    end
    
%     figure; im_tmp = DrawLL(ll_cand, size_im); imagesc(im_tmp); hold on;
%     x1 = line_seg(1:2) + line_seg(4)/2*[cos(line_seg(3)) sin(line_seg(3))];
%     x2 = line_seg(1:2) - line_seg(4)/2*[cos(line_seg(3)) sin(line_seg(3))];
%     plot([x1(1) x2(1)], [x1(2) x2(2)], 'r-', 'linewidth', 2); % updated
                  
end

