function [ idx_split ] = get_line_splitting_location( ll_cand, line_seg, thres_split_dist )
    %GET_LINE_SPLITTING_LOCATION Summary of this function goes here
    %   Detailed explanation goes here
    
    idx_split = [];
      
                    
    % Calculate the point-to-line distance
%     pts = [];
%     for i = 1:size(ll_cand, 1)
%         if i == 1
%             pts = [(ll_cand(i,1):ll_cand(i,2)); repmat(ll_cand(i,3), 1, ll_cand(i,4))];
%         else
%             pts = [pts, [(ll_cand(i,1):ll_cand(i,2)); repmat(ll_cand(i,3), 1, ll_cand(i,4))]];
%         end        
%     end
    pts = [ll_cand(:,1)' + ll_cand(:,4)'/2; ll_cand(:,3)'];% [xx'; yy'];
                
    num_pts = size(pts, 2);
    
    n = [-sin(line_seg(3)), cos(line_seg(3))];
    n = n / norm(n);
    
    dist1 = n*(pts-repmat(line_seg(1:2)',1,num_pts));
    
    [val, idx] = max(abs(dist1));
    
    %figure; bar(dist1); axis image
    
    % case 1)
    if val >= thres_split_dist
        idx_split = idx;
    end
    
    % case 2) Find the location where the sign of distance changes
    [val_min, idx_min] = min(dist1);
    [val_max, idx_max] = max(dist1);
    [val_absmin, idx_absmin] = min(abs(dist1));
    
    if val_min * val_max < 0 && val >= thres_split_dist
        idx_split = [idx_split; idx_absmin];        
    end
           
    if ~isempty(idx_split), idx_split = unique(idx_split); end
    
end

