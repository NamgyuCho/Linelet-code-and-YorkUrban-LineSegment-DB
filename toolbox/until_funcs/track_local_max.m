function [ ls, pixel_visited ] = track_local_max( im_grad, init_pos, init_dir )
%TRACK_LOCAL_MAX Summary of this function goes here
%   Detailed explanation goes here
ls = [];
size_im = size(im_grad);
pixel_stack = init_pos;
pixel_visited = false(size_im);
bInit_dir_modified = 0;

% one direction
while ~isempty(pixel_stack)
    % Stack pop
    cur_pixel = pixel_stack(1,:);
    pixel_stack(1,:) = [];
    
    pixel_visited(cur_pixel(2), cur_pixel(1)) = true;
    ls = [ls; cur_pixel];        
    
    next_loc = repmat(cur_pixel, size(init_dir,1), 1) + init_dir;
    ind = next_loc(:,1) > 0 & next_loc(:,1) <= size_im(2) & next_loc(:,2) > 0 & next_loc(:,2) <= size_im(1);        
    next_loc = next_loc(ind,:);    
    
    next_grad = im_grad(sub2ind(size_im, next_loc(:,2), next_loc(:,1)));
    [~, ind] = max(next_grad);
    
    if ~isempty(ind) &&...
            (next_loc(ind,1) > 0 && next_loc(ind,1) <= size_im(2)) && ...
            (next_loc(ind,2) > 0 && next_loc(ind,2) <= size_im(1)) && ...
            ~pixel_visited(next_loc(ind,2), next_loc(ind,1)) &&...
            im_grad(next_loc(ind,2), next_loc(ind,1))
        pixel_stack = [pixel_stack; next_loc(ind,:)];
        
        
        % Remove one of marginal tracking directions if either of them are selected
%         if ~bInit_dir_modified
%             moving_dir = next_loc(ind,:) - cur_pixel;            
%             [~, indA] = intersect(init_dir, moving_dir, 'Rows');            
%             if indA ~= 2
%                 init_dir = init_dir([indA 2],:); 
%                 bInit_dir_modified = 1;
%             end
%         end
    end
end
