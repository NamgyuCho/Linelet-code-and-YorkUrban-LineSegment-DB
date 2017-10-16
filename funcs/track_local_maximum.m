function [ pixel_visited, cur_pixel, track_pattern ] = track_local_maximum( im_target, pixel_visited, size_im, init_pos, offsets, Bresenham_bin )
%TRACK_LOCAL_MAXIMUM Summary of this function goes here
%   Detailed explanation goes here

pixel_stack = init_pos;
prev_pos = init_pos;
curr_pos = init_pos;
track_pattern = [];
x_offset = 0;

% one direction
while ~isempty(pixel_stack)
    cur_pixel = pixel_stack(1,:);
    pixel_stack(1,:) = [];
    pixel_visited(cur_pixel(2), cur_pixel(1)) = true;
    
    next_loc = repmat(cur_pixel, size(offsets,1), 1) + offsets;
    ind = next_loc(:,1) > 0 & next_loc(:,1) <= size_im(2) & next_loc(:,2) > 0 & next_loc(:,2) <= size_im(1);        
    next_loc = next_loc(ind,:);    
    
    next_grad = im_target(sub2ind(size_im, next_loc(:,2), next_loc(:,1)));
    [~, ind] = max(next_grad);
    
    if ~isempty(ind) &&...
            (next_loc(ind,1) > 0 && next_loc(ind,1) <= size_im(2)) && ...
            (next_loc(ind,2) > 0 && next_loc(ind,2) <= size_im(1)) && ...
            ~pixel_visited(next_loc(ind,2), next_loc(ind,1)) && im_target(next_loc(ind,2), next_loc(ind,1))
        pixel_stack = [pixel_stack; next_loc(ind,:)];
        %             im_target(sub2ind(size_im, next_loc(:,2), next_loc(:,1)) = 0;
        
        curr_pos = next_loc(1,:);
        x_offset = x_offset + 1;
        
        if curr_pos(2) ~= prev_pos(2)
            track_pattern = [track_pattern x_offset];
            x_offset = 0;
        end
    end
end
