function [ ls_map, line_segment ] = Direct_LS_Track( im_grad, im_dir, im_nms, param )
%DIRECT_LS_TRACK Summary of this function goes here
%   Detailed explanation goes here
tic
size_im = size(im_grad);
pixel_used = zeros(size_im);
ls_map = pixel_used;

tmp = im_grad < 100;
im_grad(tmp) = 0;

% Remove less significant pixels wrt the gradient magnitue
[pixel_grad_val, ind] = sort(im_grad(:), 'descend');
[ind_r, ind_c] = ind2sub(size_im, ind);
angle_val = im_dir(ind);
ls_map = zeros(size_im);

tmp = pixel_grad_val > 0;
ind_r = ind_r(tmp);
ind_c = ind_c(tmp);

line_segment = [];
num_lsc = 1;

nb_offset = [1 0; 1 1; 0 1; -1 1; -1 0; -1 -1; -1 0; -1 1];
ind_lc = 1;
    
for i = 1:length(ind_r)
    if pixel_used(ind_r(i), ind_c(i)), continue; end
    
    % Find the locally maximum pixel from neighbors
    lc = [ind_c(i), ind_r(i)];
    
    nb_ind = nb_offset + repmat(lc, 8, 1);
    [ pts, valid_idx ] = return_valid_loc( nb_ind, size_im, [0 0]);
    [max_val, max_ind] = max(im_grad(sub2ind(size_im, pts(:,2), pts(:,1))));
    nb_dir = pts(max_ind,:) - lc;
    
    if nb_dir == [1 0],       search_offset = [1 -1; 1 0; 1 1];
    elseif nb_dir == [1 1],   search_offset = [1 0; 1 1; 0 1];
    elseif nb_dir == [0 1],   search_offset = [1 1; 0 1; -1 1];
    elseif nb_dir == [-1 1],  search_offset = [0 1; -1 1; -1 0];
    elseif nb_dir == [-1 0],  search_offset = [-1 1; -1 0; -1 -1];
    elseif nb_dir == [-1 -1], search_offset = [-1 0; -1 -1; 0 -1];
    elseif nb_dir == [0 -1],  search_offset = [-1 -1; 0 -1;1 -1];
    elseif nb_dir == [1 -1],  search_offset = [0 -1;1 -1; 1 0];
    end
    
    pixel_stack = lc;
    pixel_used(ind_r(i), ind_c(i)) = ind_lc;
    
    reg_angle = angle_val(i);
    sumdx = cos( reg_angle );
    sumdy = sin( reg_angle );
    
    % one direction
    while ~isempty(pixel_stack)
        cur_pixel = pixel_stack(1,:);
        pixel_stack(1,:) = [];
        pixel_used(cur_pixel(2), cur_pixel(1)) = ind_lc;
        ls_map(cur_pixel(2), cur_pixel(1)) = ind_lc;
        
        next_loc = repmat(cur_pixel, size(search_offset,1), 1) + search_offset;
        ind = next_loc(:,1) > 0 & next_loc(:,1) <= size_im(2) & next_loc(:,2) > 0 & next_loc(:,2) <= size_im(1);
        next_loc = next_loc(ind,:);
        nb_pix = next_loc;
        
        if isempty(next_loc), continue; end
        
        next_grad = im_grad(sub2ind(size_im, next_loc(:,2), next_loc(:,1)));
        [next_grad, ind] = max(next_grad);
        
        next_loc = next_loc(ind,:);
        next_ang = im_dir(next_loc(2), next_loc(1));
        
        nb_dir2 = next_loc - cur_pixel;
        inter_vec = intersect(nb_dir2, search_offset, 'rows');
        
        if ~isempty(ind) && ~isempty(inter_vec) &&...
                (next_loc(1) > 0 && next_loc(1) <= size_im(2)) && ...
                (next_loc(2) > 0 && next_loc(2) <= size_im(1)) && ...
                ~pixel_used(next_loc(2), next_loc(1)) && next_grad &&...
                angle_alignment(reg_angle, next_ang, param.thres_angle_diff + 45/180*pi) &&...
                next_grad >= param.thres_1st_grad_mag
    
            sumdx = sumdx + cos(next_ang);
            sumdy = sumdy + sin(next_ang);
            reg_angle = atan2(sumdy, sumdx); % Updating angle performs better slightly, not significantly                    
            pixel_stack = [pixel_stack; next_loc];
            pixel_used(sub2ind(size_im, nb_pix(:,2), nb_pix(:,1))) = ind_lc;
        end
    end
    
    ind_lc = ind_lc + 1;
end

toc
figure; imshowpair(im_gray, ls_map > 0, 'blend')
figure; imagesc(ls_map)
end