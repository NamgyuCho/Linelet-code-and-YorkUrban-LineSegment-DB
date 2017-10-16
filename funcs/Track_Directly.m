function [px_visit_map, line_segment] = Track_Directly(im_gray, param)
size_im = size(im_gray);
[im_grad, im_dir, im_dird] = get_gradient(im_gray, 'Sobel');
im_grad(:,1) = 0; im_grad(:,end) = 0; im_grad(1,:) = 0; im_grad(end,:) = 0;
im_dird(im_dird < 0) = im_dird(im_dird < 0) + 180;
im_nms = nonmaxsup(im_grad, im_dird, 1.2, false);

im_grad = im_grad .* single(im_nms > 0);

% ind = im_grad <= param.thres_1st_grad_mag;
% im_grad(ind) = 0;
% im_dir(ind) = 0;

% Remove less significant pixels wrt the gradient magnitue
[pixel_grad_val, ind] = sort(im_grad(:), 'descend');
[ind_r, ind_c] = ind2sub(size_im, ind);
px_visit_map = single(zeros(size_im));

line_segment = [];
pixel_used = false(size_im);
num_lsc = 1;

[im_gx, im_gy] = imgradientxy(im_gray);

if strcmp(param.lsc.method, 'LOCAL')
    for i = 1:length(ind_r)
        try
        if pixel_used(ind_r(i), ind_c(i)), continue; end
        
        [xx, yy] = get_patch_indices(ind_c(i), ind_r(i), param.lsc.radius, size_im);
                
        tmp = im_grad(yy, xx);  [ ~, ind_max, max_center ] = get_matrix_extreme_val( tmp, 'max' );
        tmp(ind_max) = 0;        
        
        if ind_c(i) == 172 && ind_r(i) == 95
            max_center;
        end
        
        [im_lsc_map, pixel_used_local] = cluster_pixel_gradient_dir([max_center(2) max_center(1)], im_grad(yy,xx), im_dir(yy,xx), pixel_used(yy,xx), num_lsc, param);
        if isempty(im_lsc_map), continue; end
        pixel_used(yy,xx) = pixel_used(yy,xx) | pixel_used_local;
        im_grad_local = im_grad(yy,xx) .* single(im_lsc_map > 0);
        
        % The first direction to near local maximum        
        [ ~, ind_max, sub_max ] = get_matrix_extreme_val( tmp, 'max' ); tmp(ind_max) = 0;
        init_dir = get_move_offset(atan2(sub_max(1) - max_center(1), sub_max(2) - max_center(2))); % (x,y) order
        [local_ls, pixel_visited1] = track_local_max(im_grad_local, [max_center(2) max_center(1)], init_dir);
        px_visit_map(yy,xx) = px_visit_map(yy,xx) | pixel_visited1;
        
        % The next direction to near local maximum        
        [ ~, ind_max, sub_max ] = get_matrix_extreme_val( tmp, 'max' );
        init_dir = get_move_offset(atan2(sub_max(1) - max_center(1), sub_max(2) - max_center(2))); % (x,y) order
        [local_ls, pixel_visited2] = track_local_max(im_grad_local, [max_center(2) max_center(1)], init_dir);
        px_visit_map(yy,xx) = px_visit_map(yy,xx) | pixel_visited2;
        catch err
            fprintf('error occured at %d, ', i);
            rethrow(err);
        end
    end
end
