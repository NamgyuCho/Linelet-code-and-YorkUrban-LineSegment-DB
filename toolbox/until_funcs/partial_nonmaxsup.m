function [ im_nms_x, im_nms_y, mem_x, mem_y ] = partial_nonmaxsup( im_gx, im_gy, bSkeletonize, nms_radius )
%PARTIAL_NONMAXSUP Summary of this function goes here
%   [ im_nms_x, im_nms_y, mem_x, mem_y ] = partial_nonmaxsup( im_gx, im_gy, bSkeletonize, nms_radius )

size_im = size(im_gx);

im_nms_x = im_gx;
im_nms_y = im_gy;
mem_x = zeros(size_im);
mem_y = zeros(size_im);


for i=1+nms_radius:size_im(1)-nms_radius
    for j=1+nms_radius:size_im(2)-nms_radius
        % x direction
        if im_gx(i,j)
            x_ind = [j-nms_radius:j+nms_radius]; x_ind(nms_radius+1) = [];
            if any(im_nms_x(i,j) <= im_gx(i, x_ind))
                im_nms_x(i,j) = 0; 
                mem_x(i,j) = 0;
            else
                mem_x(i,j) = 1;
            end

            if any(im_nms_x(i,j) == im_gx(i, x_ind))
                mem_x(i,j) = .5; 
            end
        end
        
        % y direction
        if im_gy(i,j)
            y_ind = [i-nms_radius:i+nms_radius]; y_ind(nms_radius+1) = [];
            if any(im_nms_y(i,j) <= im_gy(y_ind, j))
                im_nms_y(i,j) = 0;
                mem_y(i,j) = 0;
            else
                mem_y(i,j) = 1;
            end
            
            if any(im_nms_y(i,j) == im_gy(y_ind, j))
                mem_y(i,j) = .5;
            end
        end
    end
end
% 
nms_radius = 1;
im_nms_x(1, :) = 0; im_nms_x(end, :) = 0;  
im_nms_x(:, 1) = 0; im_nms_x(:, end) = 0;  
im_nms_y(1, :) = 0; im_nms_y(end, :) = 0;  
im_nms_y(:, 1) = 0; im_nms_y(:, end) = 0;  

skel_x = bwmorph(im_nms_x>0,'thin');
skel_y = bwmorph(im_nms_y>0,'thin');

im_clean_x = bwmorph(im_nms_x > 0, 'clean');
im_clean_y = bwmorph(im_nms_y > 0, 'clean');

im_nms_x = im_nms_x .* double(im_clean_x);
im_nms_y = im_nms_y .* double(im_clean_y);

if bSkeletonize
    im_nms_x = im_nms_x.*skel_x; 
    im_nms_y = im_nms_y.*skel_x; 
end