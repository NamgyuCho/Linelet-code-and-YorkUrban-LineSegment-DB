function [im_lsc_map, pixel_used] = cluster_pixel_gradient_dir(p0, im_grad, im_dir, pixel_used, lsc_id, param)
%CLUSTER_PIXEL_GRADIENT_DIR Summary of this function goes here
%   Detailed explanation goes here

size_im = size(im_grad);
im_lsc_map = zeros(size_im);
lc = p0;
reg_angle = im_dir(lc(2), lc(1));
sumdx = cos( reg_angle );
sumdy = sin( reg_angle );

ii = 1;
ind_lc = 1;

if pixel_used(lc(2), lc(1))
    im_lsc_map = [];
    return; 
end


% Cluster a set of pixels which have similar gradient direction values
while ii <= ind_lc
    for xx = lc(ii, 1)-1:lc(ii, 1)+1
        for yy = lc(ii, 2)-1:lc(ii, 2)+1
            if xx <= size_im(2) && xx >= 1 && yy <= size_im(1) && yy >= 1 &&...
               pixel_used(yy,xx) == 0 &&...
               angle_alignment(reg_angle, im_dir(yy,xx), param.thres_angle_diff) &&...
               im_grad(yy,xx) > param.thres_1st_grad_mag
                
                pixel_used(yy, xx) = lsc_id;
                im_lsc_map(yy, xx) = lsc_id;
                
                sumdx = sumdx + cos(im_dir(yy, xx));
                sumdy = sumdy + sin(im_dir(yy, xx));
                reg_angle = atan2(sumdy, sumdx); % Updating angle performs better slightly, not significantly
                
                lc = [lc; xx, yy];
                ind_lc = ind_lc + 1;
            end
        end
    end
    ii = ii + 1;
end

% Remove too small set
if size(lc,1) < param.thres_cluster_size
    im_lsc_map = [];
end

end