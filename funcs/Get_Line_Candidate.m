function [region_candi_map, im_lsc_map, lc_st, line_segment, pixel_cluster] = Get_Line_Candidate(...
    im_gray, im_grad, im_nms, im_dir, im_dird, param, bEstimateLSc, bDrawRet)
size_im = size(im_grad);

% Remove less significant pixels wrt the gradient magnitue
[pixel_grad_val, ind] = sort(im_grad(:), 'descend');
[ind_r, ind_c] = ind2sub(size_im, ind);
angle_val = im_dir(ind);
im_lsc_map = zeros(size_im);

line_segment = [];
lc_st = [];
pixel_cluster = [];
pixel_used = zeros(size_im);
num_lsc = 1;

x_offset = -1:1;
y_offset = -1:1;

for i = 1:length(ind_r)
    if pixel_used(ind_r(i), ind_c(i)), continue; end
    
    ind_lc = 1;
    lc = [ind_c(i), ind_r(i)];
    weight = im_grad(ind_r(i), ind_c(i));
    reg_angle = angle_val(i);
    sumdx = cos( reg_angle );
    sumdy = sin( reg_angle );
    
    ii = 1;
    
    if weight <= param.thres_grad_mag, continue; end
    
    % Cluster a set of pixels which have similar gradient direction values
    while ii <= ind_lc
        for ii_x = 1:length(x_offset)%xx = lc(ii, 1) + x_offset%-1:lc(ii, 1)+1
            xx = lc(ii, 1) + x_offset(ii_x);
            for ii_y = 1:length(y_offset)%   yy = lc(ii, 2) + y_offset%-1:lc(ii, 2)+1
                yy = lc(ii, 2) + y_offset(ii_y);
                if xx <= size_im(2) && xx >= 1 &&...
                        yy <= size_im(1) && yy >= 1 &&...
                        pixel_used(yy,xx) == 0 &&...
                        angle_alignment(reg_angle, im_dir(yy,xx), param.thres_angle_diff) &&...
                        im_grad(yy,xx) > param.thres_grad_mag
                    
                    pixel_used(yy, xx) = num_lsc;
                    im_lsc_map(yy, xx) = num_lsc;
                    
                    weight = [weight; im_grad(yy, xx)];
                    
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
    
    % Remove too small set -- to reduce computation
    if size(lc,1) <= 3
        im_lsc_map(sub2ind(size_im, lc(:,2), lc(:,1))) = 0;
        continue;
    end
    
    lc_st(num_lsc,1).pos = lc;
    num_lsc = num_lsc + 1;
end
region_candi_map = pixel_used;
