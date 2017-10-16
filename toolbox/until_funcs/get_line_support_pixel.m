function [sp_pix_top, sp_pix_mid, sp_pix_bot] = get_line_support_pixel(line_seg, size_im)

% Get supporter direction
% if (line_segment(5) >= -pi/4 && line_segment(5) <= pi/4) || (line_segment(5) <= -pi*3/4) || (line_segment(5) >= pi*3/4)
%     xy_offset = [0 1];
% else
%     xy_offset = [1 0];
% end
xy_offset = round(3*[cos(line_seg(5)+pi/2) sin(line_seg(5)+pi/2)]);

% Get pixel localtion
pts = [line_seg(1:2) + line_seg(3)*[cos(line_seg(5)) sin(line_seg(5))] ...
       line_seg(1:2) - line_seg(3)*[cos(line_seg(5)) sin(line_seg(5))] ];
[nx, ny] = bresenham(pts(1), pts(2), pts(3), pts(4));

sp_pix = [ [nx ny] + repmat(xy_offset, length(nx), 1);
            nx ny;
           [nx ny] - repmat(xy_offset, length(nx), 1) ];
       
[valid_pix_ind_tot, valid_pix_ind_main] = remove_boundary_violation_pix(size_im, sp_pix, 3);

if isempty(valid_pix_ind_tot)
    sp_pix_top = [];
    sp_pix_mid = [];
    sp_pix_bot = [];
    return;
end

sp_pix_top = [nx ny] + repmat(xy_offset, length(nx), 1);
sp_pix_mid = [nx ny];
sp_pix_bot = [nx ny] - repmat(xy_offset, length(nx), 1);

sp_pix_top = sp_pix_top(valid_pix_ind_main,:);
sp_pix_mid = sp_pix_mid(valid_pix_ind_main,:);
sp_pix_bot = sp_pix_bot(valid_pix_ind_main,:);