function [valid_pix_ind, valid_pix_ind_main] = remove_boundary_violation_pix(size_im, sp_pix, num_set)

num_pixel = size(sp_pix,1) / num_set;
valid_pix_ind = [];

ind = [];
for i = 1:num_set
    idx = [1:num_pixel] + (i-1) * num_pixel;
    
    ind{i,1} = find(sp_pix(idx,1) >= 1 & sp_pix(idx,1) <= size_im(2) &...
                    sp_pix(idx,2) >= 1 & sp_pix(idx,2) <= size_im(1));
end
    
ind = intersect(intersect(ind{1}, ind{2}), ind{3});
valid_pix_ind_main = ind;

for i = 1:num_set
    idx = [1:num_pixel] + (i-1) * num_pixel;
    valid_pix_ind = [valid_pix_ind; idx(ind)'];
end