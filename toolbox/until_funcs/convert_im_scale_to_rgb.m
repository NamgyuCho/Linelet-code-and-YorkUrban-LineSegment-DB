function [ im_rgb ] = convert_im_scale_to_rgb( im_sc, cmap )
%CONVERT_IM_SCALE_TO_RGB Summary of this function goes here
%   Detailed explanation goes here

im_r = zeros(size(im_sc));
im_g = zeros(size(im_sc));
im_b = zeros(size(im_sc));

val_sc = unique(im_sc);

for i = 1:max(im_sc(:))
    ind = im_sc == i;
    
    im_r(ind) = uint8(cmap(i,1) * 255);
    im_g(ind) = uint8(cmap(i,2) * 255);
    im_b(ind) = uint8(cmap(i,3) * 255);
end

im_rgb = cat(3, im_r, im_g, im_b);
% im_rgb = uint8(im_rgb);

end

