function [xx, yy] = get_patch_indices(x0, y0, rad, size_im)

xx = x0 - rad : x0 + rad;
if x0 - rad < 1,              xx = 1:x0+rad;
elseif x0 + rad > size_im(2), xx = x0 - rad : size_im(2);
end

yy = y0 - rad : y0 + rad;
if y0 - rad < 1,              yy = 1:y0+rad;
elseif y0 + rad > size_im(1), yy = y0 - rad : size_im(1);
end

