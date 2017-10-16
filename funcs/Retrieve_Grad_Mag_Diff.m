function [grad_mag_diff, grad_mag_std, grad_dir_std, overlap_ratio] = Retrieve_Grad_Mag_Diff(im_grad, im_dir, lc, xys, thetas, lens)

if size(xys,1) ~= size(thetas,1) || size(xys,1) ~= size(lens,1) || size(thetas,1) ~= size(lens,1)
    fprintf('Length of input parameter must be same!\n');
end
thetas = thetas - pi;
num_input = size(xys, 1);
grad_mag_diff = [];
grad_mag_std = [];
grad_dir_std = []; 
overlap_ratio = [];

im_size = size(im_grad);

for i = 1:num_input
    % Make end points
    p1 = floor(xys(i,:) + lens(i)*[cos(thetas(i)) sin(thetas(i))]);
    p2 = floor(xys(i,:) - lens(i)*[cos(thetas(i)) sin(thetas(i))]);
    [xx, yy] = bresenham(p1(1), p1(2), p2(1), p2(2));
    
    if thetas(i) > pi/4 && thetas(i) <= pi*3/4
        xx_p = xx - 1;
        yy_p = yy;
        xx_n = xx + 1 ;        
        yy_n = yy;
    elseif thetas(i) > pi*5/4 && thetas(i) <= pi*7/4
        xx_p = xx - 1;
        yy_p = yy;
        xx_n = xx + 1;        
        yy_n = yy;
    else
        xx_p = xx;
        yy_p = yy - 1;
        xx_n = xx;        
        yy_n = yy + 1;
    end
    
    % remove pixel violate image region
    ind1 = find(xx_p < 1 | xx_p >= im_size(2) | xx_n < 1 | xx_n >= im_size(2) |...
                yy_p < 1 | yy_p >= im_size(1) | yy_n < 1 | yy_n >= im_size(1));
            
    xx(ind1) = [];  xx_p(ind1) = [];  xx_n(ind1) = [];
    yy(ind1) = [];  yy_p(ind1) = [];  yy_n(ind1) = [];
    
    if isempty(xx)
        continue;
    end
    
    indO = sub2ind(im_size, yy, xx);
    indP = sub2ind(im_size, yy_p, xx_p);
    indN = sub2ind(im_size, yy_n, xx_n);
        
    C = intersect(lc, [xx yy], 'rows');
    
    grad_mag_diff = [grad_mag_diff; sum( abs(2*im_grad(indO) - im_grad(indP) - im_grad(indN)) )];
    grad_mag_std  = [grad_mag_std; std(im_grad(indO))];
    grad_dir_std  = [grad_dir_std; std(im_dir(indO))];
    overlap_ratio = [overlap_ratio; length(C) / size(xx, 1)];
end


% imtmp = 0 * im_grad;
% imtmp(indO) = 10;
% imtmp(indP) = 20;
% imtmp(indN) = 30;
% figure; imagesc(imtmp);