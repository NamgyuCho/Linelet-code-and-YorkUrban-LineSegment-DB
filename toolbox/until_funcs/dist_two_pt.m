function dist = dist_two_pt(x1, x2)
if (size(ang1,1) ~= size(ang2,1)) || (size(ang1,2) ~= size(ang2,2))
    fprintf('Dimension of two input should be equal.\n'); 
    return;
end

dist = sqrt(sum((x2-x1).^2));
end