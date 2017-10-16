function diff_val = angle_diff(ang1, ang2)

diff_val = [];

if (size(ang1,1) ~= size(ang2,1)) || (size(ang1,2) ~= size(ang2,2))
    fprintf('Dimension of two input should be equal.\n'); 
    return;
end

diff_val = abs(ang1 - ang2);
idx = diff_val > 90;
diff_val(idx) = diff_val(idx) - 90;

end