function [ ret ] = GetAngleDiff( ang1, ang2 )
    %GETANGLEDIFF Summary of this function goes here
    %   Detailed explanation goes here
    idx = ang1 < 0; ang1(idx) = ang1(idx) + pi;
    ang1 = mod(ang1, pi);
    
    idx = ang2 < 0; ang2(idx) = ang2(idx) + pi;
    ang2 = mod(ang2, pi);
    
    ang2 = ang2 - ang1;
    idx = ang2 < 0;         ang2(idx) = -ang2(idx);
    idx = ang2 > pi*3/4;    ang2(idx) = pi - ang2(idx);

    ret = ang2;
end

