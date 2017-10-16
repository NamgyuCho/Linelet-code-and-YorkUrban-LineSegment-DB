function [ out_ang ] = get_opposite_angle( in_ang, CordType )
%GET_OPPOSITE_ANGLE Summary of this function goes here
%   Detailed explanation goes here

if strcmp(CordType, 'posneg')
    if in_ang > 0 && in_ang < pi
        out_ang = in_ang - pi;
    elseif in_ang < 0 && in_ang > -pi
        out_ang = in_ang + pi;
    elseif in_ang == 0
        out_ang = pi;
    elseif in_ang == -pi || in_ang == pi
        out_ang = 0;
    end
elseif strcmp(CordType, 'posonly')
    out_ang = mod(in_ang + pi, pi*2);
end