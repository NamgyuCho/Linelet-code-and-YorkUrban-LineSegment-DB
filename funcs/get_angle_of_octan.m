function [ ret_ang ] = get_angle_of_octan( octan_number, ang_val )
%GET_ANGLE_OF_OCTAN Summary of this function goes here
%   Detailed explanation goes here

if octan_number > 0
    ang_offset = 45 * (octan_number - 1);
    ret_ang = ang_val + ang_offset;
else
    ang_offset = -45 * (octan_number + 1);
    ret_ang = -ang_val + ang_offset;
end

ret_ang = ret_ang * pi / 180;

end

