function [offset, sub_offset, bX_major_oct, octan_section_num] = get_octan_offset_and_target_gradient(init_angle)


if init_angle >= 0 && init_angle < pi/4
    offset = [1 0; 1 1];
    sub_offset = [0 1];
    octan_section_num = 1;
elseif init_angle >= pi/4 && init_angle < pi/2
    offset = [0 1; 1 1];
    sub_offset = [1 0];
    octan_section_num = 2;
elseif init_angle >= pi/2 && init_angle < pi*3/4
    offset = [0 1; -1 1];    
    sub_offset = [1 0];
    octan_section_num = 3;
elseif init_angle >= pi*3/4 && init_angle <= pi
    offset = [-1 0; -1 1];    
    sub_offset = [0 1];
    octan_section_num = 4;
elseif init_angle > -pi/4 && init_angle <= 0
    offset = [1 0; 1 -1];
    sub_offset = [0 -1];
    octan_section_num = -1;
elseif init_angle > -pi/2 && init_angle <= -pi/4
    offset = [0 -1; 1 -1];
    sub_offset = [1 0];
    octan_section_num = -2;
elseif init_angle > -pi*3/4 && init_angle <= -pi/2
    offset = [0 -1; -1 -1];
    sub_offset = [1 0];
    octan_section_num = -3;
elseif init_angle > -pi && init_angle <= -pi*3/4
    offset = [-1 0; -1 -1];
    sub_offset = [0 -1];
    octan_section_num = -4;
elseif init_angle <= -pi
    offset = [-1 0; -1 -1];
    sub_offset = [0 1];
    octan_section_num = -4;
end


if offset(1,2) == 0
    bX_major_oct = false;
else
    bX_major_oct = true;
end