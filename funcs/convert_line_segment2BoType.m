function [ ls_dst ] = convert_line_segment2BoType( ls_src )
    %CONVERT_LINE_SEGMENT2BOTYPE Summary of this function goes here
    %   Detailed explanation goes here
    % ls_src = line_gnd;
    % Convert to the required line type
    
    % e.g., gnd is in a form, (x1, y1, x2, y2, center_x, center_y, length, angle)
    ls_dst = [];
    
    mean_x = ls_src(:,5);
    mean_y = ls_src(:,6);
    theta = ls_src(:,8);
    l = ls_src(:, 7);
    
    x1 = ls_src(:,1);
    x2 = ls_src(:,3);
    y1 = ls_src(:,2);
    y2 = ls_src(:,4);
    
    % theta = theta - pi / 2;
    
    r = mean_x.*cos(theta)+mean_y.*sin(theta);
    
    ls_dst = [y1 y2 x1 x2 theta r l mean_y mean_x ];
end