function [a, b, c, nom] = get_line_param(line_segment)

a = tan(line_segment(5));
b = -1;
c = -a*line_segment(1) + line_segment(2);
nom = sqrt(a^2 + b^2);