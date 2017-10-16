function [ lsc_map_ret ] = GetWholeLsc( lsc_map, lsc_org_map, size_im )
%GETWHOLELSC Summary of this function goes here
%   Detailed explanation goes here

unq_x = unique(lsc_map);
lsc_map_ret = zeros(size_im);
for i = 1 :length(unq_x)
    lsc_map_ret(lsc_org_map == unq_x(i)) = unq_x(i);
end

end

