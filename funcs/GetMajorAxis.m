function [ iMajor, iMinor ] = GetMajorAxis( ll_type )
%GETAXISINDEX Summary of this function goes here
%   Return indices of axis, whether x- or -y, according to the type of linelet, horizontal, vertical, or diagonal

if      ll_type == 1 % horizontal
    iMajor = [1 3];
    iMinor = [2 4];
elseif  ll_type == 3 % vertical
    iMajor = [2 4];
    iMinor = [1 3];
else                 % diagonal
    iMajor = [1 3];
    iMinor = [2 4];
end

end

