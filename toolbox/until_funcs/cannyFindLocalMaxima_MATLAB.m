function [ idxLocalMax ] = cannyFindLocalMaxima_MATLAB( direction,ix,iy,mag )
%CANNYFINDLOCALMAXIMA_MATLAB Summary of this function goes here
%   Detailed explanation goes here
 
%
% This sub-function helps with the non-maximum suppression in the Canny
% edge detector.  The input parameters are:
%
%   direction - the index of which direction the gradient is pointing,
%               read from the diagram below. direction is 1, 2, 3, or 4.
%   ix        - input image filtered by derivative of gaussian along x
%   iy        - input image filtered by derivative of gaussian along y
%   mag       - the gradient magnitude image
%
%    there are 4 cases:
%
%                         The X marks the pixel in question, and each
%         3     2         of the quadrants for the gradient vector
%       O----0----0       fall into two cases, divided by the 45
%     4 |         | 1     degree line.  In one case the gradient
%       |         |       vector is more horizontal, and in the other
%       O    X    O       it is more vertical.  There are eight
%       |         |       divisions, but for the non-maximum suppression
%    (1)|         |(4)    we are only worried about 4 of them since we
%       O----O----O       use symmetric points about the center pixel.
%        (2)   (3)


[m,n] = size(mag);

% Find the indices of all points whose gradient (specified by the
% vector (ix,iy)) is going in the direction we're looking at.

switch direction
    case 1
        idx = find((iy<=0 & ix>-iy)  | (iy>=0 & ix<-iy));
    case 2
        idx = find((ix>0 & -iy>=ix)  | (ix<0 & -iy<=ix));
    case 3
        idx = find((ix<=0 & ix>iy) | (ix>=0 & ix<iy));
    case 4
        idx = find((iy<0 & ix<=iy) | (iy>0 & ix>=iy));
end

% Exclude the exterior pixels
if ~isempty(idx)
    v = mod(idx,m);
    extIdx = (v==1 | v==0 | idx<=m | (idx>(n-1)*m));
    idx(extIdx) = [];
end

ixv = ix(idx);
iyv = iy(idx);
gradmag = mag(idx);

% Do the linear interpolations for the interior pixels
switch direction
    case 1
        d = abs(iyv./ixv);
        gradmag1 = mag(idx+m).*(1-d) + mag(idx+m-1).*d;
        gradmag2 = mag(idx-m).*(1-d) + mag(idx-m+1).*d;
    case 2
        d = abs(ixv./iyv);
        gradmag1 = mag(idx-1).*(1-d) + mag(idx+m-1).*d;
        gradmag2 = mag(idx+1).*(1-d) + mag(idx-m+1).*d;
    case 3
        d = abs(ixv./iyv);
        gradmag1 = mag(idx-1).*(1-d) + mag(idx-m-1).*d;
        gradmag2 = mag(idx+1).*(1-d) + mag(idx+m+1).*d;
    case 4
        d = abs(iyv./ixv);
        gradmag1 = mag(idx-m).*(1-d) + mag(idx-m-1).*d;
        gradmag2 = mag(idx+m).*(1-d) + mag(idx+m+1).*d;
end
idxLocalMax = idx(gradmag>=gradmag1 & gradmag>=gradmag2);