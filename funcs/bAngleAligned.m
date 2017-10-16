function [ bAligned ] = bAngleAligned( ang1, ang2, thres )
%BANGLEALIGNED Summary of this function goes here
%   Detailed explanation goes here
% ang1 = line_gnd(i_gnd,8); ang2 = line_est(idx_cand,8);
% ang1 = ang_pivot; ang2 = ang_pt-pi/2;
% thres = param.thres_angle_diff;
%   /* it is assumed that 'theta' and 'a' are in the range [-pi,pi] */
%   theta -= a;
%   if( theta < 0.0 ) theta = -theta;
%   if( theta > M_3_2_PI )
%     {
%       theta -= M_2__PI;
%       if( theta < 0.0 ) theta = -theta;
%     }
% 
%   return theta < prec;

n1 = numel(ang1);
n2 = numel(ang2); % assume that n2 >= n1 and n1 == 1;

if n1 == 1 && n2 >= n1
    % More newer one
    idx = ang1 < 0; ang1(idx) = ang1(idx) + pi;
    ang1 = mod(ang1, pi);
    
    idx = ang2 < 0; ang2(idx) = ang2(idx) + pi;
    ang2 = mod(ang2, pi);
        
    ang2 = ang2 - ang1;
    idx = ang2 < 0;         ang2(idx) = -ang2(idx);
    idx = ang2 > pi*3/4;    ang2(idx) = pi - ang2(idx);
    
%     % Newer version
%     ang2 = ang2 - ang1;
%     idx = ang2 < 0;         ang2(idx) = -ang2(idx);
%     idx = ang2 > pi*3/2;    ang2(idx) = ang2(idx) - pi*2;
%     idx = ang2 < 0;         ang2(idx) = -ang2(idx);
    
    bAligned = ang2 <= thres;
elseif n1 == n2
    % More newer one
    idx = ang1 < 0; ang1(idx) = ang1(idx) + pi;
    ang1 = mod(ang1, pi);
    
    idx = ang2 < 0; ang2(idx) = ang2(idx) + pi;
    ang2 = mod(ang2, pi);
    
    ang2 = ang2 - ang1;
    idx = ang2 < 0;         ang2(idx) = -ang2(idx);
    idx = ang2 > pi*3/4;    ang2(idx) = pi - ang2(idx);
    
%     ang2 = ang2 - ang1;
%     idx = ang2 < 0;         ang2(idx) = -ang2(idx);
%     idx = ang2 > pi*3/2;    ang2(idx) = ang2(idx) - pi*2;
%     idx = ang2 < 0;         ang2(idx) = -ang2(idx);
    
    bAligned = ang2 <= thres;
else
    error('The length of angle 1 should be one or both length should be equal.');
end



% Previous methods
%
% ang2 = ang2 - ang1;
% 
% idx = ang2 < 0.0;
% ang2(idx) = -ang2(idx);
% 
% idx = ang2 > pi*3/2;
% ang2(idx) = ang2(idx) - pi*2;
% 
% idx = ang2 < 0.0;
% ang2(idx) = -ang2(idx);
% 
% bAligned = ang2 < thres;
    

%     #define M_3_2_PI 4.71238898038
% 
% /** 2 pi */
% #define M_2__PI  6.28318530718
% int x, int y, image_double angles, double theta,
%                       double prec

%   double a;
% 
%   /* check parameters */
%   if( angles == NULL || angles->data == NULL )
%     error("isaligned: invalid image 'angles'.");
%   if( x < 0 || y < 0 || x >= (int) angles->xsize || y >= (int) angles->ysize )
%     error("isaligned: (x,y) out of the image.");
%   if( prec < 0.0 ) error("isaligned: 'prec' must be positive.");
% 
%   /* angle at pixel (x,y) */
%   a = angles->data[ x + y * angles->xsize ];
% 
%   /* pixels whose level-line angle is not defined
%      are considered as NON-aligned */
%   if( a == NOTDEF ) return FALSE;  /* there is no need to call the function
%                                       'double_equal' here because there is
%                                       no risk of problems related to the
%                                       comparison doubles, we are only
%                                       interested in the exact NOTDEF value */
% 
%   /* it is assumed that 'theta' and 'a' are in the range [-pi,pi] */
%   theta -= a;
%   if( theta < 0.0 ) theta = -theta;
%   if( theta > M_3_2_PI )
%     {
%       theta -= M_2__PI;
%       if( theta < 0.0 ) theta = -theta;
%     }
% 
%   return theta < prec;

end

