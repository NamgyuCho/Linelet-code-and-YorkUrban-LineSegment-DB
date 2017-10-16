function [x, y, c] = Generate_AAB_Proposal(x1, x2, y1, y2, pt1, pt2)
% pt_sam1(1), pt_sam2(1), pt_sam1(2), pt_sam2(2), pt1(1), pt2(1)
% x1 = pt_sam1(1); x2 = pt_sam2(1);
% y1 = pt_sam1(2); y2 = pt_sam2(2);
% x10 = pt1(1); x20 = pt2(1);

dx = x2 - x1;
dy = y2 - y1;
%preallocate memory for x,y, and c
x=zeros(floor(sqrt(dx^2+dy^2)),2);
y=zeros(size(x));
c=zeros(size(x));

swapped=false;
if abs(dx) < abs(dy)
    [x1 y1]=swap (x1, y1);
    [x2 y2]=swap (x2, y2);
    [dx dy]=swap (dx, dy);    
    swapped=true;
end
if x2 < x1
    [x1 x2]=swap (x1, x2);
    [y1 y2]=swap (y1, y2);
end

if pt2(1) < pt1(1) 
    [pt1, pt2] = swap_vec(pt1, pt2);
end

dx_l = abs(x1 - pt1(1));
dx_r = abs(pt2(1) - x2);

gradient = dy / dx;

% handle first endpoint
xend = round(x1);
yend = y1 + gradient * (xend - x1);
xgap = rfpart(x1 + 0.5);
xpxl1 = xend;  % this will be used in the main loop
ypxl1 = ipart(yend);
% x(1)=xpxl1; y(1)=ypxl1; c(1)=rfpart(yend) * xgap;
% x(2)=xpxl1; y(2)=ypxl1 + 1; c(2)=fpart(yend) * xgap;
intery = yend + gradient; % first y-intersection for the main loop

% handle second endpoint
xend = round (x2);
yend = y2 + gradient * (xend - x2);
xgap = fpart(x2 + 0.5);
xpxl2 = xend;  % this will be used in the main loop
ypxl2 = ipart (yend);
% x(3)=xpxl2; y(3)=ypxl2; c(3)=rfpart (yend) * xgap;
% x(4)=xpxl2; y(4)=ypxl2 + 1; c(4)=fpart (yend) * xgap;

% main loop
k=1;
for i =xpxl1:xpxl2
    x(k,1)=i; y(k,1)=ipart (intery); c(k,1)=rfpart (intery);
%     k=k+1;
    x(k,2)=i; y(k,2)=ipart (intery) + 1; c(k,2)= fpart (intery);
    intery = intery + gradient;
    k=k+1;
end

%truncate the vectors to proper sizes
x=x(1:k-1,:);
y=y(1:k-1,:);
c=c(1:k-1,:);

% Append additional points to the left and right end points
if dx_l > 0
    xl_add = zeros(dx_l, 2);
    yl_add = xl_add;
    cl_add = xl_add;
    intery = y1 - gradient;
    k = 1;
    for i = xpxl1-1:-1:round(pt1(1))
        xl_add(dx_l - k + 1, 1) = i;
        yl_add(dx_l - k + 1, 1) = ipart(intery);
        cl_add(dx_l - k + 1, 1) = rfpart(intery);
        xl_add(dx_l - k + 1, 2) = i;
        yl_add(dx_l - k + 1, 2) = ipart(intery) + 1;
        cl_add(dx_l - k + 1, 2) = fpart(intery);
        intery = intery - gradient;
        k = k + 1;
    end
    xl_add = xl_add(1:k-1,:);
    yl_add = yl_add(1:k-1,:);
    cl_add = cl_add(1:k-1,:);
    
    x = [xl_add; x];
    y = [yl_add; y];
    c = [cl_add; c];
end

if dx_r > 0
    xr_add = zeros(dx_r, 2);
    yr_add = xr_add;
    cr_add = xr_add;
    intery = ypxl2 + gradient;
    k = 1;
    for i = xpxl2+1:round(pt2(1))
        xr_add(k,1) = i;
        yr_add(k,1) = ipart(intery);
        cr_add(k,1) = rfpart(intery);
        xr_add(k, 2) = i;
        yr_add(k, 2) = ipart(intery) + 1;
        cr_add(k, 2) = fpart(intery);
        intery = intery + gradient;
        k = k + 1;
    end
    
    xr_add = xr_add(1:k-1,:);
    yr_add = yr_add(1:k-1,:);
    cr_add = cr_add(1:k-1,:);
    
    x = [x; xr_add];
    y = [y; yr_add];
    c = [c; cr_add];
end

if swapped
    [x y]=swap (x, y);
end


%integer part
function i=ipart(x)
if x>0
    i=floor(x);
else
    i=ceil(x);
end

function r=round(x)
r= ipart(x + 0.5);

%fractional part
function f=fpart(x)
f=x-ipart(x);

function rf=rfpart(x)
rf= 1 - fpart(x);

function [x y]=swap(x,y)
a=x;
x=y;
y=a;

