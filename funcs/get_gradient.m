function [ im_grad, im_dir, im_dird, im_gx, im_gy ] = get_gradient( im_gray, strMethod, SigmaSize )
%GET_GRADIENT Summary of this function goes here
%   Detailed explanation goes here

if find(strcmp(strMethod, {'Sobel', 'Prewitt', 'CentralDifference', 'IntermediateDifference'}))
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % IT IS NECESSARY TO CHECK THE INFLUENCE OF THE KERNEL SIZE AND SIGMA
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    myfilter = fspecial('gaussian',[3 3], SigmaSize); 
    im_gray = imfilter(im_gray, myfilter, 'replicate');
    
    [im_grad, im_dird] = imgradient(im_gray, strMethod);
    [im_gx, im_gy] = imgradientxy(im_gray, strMethod);
    im_dir = im_dird / 180 * pi;
    
    im_gx = abs(im_gx); im_gy = abs(im_gy);
    
elseif strcmp(strMethod, 'LSD')
    % Smooth image
    myfilter = fspecial('gaussian',[3 3], SigmaSize);
    im_gray = imfilter(im_gray, myfilter, 'replicate');
    %figure; imshowpair(im_gray, im_gray1, 'diff'); colormap jet
    
    % ------
    % Norm 2 computation using 2x2 pixel window:
    %  A B
    %  C D
    % and
    %  com1 = D-A,  com2 = B-C.
    % Then
    %  gx = B+D - (A+C)   horizontal difference
    %  gy = C+D - (A+B)   vertical difference
    % com1 and com2 are just to avoid 2 additions.
    
    hx1 = [-1 1; -1 1];% hx2 = hx1 * -1;
    hy1 = [-1 -1; 1 1];% hy2 = hy1 * -1;
    gx1 = conv2(double(im_gray), hx1);% gx2 = imfilter(im_gray, hx2); gx = single((gx1 + gx2)/2);
    gy1 = conv2(double(im_gray), hy1);% gy2 = imfilter(im_gray, hy2); gy = single((gy1 + gy2)/2);
    im_grad = gx1.*gx1+gy1.*gy1;
    im_grad = sqrt( im_grad / 4.0 ); % gradient norm
    im_grad = im_grad(2:end, 2:end);
    im_grad(1:end,end) = 0;
    im_grad(end,1:end) = 0;
    im_dir = atan2(gx1, -gy1);
    im_dir = im_dir(2:end, 2:end);    
    im_dird = im_dir * 180 / pi;
    
    im_gx = abs(gx1) /2;
    im_gy = abs(gy1) /2;
    
    im_gx = im_gx(2:end, 2:end);
    im_gy = im_gy(2:end, 2:end);    
    im_grad(:,1) = 0; im_grad(:,end) = 0; im_grad(1,:) = 0; im_grad(end,:) = 0;
    im_gx(:,1) = 0; im_gx(:,end) = 0; im_gx(1,:) = 0; im_gx(end,:) = 0;
    im_gy(:,1) = 0; im_gy(:,end) = 0; im_gy(1,:) = 0; im_gy(end,:) = 0;
elseif strcmp(strMethod, 'OneNB')
    % Smooth image
    myfilter = fspecial('gaussian',[3 3], 0.5);
    im_gray = imfilter(im_gray, myfilter, 'replicate');
    %figure; imshowpair(im_gray, im_gray1, 'diff'); colormap jet
    
    % ------
    % Norm 2 computation using 2x2 pixel window:
    %  A B
    %  C D
    % and
    %  com1 = D-A,  com2 = B-C.
    % Then
    %  gx = B+D - (A+C)   horizontal difference
    %  gy = C+D - (A+B)   vertical difference
    % com1 and com2 are just to avoid 2 additions.
    
    hx1 = [-1 1];% hx2 = hx1 * -1;
    hy1 = [-1; 1];% hy2 = hy1 * -1;
    
    gx1 = conv2(double(im_gray), hx1);% gx2 = imfilter(im_gray, hx2); gx = single((gx1 + gx2)/2);
    gx1 = gx1(:,1:end-1);
    gx1(:,1) = 0;
    
    gy1 = conv2(double(im_gray), hy1);% gy2 = imfilter(im_gray, hy2); gy = single((gy1 + gy2)/2);
    gy1 = gy1(1:end-1,:);
    gy1(1,:) = 0;
%     figure; imshowpair(gx1, gy1, 'montage'); 
    
    im_grad = gx1.*gx1+gy1.*gy1;
    im_grad = sqrt( im_grad / 2 ); % gradient norm
    %im_grad = im_grad(2:end, 2:end);
%     im_grad(1:end,1) = 0;
%     im_grad(1,1:end) = 0;
    im_dir = atan2(gx1, -gy1);
    %im_dir = im_dir(2:end, 2:end);    
    im_dird = im_dir * 180 / pi;
    
    im_gx = gx1;
    im_gy = gy1;
elseif strcmp(strMethod, 'OwnCentral')
    myfilter = fspecial('gaussian',[3 3], SigmaSize);
    im_gray = imfilter(im_gray, myfilter, 'replicate');
    
    hx1 = [1 -1];% hx2 = hx1 * -1;
    hy1 = [1; -1];% hy2 = hy1 * -1;
    
    gx1 = imfilter(double(im_gray), hx1);% gx = single((gx1 + gx2)/2);
    %gx1(:,end) = []; gx1(:,1) = [];
%     figure; imagesc(gx1)
    
    gy1 = imfilter(double(im_gray), hy1); %gy = single((gy1 + gy2)/2);
    %gy1(end,:) = []; gy1(1,:) = [];
%     figure; imshowpair(gx1, gy1, 'montage'); colormap hot
    im_grad = gx1.*gx1+gy1.*gy1;    
    im_grad = sqrt( im_grad ); % gradient norm
%     figure; imagesc(im_grad)
    %im_grad = im_grad(2:end, 2:end);
%     im_grad(1:end,1) = 0;
%     im_grad(1,1:end) = 0;
    im_dir = atan2(gx1, -gy1);
    %im_dir = im_dir(2:end, 2:end);    
    im_dird = im_dir * 180 / pi;
    
    im_gx = gx1;
    im_gy = gy1;
    
%     im_gx = imfilter(im_gray, [-1 2 -1], 'replicate');
elseif strcmp(strMethod, 'SobelWeight')
    myfilter = fspecial('gaussian',[3 3], SigmaSize);
    I = imfilter(double(im_gray), myfilter, 'replicate');
    
    h = [-1 -4 -1; 0 0 0; 1 4 1];
    gx1 = imfilter(I,h','replicate');
    gy1 = imfilter(I,h,'replicate');
    
    im_grad = gx1.*gx1+gy1.*gy1;    
    im_grad = sqrt( im_grad )/6; 
    im_dir = atan2(gx1, -gy1);
    im_dird = im_dir * 180 / pi;
        
    gx1 = abs(gx1)/6; gy1 = abs(gy1)/6; 
%     [ nmsX, nmsY, memX, memY] = partial_nonmaxsup( gx1, gy1, false,1);
%     figure; imagesc([nmsX nmsY]); axis image;
%     figure; imagesc([memX memY]); axis image;
    
%     im_gx = imfilter(im_gray, [-1 2 -1], 'replicate');
else
    fprintf('Method is wrong!\n');
end

end