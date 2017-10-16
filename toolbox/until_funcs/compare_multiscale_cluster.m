function [ output_args ] = compare_multiscale_cluster( I0, param )
%GET_MULTISCALE_IMAGE Summary of this function goes here
%   Detailed explanation goes here
%I0 = im_gray;
I1 = conv2(double(I0), fspecial('gaussian', 3, .4)); I1(1,:) = []; I1(end,:) = []; I1(:,1) = []; I1(:,end) = [];
I2 = conv2(double(I0), fspecial('gaussian', 7, .8)); I2(1:3,:) = []; I2(end-2:end,:) = []; I2(:,1:3) = []; I2(:,end-2:end) = [];
I3 = conv2(double(I0), fspecial('gaussian', 11, 1.6)); I3(1:5,:) = []; I3(end-4:end,:) = []; I3(:,1:5) = []; I3(:,end-4:end) = [];

% ls0 = lsd(double(I0)); ls1 = lsd(double(I1)); ls2 = lsd(double(I2)); ls3 = lsd(double(I3)); 

% figure;
% subplot(2,2,1); imagesc(I0); axis image; colormap gray;
% subplot(2,2,2); imagesc(I1); axis image; colormap gray;
% subplot(2,2,3); imagesc(I2); axis image; colormap gray;
% subplot(2,2,4); imagesc(I3); axis image; colormap gray;
% 
% figure;
% subplot(2,2,1); imagesc(I0); hold on; for i = 1:size(ls0,1), plot(ls0(i, 1:2:3), ls0(i, 2:2:4), 'r-'); end ;axis image; colormap gray;
% subplot(2,2,2); imagesc(I1); hold on; for i = 1:size(ls1,1), plot(ls1(i, 1:2:3), ls1(i, 2:2:4), 'r-'); end ;axis image; colormap gray;
% subplot(2,2,3); imagesc(I2); hold on; for i = 1:size(ls2,1), plot(ls2(i, 1:2:3), ls2(i, 2:2:4), 'r-'); end ;axis image; colormap gray;
% subplot(2,2,4); imagesc(I3); hold on; for i = 1:size(ls3,1), plot(ls3(i, 1:2:3), ls3(i, 2:2:4), 'r-'); end ;axis image; colormap gray;
% 
% figure;
% subplot(2,2,1); imagesc(imgradient(I0));axis image
% subplot(2,2,2); imagesc(imgradient(I1));axis image
% subplot(2,2,3); imagesc(imgradient(I2));axis image
% subplot(2,2,4); imagesc(imgradient(I3));axis image

[img0, imd0] = imgradient(I0); imdr0 = imd0 / 180 * pi; imd0(imd0 < 0) = imd0(imd0 < 0) + 180;
[img1, imd1] = imgradient(I1); imdr1 = imd1 / 180 * pi; imd1(imd1 < 0) = imd1(imd1 < 0) + 180;
[img2, imd2] = imgradient(I2); imdr2 = imd2 / 180 * pi; imd2(imd2 < 0) = imd2(imd2 < 0) + 180;
[img3, imd3] = imgradient(I3); imdr3 = imd3 / 180 * pi; imd3(imd3 < 0) = imd3(imd3 < 0) + 180;

im_nms0 = nonmaxsup(img0, imd0, 1.2, false);
im_nms1 = nonmaxsup(img1, imd1, 1.2, false);
im_nms2 = nonmaxsup(img2, imd2, 1.2, false);
im_nms3 = nonmaxsup(img3, imd3, 1.2, false);

[pixel_usage, lsc_id_map0, lsc_set, line_segment, pixel_cluster] = Get_Line_Candidate(I0, img0, im_nms0, imdr0, imd0, param, true, false);
figure; imshow(I0); hold on; colormap gray;
for k = 1:size(line_segment, 1)
    x1 = line_segment(k,1:2) + line_segment(k,3)/2*[cos(line_segment(k,4)) sin(line_segment(k,4))];
    x2 = line_segment(k,1:2) - line_segment(k,3)/2*[cos(line_segment(k,4)) sin(line_segment(k,4))];
    plot([x1(1) x2(1)], [x1(2) x2(2)], 'r-'); % updated
end
    
[pixel_usage, lsc_id_map1, lsc_set, line_segment, pixel_cluster] = Get_Line_Candidate(I1, img1, im_nms1, imdr1, imd1, param, true, false);
figure; imshow(I1); hold on; colormap gray;
for k = 1:size(line_segment, 1)
    x1 = line_segment(k,1:2) + line_segment(k,3)/2*[cos(line_segment(k,4)) sin(line_segment(k,4))];
    x2 = line_segment(k,1:2) - line_segment(k,3)/2*[cos(line_segment(k,4)) sin(line_segment(k,4))];
    plot([x1(1) x2(1)], [x1(2) x2(2)], 'r-'); % updated
end

[pixel_usage, lsc_id_map2, lsc_set, line_segment, pixel_cluster] = Get_Line_Candidate(I2, img2, im_nms2, imdr2, imd2, param, true, false);
figure; imshow(I2); hold on; colormap gray;
for k = 1:size(line_segment, 1)
    x1 = line_segment(k,1:2) + line_segment(k,3)/2*[cos(line_segment(k,4)) sin(line_segment(k,4))];
    x2 = line_segment(k,1:2) - line_segment(k,3)/2*[cos(line_segment(k,4)) sin(line_segment(k,4))];
    plot([x1(1) x2(1)], [x1(2) x2(2)], 'r-'); % updated
end

[pixel_usage, lsc_id_map3, lsc_set, line_segment, pixel_cluster] = Get_Line_Candidate(I3, img3, im_nms3, imdr3, imd3, param, true, false);
figure; imshow(I3); hold on; colormap gray;
for k = 1:size(line_segment, 1)
    x1 = line_segment(k,1:2) + line_segment(k,3)/2*[cos(line_segment(k,4)) sin(line_segment(k,4))];
    x2 = line_segment(k,1:2) - line_segment(k,3)/2*[cos(line_segment(k,4)) sin(line_segment(k,4))];
    plot([x1(1) x2(1)], [x1(2) x2(2)], 'r-'); % updated
end

figure;
subplot(2,2,1); imagesc(lsc_id_map0);axis image
subplot(2,2,2); imagesc(lsc_id_map1);axis image
subplot(2,2,3); imagesc(lsc_id_map2);axis image
subplot(2,2,4); imagesc(lsc_id_map3);axis image


output_args.I1 = I1;
output_args.I2 = I2;
output_args.I3 = I3;

end

