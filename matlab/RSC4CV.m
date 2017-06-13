%% initial - clear

clc; clear all; close all; imtool close all;

%% initial - image load

img_extension = '*.bmp;*.cur;*.gif;*.hdf;*.ico;*.jpeg;*.jpg;*.pbm;*.pcx;*.pgm;*.png;*.pnm;*.ras;*.tif;*.tiff;*.xwd';

[filename, pathname, filterindex] = uigetfile({img_extension, 'All Readable Files'}, 'select voting board');

if isnumeric(filename)
    disp('file is not selected');
    return;
end

img = imread(strcat(pathname, filename));

%% initial - select region

imshow(img);
h = impoly;
position = wait(h);
boundbox = [min(position(:,1)), min(position(:,2)), max(position(:,1))-min(position(:,1)), max(position(:,2))-min(position(:,2))];
BW = createMask(h);
img = imcrop(uint8(BW).*img, boundbox);

imshow(img);
h = impoly;
position = wait(h);
boundbox = [min(position(:,1)), min(position(:,2)), max(position(:,1))-min(position(:,1)), max(position(:,2))-min(position(:,2))];
BW = createMask(h);
template = imcrop(img, boundbox);

% imshowpair(img, template, 'montage');

%% process - template background vanishing part
%
%
% 
% 
% 
% 
%% process - gaussian and egde filtering

sigma = 3;
img2 = img;
img = rgb2gray(img);
template = rgb2gray(template);
img = double(img) / 255;
template = double(template) / 255;
img = double(edge(img, 'Canny', [], sigma));
template = double(edge(template, 'Canny', [], sigma));

% imshowpair(img, template, 'montage');

%% process - template matching

% resizing
sizeparam = 0.4;
img = imresize(img, sizeparam);
[n, m] = size(img);
out = zeros(n, m);
[tem_r, tem_c] = size(template);
max_r = 1;
min_r = tem_r;
max_c = 1;
min_c = tem_c;
for i = 1:tem_r
    for j = 1:tem_c
        if template(i, j) > 0.5
            max_r = max(max_r, i);
            min_r = min(min_r, i);
            max_c = max(max_c, j);
            min_c = min(min_c, j);
        end
    end
end
template = imcrop(template, [min_c, min_r, max_c - min_c + 1, max_r - min_r + 1]);
% [tem_r, tem_c] = size(template);
% imshow(template);
% scale invariant
j = 1;
% for j = 0.2:0.1:2
    scale = imresize(template, sizeparam * j);
    [r, c] = size(scale);
%     imshow(scale ./ max(scale(:)));
    for i = 0:15:360
        rot = imrotate(scale, i, 'nearest', 'crop');
        corr = xcorr2(img, rot);
        out = max(out, corr(round(r/2): round(n+r/2-1), round(c/2): round(m+c/2-1)));
    end
%     imshow(out);
%     out = out ./ max(out(:));
% end
out = out ./ max(out(:));
% imshow(img2);
% figure;
% imshow(out);
% figure;

windowsize = round(min(size(template))/4);
for i = 1:n
    for j = 1:m
        left = max(i-windowsize, 1);
        right = min(i+windowsize, n);
        top = max(j-windowsize, 1);
        bottom = min(j+windowsize, m);
        if out(i,j) ~= max(max(out(left:right,top:bottom)))
            out(i,j) = 0;
        end
    end
end

imshow(img2);
hold on;
matched = [];
counter = 0;
for i =1:n
    for j=1:m
        if (out(i,j) > 0.7)
            counter = counter + 1;
            matched(counter, :) = [i / sizeparam, j / sizeparam];
        end
    end
end

img3 = img2;
% imshow(img3);

cform = makecform('srgb2lab');
[n2, m2, ~] = size(img2);
lab_he = applycform(img2,cform);
img2 = double(lab_he(:,:,2:3));
img2 = reshape(img2, n2* m2 ,2);
nColors = 7;
% repeat the clustering 3 times to avoid local minima
[cluster_idx, cluster_center] = kmeans(img2, nColors, 'distance', 'sqEuclidean', 'Replicates', 3);
pixel_labels = reshape(cluster_idx, n2, m2);
imshow(pixel_labels,[]);

classified = uint32(zeros(nColors, 4));
matched = round(matched);
for i = 1:counter
    class = pixel_labels(matched(i, 1), matched(i, 2));
    classified(class, 1) = classified(class, 1) + 1;
    classified(class, 2) = classified(class, 2) + uint32(img3(matched(i, 1), matched(i, 2), 1));
    classified(class, 3) = classified(class, 3) + uint32(img3(matched(i, 1), matched(i, 2), 2));
    classified(class, 4) = classified(class, 4) + uint32(img3(matched(i, 1), matched(i, 2), 3));
end

v = sum(classified(:, 1) > 0);
imshow(img3);
title(sprintf('There are %d stickers!', counter));

for i =1:n
    for j=1:m
        if (out(i,j) > 0.7)
            counter = counter + 1;
            matched(counter, :) = [i / sizeparam, j / sizeparam];
            plot(j / sizeparam, i / sizeparam,'g.','MarkerSize', 10);
            viscircles([j / sizeparam, i / sizeparam], (max_r - min_r + max_c - min_c) / 4, 'Color', 'b');
        end
    end
end

figure;
c = 1;
for i = 1:nColors
    if classified(i, 1) == 0
        continue
    end
    color = classified(i, 2:4) ./ classified(i, 1);
    block = zeros(100, 100, 'uint8');
    block(:,:,1) = color(1);
    block(:,:,2) = color(2);
    block(:,:,3) = color(3);
    subplot(1, v, c), imshow(block);
%     figure; imshow(block);
    xlabel(sprintf('%d', classified(i, 1)));
    c = c + 1;
end