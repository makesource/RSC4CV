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

imshowpair(img, template, 'montage');

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
img = imgaussfilt(img, sigma);
template = imgaussfilt(template, sigma);
img = double(edge(img, 'Canny'));
template = double(edge(template, 'Canny'));

imshowpair(img, template, 'montage');

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
        if template(tem_r, tem_c)
            max_r = max(max_r, tem_r);
            min_r = min(min_r, tem_r);
            max_c = max(max_c, tem_c);
            min_c = min(min_c, tem_c);
        end
        endgit
end

% scale invariant
j = 1;
% for j = 0.2:0.1:2
    scale = imresize(template, sizeparam * j);
    [r, c] = size(scale);
    for i = 0:15:360
        rot = imrotate(scale, i, 'nearest', 'crop');
        corr = xcorr2(img, rot);
        out = max(out, corr(round(r/2): round(n+r/2-1), round(c/2): round(m+c/2-1)));
    end
%     imshow(out);
%     out = out ./ max(out(:));
% end
out = out ./ max(out(:));


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
for i =1:n
    for j=1:m
        if (out(i,j) > 0.7) 
            plot(j / sizeparam,i / sizeparam,'g.','MarkerSize', 10);
            viscircles([j / sizeparam, i / sizeparam], (tem_r + tem_c) / 4, 'Color', 'b');
        end
    end
end
