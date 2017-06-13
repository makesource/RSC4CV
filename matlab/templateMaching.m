clc;
clear;
close all;
imtool close all;

img = imread('../img/img4.png');
template = imread('../img/tmp4.png');

imshow(img);
h = impoly;
position = wait(h);
boundbox = [min(position(:,1)), ....
      min(position(:,2)), ....
      max(position(:,1))-min(position(:,1)), ....
      max(position(:,2))-min(position(:,2))];
BW = createMask(h);
img = imcrop(uint8(BW).*img, boundbox);
img2 = img;
img = rgb2gray(img);
I = imadjust(img);
imshow(I);

template = rgb2gray(template);
img = double(img) / 255;
template = double(template) / 255;

% img = imgaussfilt(img, 3);
% template = imgaussfilt(template, 3);
% imshow(img);

img = double(edge(img, 'Canny', [], 3));
template = double(edge(template, 'Canny', [], 3));
imshow(img);
figure;
imshow(template);
[r2, c2] = size(template);
img = imresize(img, 0.4);
template = imresize(template, 0.4);

[n, m] = size(img);
[r, c] = size(template);
corr = xcorr2(img, template);
out = corr(r/2:n+r/2-1,c/2:m+c/2-1);

for i = 0:15:360
    %rot = imrotate(template, i);
    rot = imrotate(template,i,'nearest','crop');
    %imshow(rot);
    corr = xcorr2(img, rot);
    out = max(out, corr(r/2:n+r/2-1,c/2:m+c/2-1));
end

out = out ./ max(out(:));
imshow(out);

[out_r, out_c] = size(out);
windowsize = round(min(size(template))/4);
for i = 1:out_r
    for j = 1:out_c
        left = max(i-windowsize, 1);
        right = min(i+windowsize, out_r);
        top = max(j-windowsize, 1);
        bottom = min(j+windowsize, out_c);
        if out(i,j) ~= max(max(out(left:right,top:bottom)))
            out(i,j) = 0;
        end
    end
end
imshow(img2);
hold on;
for i =1:out_r
    for j=1:out_c
        if (out(i,j) > 0.5) 
            plot(j / 0.4,i / 0.4,'g.','MarkerSize', 5);
            viscircles([j/0.4 i/0.4],(r2+c2)/4,'Color','b');
        end
    end
end

