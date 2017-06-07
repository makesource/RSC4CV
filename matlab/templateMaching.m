clc;
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
imshow(img)


img = rgb2gray(img);
template = rgb2gray(template);
img = double(img) / 255;
template = double(template) / 255;

img = double(edge(img, 'Sobel'));
template = double(edge(template, 'Sobel'));
imshow(img);
figure;
imshow(template);
% 
% sigma     = 2;
% threshold = 0.1;
% [Im, Io, Ix, Iy] = myEdgeFilter(img, sigma);
% img = double(Im > threshold);
% imshow(Im > threshold);
% [Im, Io, Ix, Iy] = myEdgeFilter(template, sigma);
% template = double(Im > threshold);
% figure;
% imshow(Im > threshold);


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
imshow(img);
hold on;
for i =1:out_r
    for j=1:out_c
        if (out(i,j) > 0.7) 
            plot(j,i,'g.','MarkerSize', 5);
            %viscircles([j i],(r+c)/4,'Color','b');
        end
    end
end



