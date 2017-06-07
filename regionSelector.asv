img = imread('circuit.tif');
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