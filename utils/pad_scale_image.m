function [img, paddings, ratio] = pad_scale_image(img, sz)

if ~exist('sz', 'var')
    sz = [1024, 2048];
end

% By default, we don't scale image.
ratio = 1;

if size(img, 1) > sz(1)
    ratio = sz(1) / size(img, 1);
    img = imresize(img, ratio);
    img = img(1: sz(1), :, :);
end

if size(img, 2) > sz(2)
    ratio = sz(2) / size(img, 2);
    img = imresize(img, ratio);
    img = img(:, 1: sz(2), :);
end

[w, h, ~] = size(img);
paddings.x0 = floor((sz(1)-w)/2);     % number of pixel padding to the top
paddings.x1 = sz(1)-w-paddings.x0;    % number of pixel padding to the bottom
paddings.y0 = floor((sz(2)-h)/2);     % number of pixel padding to the left
paddings.y1 = sz(2)-h-paddings.y0;    % number of pixel padding to the right

img = padarray(img, [paddings.x0, paddings.y0], 'pre');
img = padarray(img, [paddings.x1, paddings.y1], 'post');
end