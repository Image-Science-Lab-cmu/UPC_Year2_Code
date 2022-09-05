function dn_array = downSample(array, ratio)

N = size(array, 1);
N = N - mod(N, ratio);
array = array(1: N, 1: N);

array = imfilter(array, ones(ratio, ratio)/ ratio^2);
mid = floor(ratio / 2) + 1;
dn_array = array(mid:ratio:end, mid:ratio:end, :);
% Note: python conv2d samples the center point of a 
% convolution filter, while we take the first sample
% in Matlab previously.