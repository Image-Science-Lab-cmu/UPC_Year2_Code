function [lin_HDR, noise_var] = mergeHDR(bursts, low, high, sensor)

% imgNames = dir(sprintf('%s/*.dng', folder));
N = length(bursts);
% load radioCalibrate.mat
cumu_valid = 0;
cumu_var = 0;   % variance for regions with some valids
mean_var = 0;  % variance for regions with no valids

lin_HDR = 0;
mean_HDR = 0;
for i = 1: N
    
    % --- process raw image ---
    lin_rgb = bursts(i).burst;
    
    % Segment valid pixel according to gray value
    lin_gray = rgb2gray(lin_rgb);
    if i == 1 % Assume the shortest exposure
        valid = (lin_gray >= low);
    elseif i == N % Assume the longest exposure
        valid = (lin_gray <= high);
    else
        valid = (lin_gray >= low) .* (lin_gray <= high);
    end
    
    % Accumulate valid values and average values
    r = 1;
    lin_HDR = lin_HDR + lin_rgb .* valid ./ bursts(i).exposureTime * r;
    mean_HDR = mean_HDR + lin_rgb ./ bursts(i).exposureTime * r;
        
%         noise_var = noise_var + lin_rgb .* valid ./ (bursts(i).exposureTime)^2 * r^2;
%         mean_noise_var = mean_noise_var + lin_rgb ./ (bursts(i).exposureTime)^2 * r^2;
    
    cumu_var = cumu_var + valid ./ bursts(i).exposureTime;
    mean_var = mean_var + 1 ./ bursts(i).exposureTime;
    
    cumu_valid = cumu_valid + valid;  % N
    
end
cumu_mask = cumu_valid > 0;

lin_HDR = cumu_mask .* lin_HDR ./ max(1, cumu_valid) + ...
          (1 - cumu_mask) .* mean_HDR ./ N;

% Affine noise model
noise_var = cumu_mask .* cumu_var ./ (max(1, cumu_valid).^2) + ...
            (1 - cumu_mask) .* mean_var ./ (N^2);
noise_var = lin_HDR .* sensor.gain .* noise_var;

% TODO: max HDR value as 1
% lin_HDR = lin_HDR ./ max(lin_HDR(:));
% noise_var = noise_var ./ (max(lin_HDR(:))^2);
end