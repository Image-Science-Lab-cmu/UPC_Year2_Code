function [loss] = wiener_loss(PSF, ratio)

% Compute PSF-induced loss defined in Yang & Sankaranayan ICCP'21
% PSF:  h x w x c 
% loss: scalar
if ~exist('ratio', 'var')
    ratio = 0.3;
end
    
[h, w, c] = size(PSF);
system_filters = zeros(1, h*w*c);

% compute system frequency response and save all values
for cc = 1: c
    K = PSF(:,:,cc) / sum(sum(PSF(:,:,cc)));
%     K = PSF(:,:,cc);
    fft_K = ifftshift(fft2(fftshift(K)));
    
    % wiener filter
    wiener_filter = conj(fft_K) ./ (abs(fft_K).^2 + 0.015);
    system_filter = abs(fft_K .* wiener_filter);
    
    system_filters((cc-1)*h*w+1: (cc-1)*h*w+h*w) = system_filter(:);
end

% compute loss function
len = length(system_filters);
sorted_system_filters = sort(system_filters(:), 'ascend');
score = mean(sorted_system_filters(1:ceil(ratio*len)));

loss = -score;

% Note: Wiener score needs to be maximized, loss nedds to be minimized.

end