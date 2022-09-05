function [imgBlurnoisy] = capture(img, PSF, operator, SNR)

% blur the image
% July 12: We skip blur in y-direction for now
% since it multiplies quantum efficiency to each
% channel again.

switch operator
    case 'conv'
        [kH, kW, ~] = size(PSF);        % PSF size
        [H, W, ~] = size(img);          % image size
        
        % since PSF is 2D, we pad in both direction.
        img = padarray(img, [ceil(kH/2), ceil(kW/2)], 0, 'both');
        img = img(1:kH+H-1, 1: kW+W-1, :);
        
        imgBlur = zeros(H, W, 3);
        for cc = 1: 3
            imgBlur(:,:,cc) = conv2(img(:,:,cc), PSF(:,:,cc), 'valid');
        end
        
    case 'matmul'
        [nH, ~, ~] = size(PSF);        % PSF length  %todo
        [H, W, ~] = size(img);         % image size
        imgBlur = zeros(nH, W, 3);
        for cc = 1: 3
            imgBlur(:,:,cc) = PSF(:,:,cc) * img(:,:,cc);
        end
        
    otherwise
        msg('Unknown operator');
end

%% add noise
if SNR == Inf
    imgBlurnoisy = imgBlur;
    return;
end
SNRs = 24:4:40;
Ls = [273, 654, 1608, 4005, 10024];
noise_vars = [0.005, 0.002, 0.002, 0.001, 0.0002];

sensor.capacity = 15506;
sensor.noise_std = 4.87;
L = Ls(SNRs == SNR);
sensor.gain = 1/L;

imgBlurnoisy = add_noise(imgBlur * L, 1, sensor);
