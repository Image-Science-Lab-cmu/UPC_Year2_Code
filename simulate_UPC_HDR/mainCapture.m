clear;clc;close all;

addpath ../utils;

%% scenes and blur kernels
%  TODO: hdr scene and display pattern can be changed in mainCaptureConfig
[PSFConfigs, SceneConfigs] = mainCaptureConfig();

for psfId = 1: length(PSFConfigs)

    fprintf('Load PSF from %s\n', PSFConfigs(psfId).srcName);
    % load PSF
    load([PSFConfigs(psfId).srcDir, PSFConfigs(psfId).srcName]);
    % normlize RGB channels
    for cc = 1: 3; PSFs(:,:,cc) = PSFs(:,:,cc) / sum(PSFs(:,:,cc), [1,2]); end
    % Rot PSF if specified
    if PSFConfigs(psfId).rot; PSFs = permute(PSFs, [2,1,3]); end
    
    results_dir = PSFConfigs(psfId).dstDir;
    
    for sceneId = 1: length(SceneConfigs)
        fprintf('   Simulate captured HDR for scene %s\n', ...
                SceneConfigs(sceneId).name);
        close all;
        dataset = SceneConfigs(sceneId).srcDir;
        scene = SceneConfigs(sceneId).name;
        suffix = SceneConfigs(sceneId).suffix;
        mkdir([results_dir, scene]);

        %% Simulate exposure stack captured by UDC

        %  1. Compute latent scene (high dynamic range)
        hdrGtIm = hdrread([dataset, scene, suffix]);
        [H, W, ~] = size(hdrGtIm);
        newH = ceil(1024 / min(H, W) * H);
        newW = ceil(1024 / min(H, W) * W);
        hdrGtIm = imresize(hdrGtIm, [newH, newW], 'box');

        %  2. Capture scene under UPC (high dynamic range)
        SNR = Inf;  % no noise at this point
        hdrGtBlurIm = capture(hdrGtIm, PSFs, 'conv', SNR);

        %  3. Capture scene with different exposures
        [H, W, ~] =  size(hdrGtBlurIm);
        N = 15;
        bursts = zeros(H, W, 3, N);
        deblurredBursts = zeros(H, W, 3, N);

        sensor.capacity = 15506;        % sensor max well capacity
        sensor.noise_std = 4.87;        % sensor readout noise std
        SNRs = 24:4:40;
        Ls = [273, 654, 1608, 4005, 10024];
        SNR = 28;
        L = Ls(SNRs == SNR);             
        sensor.gain = 1/sensor.capacity;% scale intensity between [0, 1]

        exposureTime = 1/max(hdrGtBlurIm(:))/1.2;
        for imId = 1: N
            % add photon and readout noise
            burst = add_noise(hdrGtBlurIm * L, exposureTime, sensor);
            
%             % test noise model
%             gt_burst = min(hdrGtBlurIm * L * exposureTime, sensor.capacity) * sensor.gain;
%             photon_var = (gt_burst * sensor.gain);
%             readout_var = (sensor.noise_std * sensor.gain)^2;
%             normalized_burst = (burst - gt_burst) ./ sqrt(photon_var + readout_var);
%             figure(1), imagesc(normalized_burst(:,:,2));
            

            % quantize
            burst = uint16(burst * 65535);
            burst = im2double(burst);  % convert to double[0,1] for post-processing
            bursts(:,:,:,imId) = im2double(burst);
            capturedBurstsInfo(imId).burst = burst;
            capturedBurstsInfo(imId).exposureTime = exposureTime;
            
            exposureTime = exposureTime * 2;
        end

        %% Recover sharp, high dynamic range image
        %  1. Composite HDR image
        low = 0.05;
        high = 0.95;
        [hdrBlurIm, hdrBlurIm_noiseVar] = mergeHDR(capturedBurstsInfo, low, high, sensor);
        save([results_dir, scene, '/HDR_Captured.mat'], 'hdrBlurIm', 'hdrBlurIm_noiseVar');
    
%     %% Save HDRBlurIm with no noise for future reference
%     %  Save HDRBlurIm
%         exposureTime = 1/max(hdrGtBlurIm(:))/1.2;
%         for imId = 1: N % quantize
%             burst = hdrGtBlurIm * L * exposureTime * sensor.gain;
%             burst = uint16(burst * 65535);
%             burst = im2double(burst);  % convert to double[0,1] for post-processing
%             bursts(:,:,:,imId) = im2double(burst);
%             capturedBurstsInfo(imId).burst = burst;
%             capturedBurstsInfo(imId).exposureTime = exposureTime;
%             exposureTime = exposureTime * 2;
%         end
%         low = 0.05;
%         high = 0.95;
%         hdrGtBlurIm = mergeHDR(capturedBurstsInfo, low, high, sensor);
%         save([results_dir, scene, '/HDR_Captured_NoNoise.mat'], 'hdrGtBlurIm');
    end
end