%% Load captured HDR images and corresponding PSFs

addpath ../utils/;

[PSFConfigs, SceneConfigs] = mainCaptureConfig();

for psfId = 1: length(PSFConfigs)
    
    % load PSF
    load([PSFConfigs(psfId).srcDir, PSFConfigs(psfId).srcName]);
    % normlize RGB channels
    for cc = 1: 3; PSF(:,:,cc) = PSFs(:,:,cc) / sum(PSFs(:,:,cc), [1,2]); end
    % Rot PSF if specified
    if PSFConfigs(psfId).rot; PSF = permute(PSF, [2,1,3]); end
    
    for sceneId = 1: length(SceneConfigs)
        
        fprintf('%s, %s\n', PSFConfigs(psfId).srcName, SceneConfigs(sceneId).name);
        %  Captured HDR images
        folder = PSFConfigs(psfId).dstDir;
        scene = SceneConfigs(sceneId).name;

        HDR = load([folder, scene, '/HDR_Captured.mat']); 
        yVar = HDR.hdrBlurIm_noiseVar;  HDR = HDR.hdrBlurIm;
        % Scale yVar to similar scale as the regularization term (heuristic)
        sensor.capacity = 15506;
        sensor.gain = 1/sensor.capacity;
        yVar = yVar ./ (sensor.gain^2);
        yVar = max(1e-9, yVar);  % clip small values for stable division
        clear hdrBlurIm;

        %% Deblur HDR for one display pattern
        tonemap = @(x, a, b) ((x ./ (x + a)) .^ (1 / b));
        %  Reference deblur using Wiener deconvolution
        [sharpWiener, ~] = myWienerDeconv(HDR, PSF, 28);
        save([folder, scene, '/HDR_deblurred_Wiener.mat'], 'sharpWiener');

        %  Deblur using cgs with SV noise
        %  tune parameters for iterative solver
        for max_threshold = [0.1, 1, 10, 100, 1000]
            for lambda = 1/max_threshold * [1, 0.1, 0.01, 0.001]
                sharpSV = deblur_cgs_hdr(HDR, PSF, 'conv', lambda, min(max_threshold, yVar));
                save([folder, scene, '/tune2_HDR_deblurred_SVNoise_cgs_maxClip_', num2str(max_threshold), '_lambda_', num2str(lambda), '.mat'], 'sharpSV');
            end
        end

        % Deblur using cgs with spatially-invariant noise
        for lambda = [1, 0.1, 0.01, 0.001]
            sharp = deblur_cgs_hdr(HDR, PSF, 'conv', lambda, []);
            save([folder, scene, '/HDR_deblurred_cgs_lambda_' num2str(lambda), '.mat'], 'sharp');
        end
%         end

        %% imshow 
%         figure(1), hold on
%         tonemap = @(x, a, gamma) (x ./ (x+a)) .^ (1/gamma);
%         subplot(131), imshow(tonemap(sharpWiener, 0.01, 3)); title('HDR Wiener');
%         subplot(132), imshow(tonemap(sharp, 0.01, 3)); title('HDR CGS');
%         subplot(133), imshow(tonemap(sharpSV, 0.01, 3)); title('HDR CGS SV Noise');
%         hold off;
    end
end