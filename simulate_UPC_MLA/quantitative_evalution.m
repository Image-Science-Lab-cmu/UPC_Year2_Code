% set up parameters
% close all; clear; clc;
function [] = quantitative_evalution(ids)
addpath UPC_MLA_functions;
addpath ../utils
addpath ../bm3d_matlab_package
addpath ../bm3d_matlab_package/bm3d
addpath ../cgsSolver

profile off;
profile on -history;

simulationNames = {
    'TOLED_0.0047m_-0.000100m_0.000420m_Opening_4.200000e+00_Mag_4.200000e+00_pitch_336',...
    'TOLED_0.0047m_-0.000165m_0.000693m_Opening_4.200000e+00_Mag_4.200000e+00_pitch_504', ...
    };

% Reference opening
% POLED:       0.237
% L2Inv Repeat 0.227
% L2Inv Random 0.226

operators = {'matmul', 'matmul', 'matmul', 'matmul'};
openRatios = [0.060, 0.119, 0.238, 0.238]; % ref open ratio 0.119

% % TODO: use deblur parameters of TOLED
% load('sample_output/densePSF/TOLED_0.0047m_-0.000050m_0.000420m_Opening_8.400000e+00_Mag_8.400000e+00_pitch_336/tuned_parameters.mat');

for id = ids
    simulationName = simulationNames{id};
    operator = operators{id};
    openRatio = openRatios(id);
    
    simulationType = 'densePSF';
    srcImgDir = 'test_data/';
    srcImgName = dir([srcImgDir, '*.png']);
    PSFDimension = 1;
    
    %     load(sprintf('sample_output/%s/%s/tuned_parameters.mat', ...
    %         simulationType, simulationName));
    %% generate PSF matrix
    load(sprintf('sample_output/%s/%s/PSFs.mat', ...
        simulationType, simulationName));
    load(sprintf('sample_output/%s/%s/tuned_parameters.mat', ...
        simulationType, simulationName));
    % save folder
    mkdir(sprintf('sample_output/%s/%s/', ...
        simulationType, simulationName));
    
    if strcmp(operator, 'matmul')
        yPSFMatrix = construct_yPSFMatrix(PSFs_y); % use spatially-varying in y-direction
        omega = yPSFMatrix;
    else
        omega = PSFs;
    end
    
    %% for SNR; for different images
    SNRs = 24:4:40;
    curr_ssims = zeros(length(SNRs), 1);
    curr_psnrs = zeros(length(SNRs), 1);
    
    for SNR = SNRs
        
        mean_psnr = 0;
        mean_ssim = 0;
        imgIds = [1];
        
        for imgId = imgIds
            
            best_psnrVal = 0;
            best_ssimVal = 0;
            %             for lambda = [1e-2]
            
            % load image
            img = im2double(imread([srcImgDir, '/', srcImgName(imgId).name]));
            img = img ./ max(img(:));
%                             refRatio = 0.25;
            refRatio = 0.119; % TODO: all testing have 0.119 open ratio
            img = img * openRatio / refRatio;
            
            isSaved = 0;        % predefined %TODO
            lambda = best_lambdas(SNRs == SNR);
            noise_var = best_noise_vars(SNRs == SNR);
            if isSaved == 0
                % if image size if not 1024x2048, scale and pad to 1024x2048.
                %             [img_pad, paddings, ratio] = pad_scale_image(img, [1024,2048]);
                %             [img_pad, paddings, ratio] = pad_scale_image(img, [1024,2048]);
                
                % simulate blurry image and add noise
                imgBlurnoisy = capture(img, omega, operator, SNR);
                
                % denoise (optional) and deblur
                imgSharp = deblur_cgs(imgBlurnoisy, omega, operator, noise_var, lambda);
                
                % if image size if not 1024x2048, crop and scale back to
                % original size.
                %             imgSharp = imgSharp(paddings.x0+1:end-paddings.x1, ...
                %                 paddings.y0+1:end-paddings.y1, :);
                %             imgSharp = imresize(imgSharp, 1/ratio);
                %             imgSharp = imgSharp(1:size(img,1), 1:size(img,2), :);
            else
                imgSharp = im2double(imread(sprintf('sample_output/%s/%s/%s_deblurImg_L2Gradient_lambda_%.5f_SNR_%d.png', ...
                    simulationType, simulationName, srcImgName(imgId).name, 1e-2, SNR)));
            end
            
            % Intensity compensation
            img = img * refRatio / openRatio;
            imgSharp = imgSharp * refRatio / openRatio;
            
            % compute PSNR and SSIM
            psnrVal = psnr(imgSharp, img);
            [ssimVal, ssimMap] = ssim(imgSharp, img, 'Radius', 1.5);
            
            if psnrVal > best_psnrVal
                best_psnrVal = psnrVal;
                best_ssimVal = ssimVal;
                best_imgSharp = imgSharp;
                best_ssimMap = ssimMap;
                best_lambda = lambda;
            end
            
            mean_psnr = mean_psnr + best_psnrVal;
            mean_ssim = mean_ssim + best_ssimVal;
            
            fprintf('%s, SNR=%d, lambda=%.5f, imgId=%d, PSNR=%.2f, SSIM=%.2f\n', ...
                simulationName, SNR, best_lambda, imgId, best_psnrVal, best_ssimVal);
            
            if mod(imgId-1, 5) == 0 && isSaved == 0
                imwrite(imgBlurnoisy, ...
                    sprintf('sample_output/%s/%s/%s_blurImg_SNR_%d.png', ...
                    simulationType, simulationName, srcImgName(imgId).name, SNR));
                
                imwrite(best_imgSharp, ...
                    sprintf('sample_output/%s/%s/%s_deblurImg_L2Gradient_lambda_%.5f_SNR_%d.png', ...
                    simulationType, simulationName, srcImgName(imgId).name, 1e-2, SNR));
                
                %                 fprintf('%s, SNR=%d, lambda=%.5f, imgId=%d, PSNR=%.2f, SSIM=%.2f\n', ...
                %                 simulationName, SNR, lambda, imgId, psnrVal, ssimVal);
                
                close all;
                % visualize the ssimMap
                figure('Renderer', 'painters', 'Position', [10, 10, 1200, 500]);
                grid on; grid minor;
                set(gcf,'Color',[1 1 1], 'InvertHardCopy','off');
                imagesc(best_ssimMap(:,:,2)); colorbar; caxis([0, 1])
                saveas(gcf, sprintf('sample_output/%s/%s/%s_ssimMap_SNR_%d.png', ...
                    simulationType, simulationName, srcImgName(imgId).name, SNR));
                
            elseif mod(imgId-1, 5) == 0 && isSaved ~= 0
                close all;
                % visualize the ssimMap
                figure('Renderer', 'painters', 'Position', [10, 10, 1200, 500]);
                grid on; grid minor;
                set(gcf,'Color',[1 1 1], 'InvertHardCopy','off');
                imagesc(best_ssimMap(:,:,2)); colorbar;
                saveas(gcf, sprintf('sample_output/%s/%s/%s_ssimMap_SNR_%d.png', ...
                    simulationType, simulationName, srcImgName(imgId).name, SNR));
            end
            
        end
        
        mean_ssim = mean_ssim / length(imgIds);
        mean_psnr = mean_psnr / length(imgIds);
        
        curr_ssims(SNRs==SNR) = mean_ssim;
        curr_psnrs(SNRs==SNR) = mean_psnr;
        
        save(sprintf('sample_output/%s/%s/sweep_snr.mat', ...
            simulationType, simulationName), ...
            'SNR', 'curr_ssims', 'curr_psnrs');
    end
end
