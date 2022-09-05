function [lin_deblur] = deblur_cgs_hdr(lin_blur, omega, operator, lambda, noise_covariance)

addpath ../bm3d_matlab_package
addpath ../bm3d_matlab_package/bm3d
addpath ../cgsSolver

% %% denoise
% if noise_var ~= 0
%     noiseProfile = BM3DProfile();
%     noiseProfile.gamma = 0;
%     noise_type =  'gw';
%     seed = 0; % seed for pseudorandom noise realization
%     if strcmp(operator, 'two_conv')
%         [~, PSD, ~] = getExperimentNoise(noise_type, noise_var, seed, size(lin_blur(:,:,:,1)));
%         lin_blur(:,:,:,1) = CBM3D(lin_blur(:,:,:,1), PSD, noiseProfile);
%         lin_blur(:,:,:,2) = CBM3D(lin_blur(:,:,:,2), PSD, noiseProfile);
%     else
%         [~, PSD, ~] = getExperimentNoise(noise_type, noise_var, seed, size(lin_blur));
%         lin_blur = CBM3D(lin_blur, PSD, noiseProfile);
%     end
% end

    
%% deblur
% Using cgs to solve both spatially varying and invariant blur.
lin_deblur = LS_L2_Deblurring_3Ch(lin_blur, omega, [], lambda, operator, noise_covariance);