% function [lin_deblur] = deblur_cgs(lin_blur, omega, operator, SNR, lambda, tune_param)
function [lin_deblur] = deblur_cgs(lin_blur, omega, operator, noise_var, lambda)
% todo: tune noise variance

% close all; clear; clc;
% addpath utils
addpath ../bm3d_matlab_package
addpath ../bm3d_matlab_package/bm3d
addpath ../cgsSolver


%% denoise
% todo
% SNRs = 24:4:40;
if noise_var ~= 0
    noiseProfile = BM3DProfile();
    noiseProfile.gamma = 0;
    noise_type =  'gw';
    % noise_vars = [0.005, 0.002, 0.002, 0.001, 0.0002];
    % noise_vars = [0.005/4; 0.002/4; 0.002/16;0.001/16;0.0001/8];% SNR40
    % noise_var = noise_vars(SNRs==SNR); % Noise variance

    seed = 0; % seed for pseudorandom noise realization
    [~, PSD, ~] = getExperimentNoise(noise_type, noise_var, seed, size(lin_blur));
    lin_blur = CBM3D(lin_blur, PSD, noiseProfile);
end
    
%% deblur
% Using cgs to solve both spatially varying and invariant blur.
lin_deblur = LS_L2Gradient_Deblurring_3Ch(lin_blur, omega, [], lambda, operator);
