% This script computes and saves 1D PSF from -20 degrees to 20 degrees
% clear;clc;close all;
% profile off;
% profile on -history;
function [] = compute_SV_PSFs(ids)

% Simulator for UPC with double-sided Microlens arrays
% Input:
%       ids:  System parameters are saved in MLA_Params.mat
%             Specify which set of parameters to evaluate.
%
% Output:
%       The script saves PSFs of the specified setups.
%
% Written by Anqi Yang @ August 2021

addpath UPC_MLA_functions;

simulationType = 'densePSF';
SAVE = 1;
save_rootdir = 'sample_output';
mkdir(sprintf('%s/%s/', save_rootdir, simulationType));

% dense wavelength
load('spectral_response.mat');

% load geometry_ASP/MLA_Params.mat;
load MLA_Params.mat;

for id = ids

    offset1 = 0;    % deprecated parameter
    offset2 = 0;    % deprecated parameter
    f = 4.67e-3;
    
    opening = 1/MLA_Params(id).lensOpenRatio;       % open ratio
    magnification = MLA_Params(id).magnification;   % beam magnification
    f0 = MLA_Params(id).f1;                         % focal length of left MLA
    pitch = MLA_Params(id).lensPitch;               % pitch of microlens
    
    f1 = abs(f0) * magnification;                   % focal length of right MLA
    pitch2 = pitch;                                 % deprecated parameter
    z = f0 + f1;                                    % z in ideal thin lens setup
    
    save_dir = sprintf('%s/%s/TOLED_%.4fm_%.6fm_%.6fm_Opening_%d_Mag_%d_pitch_%d/', ...
        save_rootdir, simulationType, f, f0, f1, opening, magnification, pitch*1e6);
    mkdir(save_dir);
    
%     SR.spectral = SR.spectral;
%     SR.weights = SR.weights;
    tic;

    % Image system for an under-display camera
    params = init_system(f, f0, f1, 1/opening, pitch, pitch2, offset1, offset2, z, 2e-6);
    sensor = get_sensor(params);
    
    % compute spatially-varying PSF in y-direction
    % TODO: change the number of spatially-varying PSFs
    pixelIds = -512: 1: 511;  % compute PSFs centered at 1024 pixel locations
    pixelIds = [0];           % only compute PSFs in the center
    thetas_y = sensor.pixel2angle(pixelIds);
    numPixel = length(pixelIds);
    batchSize = 1;                                  % #batchsize directions each time
    PSFs_y = zeros(numPixel, sensor.M, 3);          % 1024 PSFs in x direction (short edge of the sensor)
    
    for wvl = SR.spectral

        
        weights = SR.weights(SR.spectral == wvl, :);
        
        % get sensor information under these parameters
        I3_lensOnly = propagate_mainLens_1D(params, 0, wvl);
        total_energy = sum(I3_lensOnly);
        
        % spatially-varying PSFs for wavelength 'wvl'
        I3_all = zeros(numPixel, sensor.M);
        
        t=toc;
        fprintf('opening=%d, f0=%d, wavelength=%d, used %d s...\n', ...
            opening, f0*1e6, wvl*1e9, t);
        % batch propagation
        pre = 1;
        while pre <= numPixel
            post = min(pre + batchSize - 1, numPixel);
            I3_all(pre:post, :) = propagate_1D_batch(params, thetas_y(pre:post), wvl, 'phasePlate');
            pre = post + 1;
        end
        
        % crop I3 to the window 2w centered at pixelIds
        % normalize energy by I3_lensOnly
        I3_all = I3_all / total_energy;
        
        % save weighted averaged PSFs
        PSFs_y(1:numPixel,:,:) = PSFs_y(1:numPixel,:,:) + ...
            cat(3, I3_all*weights(1), I3_all*weights(2), I3_all*weights(3));
        
       
    end
    
    % normalize by area of green channel
    PSFs_y = PSFs_y / sum(SR.weights(:,2));
    
    if SAVE
        save([save_dir, 'PSFs.mat'], 'PSFs_y');
    end
    fprintf('\n\n\n');
    

end
end