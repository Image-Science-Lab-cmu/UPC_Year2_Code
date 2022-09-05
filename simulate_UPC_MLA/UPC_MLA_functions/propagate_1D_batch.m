function [I3] = propagate_1D_batch(params, thetas, lambda, shape)
% Note: It seems that Matlab can multiply
% a gpuArray with double array. It implicitly
% transfer double array to GPU.


addpath ../utils/
addpath ../fourier_optics_package/
DEBUG_ONE = 0;
DEBUG_TWO = 1;
% clear;
% close all;

% optical system parameters
batchSize = length(thetas); 
k = 2*pi/lambda;                % wavenumber
M = floor(lambda*params.f/params.delta2/params.delta3/2)*2; % number of samples

%% Anqi @ Nov 10th, 2021
% Display plane: initialize aperture and display pattern
x1 = (-M/2: 1: M/2-1) * params.delta1;
displayAperture = rect(x1/params.displayApertureSize);
displayApertureM = sum(displayAperture);

u1 = load_display(displayApertureM, params.delta1, params.pitch, params.openRatio);
% u2_aperture = load_display(displayApertureM, params.delta1, params.pitch, 1/3);
microLensArray1 = load_microLensArray(k, displayApertureM, params.delta1, params.pitch, params.f1, 'flat_spherical', params.offset1, params.openRatio);
microLensArray2 = load_microLensArray(k, displayApertureM, params.delta2, params.pitch2, params.f2, 'spherical', params.offset2);

% pad patterns to M
padM = ceil((M-displayApertureM)/2);
u1 = padarray(u1, [0,padM], 1,'both');
% u2_aperture = padarray(u2_aperture, [0, padM], 1, 'both');u2_aperture=u2_aperture(1:M);
microLensArray1 = padarray(microLensArray1, [0,padM], 1,'both');
microLensArray2 = padarray(microLensArray2, [0,padM], 1,'both');
u1 = u1(1:M);microLensArray1=microLensArray1(1:M); microLensArray2=microLensArray2(1:M);

% prepare input field
u1 = u1 .* displayAperture;
u1 = repmat(u1, [batchSize, 1]);
% Tilted wavefront
thetas = 90-thetas;         % propagating direction (degree)
alphas = cos(thetas / 180 * pi);
tilted_phases = exp(1i*k*alphas'*x1);
u1 = u1 .* tilted_phases;

clear tilted_phases;

%% LensArray1: microlens array at display plane
u1 = u1 .* microLensArray1;
clear microLensArray1;

%% Fresnel Propagation: from display plane z (m)
[~,u2] = propASP_1D_batch(u1, lambda, params.delta1, params.delta2, params.z);
clear u1;
%% LensArray2: microlens array at lens plane
mainLensAperture = rect(x1/params.lensApertureSize);
u2 = u2 .* microLensArray2 .* mainLensAperture;    
clear microLensArray2 mainLensAperture;   
%% Lens propagation: right after positive microlens, and before the main
% lens, propagate from main lens to sensor.
[u3,~,~] = propFF_1D_batch(u2, M*params.delta2, lambda, params.f, 0);
I3=abs(u3).^2;
clear u2 u3;

%% crop PSF and downsample to sensor pitch
x3 = (-M/2: 1: M/2-1);
sensorWindow = rect(x3 / (params.sensorRes * params.scale));
I3 = I3(:, sensorWindow);
I3 = downSample1D(I3, params.scale);
I3 = gather(I3);
