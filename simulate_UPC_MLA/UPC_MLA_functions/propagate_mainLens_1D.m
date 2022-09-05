function [I3_lensOnly] = propagate_mainLens_1D(params, thetas, lambda)
addpath ../utils/
addpath ../fourier_optics_package/
% clear;
% close all;

%% optical system parameters
batchSize = length(thetas);
% % D0 = 4.2e-3;                    % diam of the disp aperture [m]
% D1 = 3.475e-3;                    % diam of the lens / display aperture [m]
% % D1 = 3.277e-3;
k = 2*pi/lambda;                % wavenumber
% f = 5.6e-3;                     % main lens focal length (m)
% 
% % delta1 = 0.05e-6;               % diam of the lens aperture [m]
% % delta2 = 0.05e-6;               % diam of the lens aperture [m]
% % M = 2^18;                       % number of samples
% 
% delta1 = 1.5e-8;               % diam of the lens aperture [m]
% delta2 = 1.5e-8;               % diam of the lens aperture [m]
% % delta3 = 1.5e-7;               % todo: spacing on sensor plane
% delta3 = 1.45e-7;               % todo: spacing on sensor plane
% % M = 2^20;                       % number of samples
M = floor(lambda*params.f/params.delta2/params.delta3/2)*2;
% 
% % repeated structure
% pitch = 168e-6;
% 
% assert(D1 == params.lensApertureSize, 'D2 Error');
% assert(f == params.f, 'main lens f Error');
% assert(delta1 == params.delta1, 'delta1 error');
% assert(delta2 == params.delta2, 'delta2 error');
% assert(delta3 == params.delta3, 'delta3 errro');
%% Display plane: initialize aperture and display pattern

% load OLED display
u1_1 = load_display(M, params.delta1, params.pitch, params.openRatio);

% main lens aperture
x1 = (-M/2: 1: M/2-1) * params.delta1;
mainLensAperture = rect(x1/params.lensApertureSize);
u1_1 = u1_1 .* mainLensAperture;
u1_1 = repmat(u1_1, [batchSize, 1]);

% Tilted wavefront
thetas = 90-thetas;         % propagating direction (degree)
alphas = cos(thetas / 180 * pi);
tilted_phase = exp(1i*k*alphas'*x1);

u1_1 = u1_1 .* tilted_phase;
clear tilted_phase;


% lens propagation
[u3_lensOnly,~,~] = propFF_1D_batch(u1_1, M*params.delta1, lambda, params.f, 0);
I3_lensOnly = abs(u3_lensOnly).^2;
clear u3_lensOnly;

%% crop PSF and find the center location

% % D3 = 0.006;                        % sensor size (2048*sensorPitch)
% % sensorPitch = 2e-6;
% % I3: crop to sensor size
% % delta3 = lambda * f / delta2 / M;  % sensor plane spacing  % todo
% 
% x3 = (-M/2: 1: M/2-1) * params.delta3;
% sensorWindow = rect(x3/params.sensorSize);
% % scale = round(sensorPitch / delta3);
% 
% I3_lensOnly = I3_lensOnly(:, sensorWindow);
% I3_lensOnly_full = I3_lensOnly;
% I3_lensOnly = downSample1D(I3_lensOnly, params.scale);
% I3_lensOnly = gather(I3_lensOnly);
% 
% % figure, hold on; plot(I3/sum(I3(:))), plot(I3_lensOnly/sum(I3_lensOnly(:))); hold off;
% % dtheta = 90 - 180/pi*acos(scale*lambda/M/delta1);

% Anqi @ Oct 28
% Since we know sensor has 1024 pixels and has a pixel pitch of 2um
sensorRes = round(params.sensorSize / params.sensorPitch);
x3 = (-M/2: 1: M/2-1);
sensorWindow = rect(x3 / (sensorRes * params.scale));
I3_lensOnly = I3_lensOnly(:, sensorWindow);
I3_lensOnly = downSample1D(I3_lensOnly, params.scale);
I3_lensOnly = gather(I3_lensOnly);

%% I3 info
% sensor.N = M;
% sensor.M = length(I3_lensOnly);
% sensor.pitch = 2e-6;
% sensor.x = (0:1:sensor.M-1) - floor(sensor.M/2);
% % sensor.dtheta = dtheta;
% sensor.scale = params.scale;
% sensor.delta1 = params.delta1;
% % sensor.pixel2angle = @(pid) 90 - 180/pi*acos(pid*params.scale*lambda/M/params.delta1);
% t = params.delta2 / params.delta1;
% sensor.pixel2angle = @(pid) 90 - 180/pi*acos(pid*params.sensorPitch*t/params.f);
% sensor.angle2pixel = @(theta) round(cos(pi/180*(90-theta))*M*params.delta1/params.scale/lambda);
