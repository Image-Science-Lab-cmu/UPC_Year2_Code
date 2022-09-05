function [params] = init_system(f, f1, f2, openRatio, pitch, pitch2, offset1, offset2, z, sensorPitch)


% main lens aperture [m]
if ~exist('f', 'var')
%     params.f = 5.6e-3;                         % main lens focal length (m) HUAWEI P30
    params.f = 4.67e-3;                          % main lens focal length (m) Nexus 6p
    params.f_scale = 1;
else
    params.f = f;
    params.f_scale = f/4.67e-3;
end
% params.lensApertureSize = params.f/2/sqrt(2);    % FNumber 2 for Nexus 6p
params.lensApertureSize = params.f/2;            % FNumber 2 for Nexus 6p
% Use the circular diameter (@ May 24, 2021)

if ~exist('offset2', 'var')
    params.offset2 = 0;
else
    params.offset2 = offset2;
end

if ~exist('offset1', 'var')
    params.offset1 = 0;
else
    params.offset1 = offset1;
end

if ~exist('pitch2', 'var')
    params.pitch2 = pitch;
else
    params.pitch2 = pitch2;
end

if ~exist('sensorPitch', 'var')
    params.sensorPitch = 2e-6;              % sensor pitch of a conventional camera 
else
    params.sensorPitch = sensorPitch;       % sensor pitch of customized camera
end
    
% k = 2*pi/lambda;                         % wavenumber
params.n = 1.515;                          % refractive index of the lens

params.delta1 = 1.5e-8;                    % spacing on display plane [m]
params.delta2 = 1.5e-8;                    % spacing on lens plane [m]
params.delta3 = 1.45e-7;                   % spacing on sensor plane [m]
params.delta3 = params.delta3 * params.f_scale;  % keep M at 2^20
% params.M = floor(lambda*f/delta2/delta3/2)*2; % number of samples 

% Consider the short edge as x-axis
% todo: try a smaller sensor pitch
% params.sensorPitch = 2e-6;                 % diam of sensor pitch [m]

params.sensorSize = 1024*2e-6;              % diam of sensor [m] (1024*sensorPitch)
params.sensorRes = round(params.sensorSize / params.sensorPitch);                     % number of pixels on sensor
% params.sensorSize = 2048*2e-6;             % diam of sensor [m] (2048*sensorPitch)
params.scale = round(params.sensorPitch / params.delta3);
                                           % downsample scale

% params.pitch = 168e-6;                   % diam of display pixel [m]
% params.pitch = 84e-6;                    % diam of display pixel [m]
params.pitch = pitch;                      % diam of display pixel [m]
params.openRatio = openRatio;              % open ratio of display pixel [m]


% f1 = -0.0002;                             % first microlens focal length [m]
% f2 = 0.0004;                              % second microlens focal length [m]
params.f1 = f1;
params.f2 = f2;

if ~exist('z', 'var')
    params.z=f1+f2;
else
    params.z = z;
end
% params.z=f1+f2;                             % Fresnel prop distance [m]
params.displayApertureSize = ...
    params.lensApertureSize + 2*params.z*tan(50/180*pi);
end