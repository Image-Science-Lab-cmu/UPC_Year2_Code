function [PSFConfigs, SceneConfigs] = mainCaptureConfig()

% Blur kernels information
PSFConfigs(1).srcDir = 'data/psf/';
PSFConfigs(1).srcName = 'TOLED_Rect_f_4.67mm_pitch_168um.mat';
PSFConfigs(1).dstDir = 'results_TOLED168/';
PSFConfigs(1).rot = false;
% 
% PSFConfigs(2).srcDir = 'data/psf/';
% PSFConfigs(2).srcName = 'POLED_f_4.67mm_pitch_336um.mat';
% PSFConfigs(2).dstDir = 'results_POLED168/';
% PSFConfigs(2).rot = false;
% 
% PSFConfigs(3).srcDir = 'data/psf/';
% PSFConfigs(3).srcName = 'L2Inv_Repeat_f_4.67mm_pitch_168um.mat';
% PSFConfigs(3).dstDir = 'results_L2InvRepeat168/';
% PSFConfigs(3).rot = false;
% 
% PSFConfigs(4).srcDir = 'data/psf/';
% PSFConfigs(4).srcName = 'L2Inv_Random_f_4.67mm_pitch_168um.mat';
% PSFConfigs(4).dstDir = 'results_L2InvRandom168/';
% PSFConfigs(4).rot = false;

% % Latent HDR Scene information
% SceneConfigs(1).srcDir = 'data/hdr/';
% SceneConfigs(1).name = 'vinesunset';
% SceneConfigs(1).suffix = '.pic';

SceneConfigs(1).srcDir = 'data/hdr/';
SceneConfigs(1).name = 'rosette';
SceneConfigs(1).suffix = '.pic';

% SceneConfigs(3).srcDir = 'data/hdr/';
% SceneConfigs(3).name = 'nave';
% SceneConfigs(3).suffix = '.pic';
% 
% SceneConfigs(4).srcDir = 'data/hdr/';
% SceneConfigs(4).name = 'memorial';
% SceneConfigs(4).suffix = '.pic';
% 
% SceneConfigs(5).srcDir = 'data/hdr/';
% SceneConfigs(5).name = 'groveD';
% SceneConfigs(5).suffix = '.pic';
% % 
% SceneConfigs(6).srcDir = 'data/hdr/';
% SceneConfigs(6).name = 'groveC';
% SceneConfigs(6).suffix = '.pic';
%       
end