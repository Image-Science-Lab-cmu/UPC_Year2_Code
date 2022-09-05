function [yPSFMatrix] = construct_yPSFMatrix_v2(PSFs_y)

% convert spatially varying PSF to a PSF matrix
% Input:
%   PSFs_y: number of PSF x PSF dim x 3 
% Output:
%   yPSFMatrix: each row is one PSF, dimension 1024 + 2048 - 1
%               number of row is number of spatially-varying PSF, ie. 1024

%number of pixels on sensor
H = 1024;                               % ground-truth image height

if size(PSFs_y, 2) < 2048
    % TOLED PSF has length of 1024.
    % Pad it to 2048 so that every PSF has the same length.
    pad = ceil((2048-size(PSFs_y, 2))/2);
    PSFs_y = padarray(PSFs_y, [0, pad, 0], 'both');
    PSFs_y = PSFs_y(:, 1:2048, :);
end

kH = size(PSFs_y, 2);                   % vertical(y) PSF length
yPSFMatrix = zeros(H+kH-1, H, 3);       % PSF matrix

% If PSF is computed for 2048 incoming directions, 
% take the center 1024 directions.
if size(PSFs_y, 1) == 2048              
    PSFs_y = PSFs_y(513: 1536, :, :);
end

% Construct spatially-varying PSF matrix
for cc = 1: 3
    % y PSF
    for iy = 1: H
        if size(PSFs_y, 1) > 1
            PSF_y = PSFs_y(iy, :, cc);
        else
            PSF_y = PSFs_y(1, :, cc);
        end
        yPSFMatrix(iy:iy+kH-1, iy, cc) = PSF_y;
    end
end
yPSFMatrix = yPSFMatrix(floor(kH/2)+1:floor(kH/2)+H, :, :);


% % for new generated PSFs_y
% yPSFMatrix = zeros(size(PSFs_y));
% for cc = 1: 3
%     for iy = 1: 1024
%         yPSFMatrix(:, iy, cc) = PSFs_y(iy, :, cc);
%     end
% end

end