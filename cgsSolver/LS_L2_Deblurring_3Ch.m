function recons = LS_L2_Deblurring_3Ch(y, omega, init, lambda, operator, noise_covariance)

[H, W, ~] = size(y);        %shape of image
recons = zeros(H, W, 3);
% todo
if strcmp(operator, 'matmul')
    [~, nH, ~] = size(omega);
    recons = zeros(nH, W, 3);
end
% ===

for cc = 1:3
    if isempty(noise_covariance)
        recons(:,:,cc) = leastSquareDeblurring(y(:,:,cc,:), omega(:,:,cc,:), ...
            init, lambda, operator, []);
    else
        recons(:,:,cc) = leastSquareDeblurring(y(:,:,cc,:), omega(:,:,cc,:), ...
            init, lambda, operator, noise_covariance(:,:,cc,:));
    end
end

end



function recons = leastSquareDeblurring(y, omega, init, lambda, operator, noise_covariance)
cgs_iters = 150;  % Number of conjugate gradient descent iterations
cgs_tol = 1e-6;  % cgs tolerance
% lambda = 1e-3;   % Regularization parameter

% todo: read in blurry image
% img = imread('cameraman.tif');
% img = double(img)/255;
% 
[H, W] = size(y); %shape of image
[kH, kW] = size(omega);
% y = AOmega(img, omega);

% If Nvidia GPU available, send y and omega
% to GPU and all relavent computation will
% be conducted on GPU.
if gpuDeviceCount > 0 
    omega = gpuArray(omega);
    y = gpuArray(y);
end

vec = @(x) x(:);

if strcmp(operator, 'matmul')
    
    % y = Omega * x
    % y:    [H, W]
    % Omega [H, kH]
    % x:    [kH, W]
    
    [H, W] = size(y); %shape of image
    [~, kH] = size(omega);
    
    AOmega = @(x) omega*x;
    AOmega_adj = @(y) omega'*y;

    Aforward = @(x) vec(AOmega(reshape(x,kH,W)));
    Aadj = @(y) vec(AOmega_adj(reshape(y,H,W)));
    
elseif strcmp(operator, 'conv')
    
    AOmega = @(x) myConv2(x, omega, 'valid');
    AOmega_adj = @(x) myConv2(x, omega(end:-1:1, end:-1:1), 'full');
    
    Aforward = @(x) vec(AOmega(reshape(x,H+kH-1,W+kW-1)));
    Aadj = @(x) vec(AOmega_adj(reshape(x,H,W)));
    

elseif strcmp(operator, 'two_conv')
    [H, W, ~, ~] = size(y);
    [kH, kW, ~, ~] = size(omega);
    omega1 = omega(:,:,:,1);
    omega2 = omega(:,:,:,2);
    y1 = y(:,:,:,1);
    y2 = y(:,:,:,2);
    
    AOmega1 = @(x) myConv2(x, omega1, 'valid');
    AOmega1_adj = @(x) myConv2(x, omega1(end:-1:1, end:-1:1), 'full');
    
    Aforward1 = @(x) vec(AOmega1(reshape(x,H+kH-1,W+kW-1)));
    Aadj1 = @(x) vec(AOmega1_adj(reshape(x,H,W)));
    
    AOmega2 = @(x) myConv2(x, omega2, 'valid');
    AOmega2_adj = @(x) myConv2(x, omega2(end:-1:1, end:-1:1), 'full');
    
    Aforward2 = @(x) vec(AOmega2(reshape(x,H+kH-1,W+kW-1)));
    Aadj2 = @(x) vec(AOmega2_adj(reshape(x,H,W)));
    
else
    error('Error: undefined operator');
end

if strcmp(operator, 'two_conv')
    if isempty(noise_covariance)
        A = @(x) Aadj1(Aforward1(x)) + Aadj2(Aforward2(x)) + lambda * x;
        b = Aadj1(y1(:)) + Aadj2(y2(:));
    else
        fprintf('Denoise with spatially-varying noise!\n');
        cov1 = vec(noise_covariance(:,:,:,1)) + 1e-12; % Add a small number to avoid numerical issues
        cov2 = vec(noise_covariance(:,:,:,2)) + 1e-12;
        A = @(x) Aadj1((1./cov1).*Aforward1(x)) + ...
                 Aadj2((1./cov2).*Aforward2(x)) + lambda * x;
        b = Aadj1((1./cov1).*y1(:)) + Aadj2((1./cov2).*y2(:));
    end
else
    if isempty(noise_covariance)
        A = @(x) Aadj(Aforward(x)) + lambda * x;
        b = Aadj(y(:));
    else
        if sum(noise_covariance < 1e-12) > 1
            error('Error: noise covariance smaller than 1e-12!');
        end
        fprintf('Denoise with spatially-varying noise!\n');
        cov = vec(noise_covariance);
        A = @(x) Aadj((1./cov).*Aforward(x)) + lambda * x;
        b = Aadj((1./cov).*y(:));
    end
end


if ~isempty(init)
    recons = cgs(A, b, cgs_tol, cgs_iters, [], [], vec(init));
else
    [recons,~,~,~,resvec] = cgs(A, b, cgs_tol, cgs_iters);
%     recons = cgs(A, b, cgs_tol, cgs_iters);
end

switch operator
    case 'matmul'
%         recons = reshape(recons, H, kH);
%         recons = recons(:, (kH-kW+1)/2: (kH-kW+1)/2+W-1);
        recons = reshape(recons, kH, W);
%         todo: for PSF matrix (H, H)
%         recons = recons((kH-H+1)/2+1:(kH-H+1)/2+H, :);
    case 'conv'
        recons = reshape(recons, H+kH-1, W+kW-1);
        recons = recons(floor(kH/2)+1: floor(kH/2)+H, floor(kW/2)+1:floor(kW/2)+W);
    case 'two_conv'
        recons = reshape(recons, H+kH-1, W+kW-1);
        recons = recons(floor(kH/2)+1: floor(kH/2)+H, floor(kW/2)+1:floor(kW/2)+W);
    otherwise
        error('Unknown operator');
end

if gpuDeviceCount > 0
    recons = gather(recons);
end
end

% figure;
% subplot(2,2,1);
% imshow(img);
% title('Original image');
% subplot(2,2,2);
% imshow(y);
% title('AOmega');
% subplot(2,2,3);
% imshow(recons);
% title('Reconstructed Image')
% subplot(2,2,4);
% imshow(abs(img-recons));
% title('Abs difference')
