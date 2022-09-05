function recons = LS_L2Gradient_Deblurring_3Ch(y, omega, init, lambda, operator)

[H, W, ~] = size(y);        %shape of image
recons = zeros(H, W, 3);
% todo
if strcmp(operator, 'matmul')
    [~, nH, ~] = size(omega);
    recons = zeros(nH, W, 3);
end
% ===


% switch operator
%     case 'matmul'
%         kW = 2048;                     % PSF length %todo
%         [H, W, ~] = size(y);           % blur image size
%         recons = zeros(H, W, 3);
%     case 'conv'
%         recons = zeros(H, W, 3);
%     otherwise
%         assert('Unknown operator');
% end

for cc = 1:3
    recons(:,:,cc) = leastSquareDeblurring(y(:,:,cc,:), omega(:,:,cc,:), init, lambda, operator);
end

end



function recons = leastSquareDeblurring(y, omega, init, lambda, operator)
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
% if strcmp(operator, 'matmul')
% %     AOmega = @(x, k) x*k;
% %     AOmega_adj = @(x, k) x*k';
% % 
% %     Aforward = @(x) vec(AOmega(reshape(x,H,kH),omega));
% %     Aadj = @(x) vec(AOmega_adj(reshape(x,H,kW), omega));
% 
%     AOmega = @(x) x*omega;
%     AOmega_adj = @(x) x*omega';
% 
%     Aforward = @(x) vec(AOmega(reshape(x,H,kH)));
%     Aadj = @(x) vec(AOmega_adj(reshape(x,H,kW)));
%     
% 
%     Dxforward = @(x) vec(Dx(reshape(x,H,kH)));
%     DxAdj = @(x) vec(Dxadj(reshape(x,H,kH-1)));
% 
%     Dyforward = @(x) vec(Dy(reshape(x,H,kH)));
%     DyAdj = @(x) vec(Dyadj(reshape(x,H-1,kH)));

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
    
    Dxforward = @(x) vec(Dx(reshape(x,kH,W)));
    DxAdj = @(x) vec(Dxadj(reshape(x,kH,W-1)));

    Dyforward = @(x) vec(Dy(reshape(x,kH,W)));
    DyAdj = @(x) vec(Dyadj(reshape(x,kH-1,W)));

elseif strcmp(operator, 'conv')
    
%     AOmega = @(x, k) myConv2(x, k, 'valid');
%     AOmega_adj = @(x, k) myConv2(x, k(end:-1:1, end:-1:1), 'full');
    
    AOmega = @(x) myConv2(x, omega, 'valid');
%     AOmega = @(x) conv2(x, omega, 'valid'); %TODO:debug
    AOmega_adj = @(x) myConv2(x, omega(end:-1:1, end:-1:1), 'full');
    
    Aforward = @(x) vec(AOmega(reshape(x,H+kH-1,W+kW-1)));
    Aadj = @(x) vec(AOmega_adj(reshape(x,H,W)));
    
    Dxforward = @(x) vec(Dx(reshape(x,H+kH-1,W+kW-1)));
    DxAdj = @(x) vec(Dxadj(reshape(x,H+kH-1,W+kW-1-1)));
    
    Dyforward = @(x) vec(Dy(reshape(x,H+kH-1,W+kW-1)));
    DyAdj = @(x) vec(Dyadj(reshape(x,H+kH-1-1,W+kW-1)));

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
    
    Dxforward = @(x) vec(Dx(reshape(x,H+kH-1,W+kW-1)));
    DxAdj = @(x) vec(Dxadj(reshape(x,H+kH-1,W+kW-1-1)));
    
    Dyforward = @(x) vec(Dy(reshape(x,H+kH-1,W+kW-1)));
    DyAdj = @(x) vec(Dyadj(reshape(x,H+kH-1-1,W+kW-1)));
else
    error('Error: undefined operator');
end

if ~strcmp(operator, 'two_conv')
    A = @(x) Aadj(Aforward(x)) + lambda * ((DxAdj(Dxforward(x))) + (DyAdj(Dyforward(x))));
    b = Aadj(y(:));
else
    A = @(x) Aadj1(Aforward1(x)) + Aadj2(Aforward2(x)) + ...
        lambda * ((DxAdj(Dxforward(x))) + (DyAdj(Dyforward(x))));
    b = Aadj1(y1(:)) + Aadj2(y2(:));
end


if ~isempty(init)
    recons = cgs(A, b, cgs_tol, cgs_iters, [], [], vec(init));
else
    recons = cgs(A, b, cgs_tol, cgs_iters);
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
