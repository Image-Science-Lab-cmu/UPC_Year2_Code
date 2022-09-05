function recons = LS_Deblurring_3Ch(y, omega, init)

[H, W, ~] = size(y); %shape of image
[kH, kW, ~] = size(omega);
recons = zeros(kH+H-1, kW+W-1, 3);

for cc = 1:3
    recons(:,:,cc) = leastSquareDeblurring(y(:,:,cc), omega(:,:,cc), init);
end

end



function recons = leastSquareDeblurring(y, omega, init)
cgs_iters = 50;  % Number of conjugate gradient descent iterations
cgs_tol = 1e-6;  % cgs tolerance

% todo: read in blurry image
% img = imread('cameraman.tif');
% img = double(img)/255;
% 
[H, W] = size(y); %shape of image
[kH, kW] = size(omega);
% y = AOmega(img, omega);

vec = @(x) x(:);
AOmega = @(x, k) myConv2(x, k, 'valid');
AOmega_adj = @(x, k) myConv2(x, k(end:-1:1, end:-1:1), 'full');

Aforward = @(x) vec(AOmega(reshape(x,H+kH-1,W+kH-1),omega));
Aadj = @(x) vec(AOmega_adj(reshape(x,H,W), omega));

Dxforward = @(x) vec(Dx(reshape(x,H+kH-1,W+kH-1)));
DxAdj = @(x) vec(Dxadj(reshape(x,H+kH-1,W+kH-1-1)));

Dyforward = @(x) vec(Dy(reshape(x,H+kH-1,W+kH-1)));
DyAdj = @(x) vec(Dyadj(reshape(x,H+kH-1-1,W+kH-1)));

A = @(x) Aadj(Aforward(x));
b = Aadj(y(:));

if ~isempty(init)
    recons = cgs(A, b, cgs_tol, cgs_iters, [], [], vec(init));
else
    recons = cgs(A, b, cgs_tol, cgs_iters);
end
recons = reshape(recons, H+kH-1, W+kW-1);
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
