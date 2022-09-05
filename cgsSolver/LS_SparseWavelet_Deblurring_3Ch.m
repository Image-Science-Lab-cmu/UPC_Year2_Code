
function recons = LS_SparseWavelet_Deblurring_3Ch(y, omega, init, beta)

[H, W, ~] = size(y); %shape of image
[kH, kW, ~] = size(omega);
recons = zeros(kH+H-1, kW+W-1, 3);

for cc = 1:3
    recons(:,:,cc) = leastSquareWaveletSparseDeblurring(y(:,:,cc), omega(:,:,cc), init, beta);
end

end


function xhat1 = leastSquareWaveletSparseDeblurring(y, k1, init, beta)

dwtmode('per');



siz= size(y); siz = siz(1:2);
siz_y = siz;
siz_x = siz+size(k1)-1;

vec = @(x) x(:);
fA = @(x) [vec(conv2(reshape(x, siz_x), k1, 'valid'))];
fAadj = @(x) ( conv2(reshape(x, siz_y), k1(end:-1:1, end:-1:1), 'full'));

xinit = fAadj(vec(y));
[~, book_keeping] = wavedec2(xinit, 8, 'db4');
psi = @(x) wavedec2(x, 8, 'db4');
psi_adj = @(y) waverec2(y, book_keeping, 'db4');


% beta = 0.5;
max_iter = 50;
xhat1 = IST_backtracking(vec(y), beta, psi, psi_adj, fA, fAadj, max_iter);
xhat1 = reshape(xhat1, siz_x);
end