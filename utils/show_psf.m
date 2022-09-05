function [] = show_psf(psf)

kernels = psf;

psf = psf ./ max(psf, [], 'all');
psf = log(psf);    
pixelSize = 2e-6;

figure, hold on; axis equal, axis tight;
psf_x = (1:size(psf,1))*pixelSize;
psf_y = (1:size(psf,2))*pixelSize;
imagesc(psf_x, psf_y, psf(:,:,2));colormap jet; colorbar;
hold off


% figure, imshow((kernels ./ max(kernels(:)) * 2000).^(1/2.2),[]);


