function [psf_out, lenN] = crop_psf_center(psf, lenN)

if ~exist('lenN', 'var')
    lenN = 64;
end

[~, maxId] = max(psf);

error = 0;
while error < 0.95
    lenN = lenN * 2;
    psf_out = psf(maxId - round(lenN/2): maxId + round(lenN/2) - 1);
    error = sum(psf_out) / sum(psf);
end
end