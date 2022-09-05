function [microLensArray, hmap] = load_microLensArray(k, M, delta, pitch, f, shape, offset, lensOpenRatio)

if ~exist('offset', 'var')
    offset = 0;
end

switch shape
    
    case 'spherical'
    
        patternN = round(pitch / delta);
        x1_2 = (-patternN/2: 1: patternN/2-1) * delta;
        
%         todo:
        nt =  1.515;   
        R = (nt-1) * f;
        hmap1 = - x1_2.^2 / 2 * (1/R);
        
        % todo: quantized by 1um !!
        quantize_delta = 0.2e-6;
        hmap1 = floor(hmap1 / quantize_delta) * quantize_delta;
        microLens1 = exp(1i*k*(nt-1)*hmap1);

%         microLens1 = exp(-1i*k/2/f * x1_2.^2);

        % create negative micro-lens array
        repNum = ceil(M * delta / pitch);
        microLensArray = repmat(microLens1, [1, repNum]);
        hmap = repmat(hmap1, [1, repNum]);
        
        microLensArray = microLensArray(:, 1:M);
        hmap = hmap(:, 1:M);

%         % todo:!! create a single lenslet, and zero amplitude else where
%         leftNum = floor(repNum / 2);
%         microLensArray(1:leftNum*patternN+1) = 0;
%         microLensArray((leftNum+1)*patternN+1:end) = 0;

    case 'flat_spherical'
        if ~exist('lensOpenRatio', 'var')
            lensOpenRatio = 1/4;
        end
        
        patternN = round(pitch / delta);
        sphericalN = round(pitch * lensOpenRatio / delta);
        flatN = round((patternN - sphericalN) / 2);
        x1_2 = (-sphericalN/2: 1: sphericalN/2-1) * delta;
        
        nt =  1.515;   
        R = (nt-1) * f;
        hmap1 = - x1_2.^2 / 2 * (1/R);
        
        % todo: quantized by 1um
%         quantize_delta = 0.3e-6;
%         quantize_ratio = round(quantize_delta / delta);
%         N = floor(sphericalN / quantize_ratio) * quantize_ratio;
%         hmap1 = hmap1(floor(quantize_ratio/2): quantize_ratio: end);
%         hmap1 = kron(hmap1, ones(1,quantize_ratio));
%         
%         if sphericalN>N; hmap1 = [hmap1, hmap1(end)*ones(1,sphericalN-N)];
%         else; hmap1=hmap1(1:sphericalN);end;
%         spherical_phase = exp(1i*k*(nt-1)*hmap1);
%         
%         microLens1 = [ones(1, flatN), spherical_phase, ones(1, flatN)];
        
        microLens1 = [ones(1, flatN), exp(-1i*k/2/f * x1_2.^2), ones(1, flatN)];
        hmap1 = [zeros(1, flatN), hmap1, zeros(1, flatN)];
        
        % create negative micro-lens array
        repNum = ceil(M * delta / pitch);
        microLensArray = repmat(microLens1, [1, repNum]);
        hmap = repmat(hmap1, [1, repNum]);
        
% %         todo: quantize at 1um
%         quantize_delta = 1e-6;
%         quantize_ratio = round(quantize_delta / delta);
%         N = floor(length(microLensArray) / quantize_ratio) * quantize_ratio;
%         microLensArray = resample(microLensArray(1:N), 1, quantize_ratio);
%         microLensArray = kron(microLensArray, ones(1,quantize_ratio));
%         
%         if M>N; microLensArray = [microLensArray, zeros(1,M-N)];end;
        microLensArray = microLensArray(:, 1:M);
        hmap = hmap(1, 1:M);
        
        % zero amplitude elsewhere
%         leftNum = floor(repNum / 2);
%         microLensArray(1:leftNum*patternN+1) = 0;
%         microLensArray((leftNum+1)*patternN+1:end) = 0;

    case 'phasePlate'
        
        patternN = round(pitch / delta);
        x1_2 = (-patternN/2: 1: patternN/2-1) * delta;
        n = 1.458;
        R = (n-1)*f;
        hmax = sign(R)*50e-6;
%         hmax = sign(R) * 0.53e-6 / (n-1) * 20;
        hmap1 = mod(x1_2.^2/2/R, hmax);
        microLens1 = exp(-1i*k*(n-1)*hmap1);
        
        % create negative micro-lens array
        repNum = ceil(M * delta / pitch);
        microLensArray = repmat(microLens1, [1, repNum]);
        microLensArray = microLensArray(:, 1:M);
        hmap = repmat(hmap1, [1, repNum]);
        hmap = hmap(:, 1:M);
        
        
    case 'constant'
        microLensArray = ones(1, M);
        hmap = ones(1, M);
        
    otherwise
        fprintf('Error: unknown MLA shape profile.\n');
        
end


if offset ~= 0
    % misalign the microlens array by offset
    assert(offset < pitch, 'Error: Microlens offset is larger than the pitch.');
    offsetN = round(offset / delta);
    microLensArray = [ones(1, offsetN), microLensArray];
    microLensArray = microLensArray(:, 1:M);
end