function [I_new] = batchCropWindow(I_all, pixelCoord, pixelIds, w)

numPixel = length(pixelIds);
I_new = zeros(numPixel, 2*w);

for idx = 1: numPixel
            
        I = I_all(idx, :);
        pixelId = pixelIds(idx); 
        
        % crop I3 to the window 2w centered at pixelId
        I_new(idx,:) = cropWindow(I, pixelCoord, pixelId, w);
        
end

end


function [I] = cropWindow(I, pixelCoord, pixelId, w)

% ids = (pixelCoord<=pixelId+w) .* (pixelCoord>=pixelId-w+1);
% I = I(ids>0);
% 
% % pad zeros if touch the boundary
% if min(pixelCoord) > pixelId-w+1
%     I = padarray(I, [0, min(pixelCoord)-(pixelId-w+1)], 'pre');
% end
% if max(pixelCoord) < pixelId+w
%     I = padarray(I, [0, (pixelId+w)-max(pixelCoord)], 'post');
% end

I = padarray(I, [0, w], 0, 'both');
N = length(I);
coord = -round(N/2): round(N/2)-1;

center = find(coord == pixelId);
I = I(center-w: center+w-1);

end