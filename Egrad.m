function E = Egrad(Im, Sigma, ROI, mask)

if nargin < 4
    mask = ones(size(Im));
end

% Energy related to image gradient
Ix = ImageDerivatives2D(Im, Sigma, 'x');
Iy = ImageDerivatives2D(Im, Sigma, 'y');
GradMag = sqrt(Ix.^2 + Iy.^2);

[Y,X]=meshgrid(ROI(2):ROI(4),ROI(1):ROI(3));
ix = sub2ind(size(GradMag),Y,X);

[M, ix] = max(GradMag(ix(:)));


Mx = ImageDerivatives2D(im2double(mask), Sigma, 'x');
My = ImageDerivatives2D(im2double(mask), Sigma, 'y');
GradMagMask = sqrt(Mx.^2 + My.^2);
GradMagMask = GradMagMask./max(GradMagMask(:));

E = 1 - GradMag.*mask.*(1-GradMagMask)/M;
