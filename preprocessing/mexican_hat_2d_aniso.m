
function H = mexican_hat_2d_aniso(sigmaX, sigmaY, halfsizeX, halfsizeY)
% MEXICAN_HAT_2D_ANISO 2D anisotropic Mexican hat filter.
%   H = mexican_hat_2d_aniso(sigmaX, sigmaY) 
%       makes kernel with half-sizes = ceil(3*sigma).
%   H = mexican_hat_2d_aniso(sigmaX, sigmaY, halfsizeX, halfsizeY) 
%       lets you set kernel extents explicitly.
%
%   sigmaX, sigmaY : std devs along X and Y
%   halfsizeX, halfsizeY : half window sizes along X and Y
%
%   Returns a (2*halfsizeY+1) × (2*halfsizeX+1) kernel.

if nargin < 3 || isempty(halfsizeX)
    halfsizeX = ceil(3*sigmaX);
end
if nargin < 4 || isempty(halfsizeY)
    halfsizeY = ceil(3*sigmaY);
end

[x,y] = meshgrid(-halfsizeX:halfsizeX, -halfsizeY:halfsizeY);

% anisotropic squared distance
r2 = (x.^2)/(sigmaX^2) + (y.^2)/(sigmaY^2);

% Mexican hat / Laplacian of anisotropic Gaussian
H = (r2 - 2) .* exp(-0.5 * r2);

% normalize: zero mean and reasonable scale
H = H - mean(H(:));          % enforce zero mean
H = -(H / sum(abs(H(:))));      % scale to unit L1 norm
end
