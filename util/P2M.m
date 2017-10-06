%==========================================================================
% function M = P2M(P)
%==========================================================================
function M = P2M(P)
% Polar decomposition parameterisation of affine transform,
% based on matrix logs

% Translations
D      = P(1:3);
D      = D(:);

% Rotation part
ind    = [2 3 6];
T      = zeros(3);
T(ind) = -P(4:6);
R      = expm(T-T');

% Symmetric part (zooms and shears)
ind    = [1 2 3 5 6 9];
T      = zeros(3);
T(ind) = P(7:12);
V      = expm(T+T'-diag(diag(T)));

M      = [V*R D ; 0 0 0 1];
return;