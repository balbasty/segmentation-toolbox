function reset_origin(P,wskew)
if nargin < 2, wskew = false; end

V   = spm_vol(P);
M   = V.mat;
dim = V.dim;

if wskew    
    % Get scale, skew and offset components
    [Scale,Skew,Offs] = decompose_mat(M,dim);

    % New orientation matrix
    M = Offs*Skew*Scale;
else
    vox  = sqrt(sum(M(1:3,1:3).^2));
    
    if det(M(1:3,1:3))<0
        vox(1) = -vox(1); 
    end;
    
    orig = (dim(1:3)+1)/2;
    off  = -vox.*orig;
    M    = [vox(1) 0      0      off(1)
               0      vox(2) 0      off(2)
               0      0      vox(3) off(3)
               0      0      0      1];
end

spm_get_space(P,M);          
%==========================================================================    

%========================================================================== 
function [Scale,Skew,Offs] = decompose_mat(M,dim)
% Inspiration: https://github.com/cran/AnalyzeFMRI/blob/master/R/niftiFMRI.R (line: 2187)

x = M(1:3,1);
y = M(1:3,2);
z = M(1:3,3);

sx = sqrt(sum(x.^2));

sy = sqrt(sum(y.^2)-(sum(x.*y))^2/(sx)^2);

a = sum(x.*y)/(sx*sy);

x0 = x/sx;
y0 = y/sy - a*x0;

sz = sqrt(sum(z.^2)-(sum(x0.*z))^2-(sum(y0.*z))^2);

b = sum(x0.*z)/sz;
c = sum(y0.*z)/sz;

% skew component
Skew = [1 a b 0;
        0 1 c 0;
        0 0 1 0;
        0 0 0 1];
    
if det(M(1:3,1:3)) < 0, 
    sx = -sx; % flip
end; 

% scale component
Scale = diag([sx sy sz 1]);
    
% offset component
orig        = (dim(1:3)+1)/2;
off         = -[sx sy sz].*orig;
Offs        = eye(4);
Offs(1:3,4) = off;    
%==========================================================================