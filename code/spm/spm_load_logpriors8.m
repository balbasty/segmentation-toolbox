function logtpm = spm_load_logpriors8(V,deg)
% Load the tissue probability maps for segmentation
% FORMAT tpm = spm_load_priors8(V)
% V   - structures of image volume information (or filenames)
% tpm - a structure for tissue probabilities
%
% This function is intended to be used in conjunction with spm_sample_priors.
% V = spm_vol(P);
% T = spm_load_priors(V);
% B = spm_sample_priors(T,X,Y,Z);
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_load_priors8.m 5962 2014-04-17 12:47:43Z spm $

if nargin<2, deg  = 2;     end

if ~isstruct(V), V = spm_vol(V); end
spm_check_orientations(V);

logtpm.V = V;
logtpm.M = logtpm.V(1).mat;
logtpm.d = logtpm.V(1).dim;

Kb = numel(logtpm.V);
logtpm.dat = cell(Kb,1);
for k1=1:(Kb)
    logtpm.dat{k1} = zeros(logtpm.V(1).dim(1:3),'single');
end

spm_progress_bar('Init',logtpm.V(1).dim(3),'Loading priors','Planes loaded');
for i=1:logtpm.V(1).dim(3)
    M = spm_matrix([0 0 i]);
    for k1=1:Kb
        logtpm.dat{k1}(:,:,i) = spm_slice_vol(logtpm.V(k1),M,logtpm.V(1).dim(1:2),0);
    end
    spm_progress_bar('Set',i);
end
logtpm.bg1 = zeros(Kb,1);
for k1=1:Kb
    logtpm.bg1(k1) = mean(mean(logtpm.dat{k1}(:,:,1)));
    logtpm.bg2(k1) = mean(mean(logtpm.dat{k1}(:,:,end)));
end

logtpm.deg = deg;

spm_progress_bar('Clear');
