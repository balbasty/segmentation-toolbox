function tpm = spm_load_logpriors8(V,tiny,deg)
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

if nargin<2, tiny = 1e-4; end
if nargin<3, deg  = 2; end

if ~isstruct(V), V  = spm_vol(V); end
spm_check_orientations(V);

tpm.V = V;
tpm.M = tpm.V(1).mat;
tpm.d = tpm.V(1).dim;

Kb = numel(tpm.V);
tpm.dat = cell(Kb,1);
for k1=1:(Kb)
    tpm.dat{k1} = zeros(tpm.V(1).dim(1:3),'single');
end

spm_progress_bar('Init',tpm.V(1).dim(3),'Loading priors','Planes loaded');
for i=1:tpm.V(1).dim(3)
    M = spm_matrix([0 0 i]);
    for k1=1:Kb
        tpm.dat{k1}(:,:,i) = spm_slice_vol(tpm.V(k1),M,tpm.V(1).dim(1:2),0);
    end
    spm_progress_bar('Set',i);
end
tpm.bg1 = zeros(Kb,1);
for k1=1:Kb
    tpm.bg1(k1) = mean(mean(tpm.dat{k1}(:,:,1)));
    tpm.bg2(k1) = mean(mean(tpm.dat{k1}(:,:,end)));
end

tpm.tiny = tiny;
tpm.deg  = deg;

spm_progress_bar('Clear');
