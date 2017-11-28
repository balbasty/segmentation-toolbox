function mom = spm_SuffStats(X,Q,code,mom)
% Generate sufficient statistics (handling missing data)
%
% FORMAT mom = spm_SuffStats(X,Q,mom)
%     X   - NxM array of data, where NaN or zeros are to be
%           ignored (N=number of samples, M=dimensionality).
%     Q   - NxK array of belonging probabilities (responsibilities).
%           Pass an empty array to compute simple sufficient stats.
%     mom - data strucure containing zeroeth, first and second moments.
%           Structure is organised such that 2^d elements are computed
%           where each element accounts for different missing data.
%_______________________________________________________________________
%
% See spm_GaussiansFromSuffStats.m
%_______________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

K = size(Q,2);
N = size(X,2);

% Check compatibility of data structure
if (numel(mom)~=2^N) || (size(mom(1).s0,2)~=K),
    error('Incorrect moment dimensions');
end

Q(isnan(Q)) = 0;

for i=2:numel(mom)
    msk0  = mom(i).ind;
    ind   = find(code==msk0*(2.^(0:(N-1))'));
    if ~isempty(ind)
        x = X(ind,msk0);
        for k=1:K,
            q                = Q(ind,k);
            mom(i).s0(1,k)   = mom(i).s0(1,k)   + sum(q);
            mom(i).s1(:,k)   = mom(i).s1(:,k)   + x'*q;
            mom(i).S2(:,:,k) = mom(i).S2(:,:,k) + bsxfun(@times,q,x)'*x;
        end
    end
end