function mom = spm_SuffStats(X,Q,mom)
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
% See spm_VBGaussiansFromSuffStats.m
%_______________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$

if nargin<2 || isempty(Q),
    Q = ones(size(X,1),1);
end

K = size(Q,2);
M = size(X,2);
if M<=8,
    cast = @uint8;
    typ  = 'uint8';
elseif M<=16,
    cast = @uint16;
    typ  = 'uint16';
elseif M<=32,
    cast = @uint32;
    typ  = 'uint32';
elseif M<=64,
    cast = @uint64;
    typ  = 'uint64';
else,
    error('Too many dimensions.');
end

if nargin<3,
    % Create empty data structure
    mom = struct('ind',[],'s0',0,'s1',[],'S2',[]);
    for i=1:2^M,
        mom(i).ind = dec2bin(i-1,M)=='1'; % Indices
        Mi         = sum(mom(i).ind);
        mom(i).s0  = zeros(1,K);     % Zeroeth moments
        mom(i).s1  = zeros(Mi,K);    % First moments
        mom(i).S2  = zeros(Mi,Mi,K); % Second moments
    end
else
    % Check compatibility of data structure
    if (numel(mom)~=2^M) || (size(mom(1).s0,2)~=K),
        error('Incorrect moment dimensions');
    end
end
Q(isnan(Q))=0;
code = zeros([size(X,1),1],typ);
for i=1:M,
    code = bitor(code,bitshift(feval(cast,isfinite(X(:,i)) & (X(:,i)~=0)),(i-1)));
end
for i=2:numel(mom),
    msk0      = mom(i).ind;
    ind       = find(code==msk0*(2.^(0:(M-1))'));
    if ~isempty(ind),
        x         = X(ind,msk0);
        for k=1:K,
            q                = Q(ind,k);
            mom(i).s0(1,k)   = mom(i).s0(1,k)   + sum(q);
            mom(i).s1(:,k)   = mom(i).s1(:,k)   + x'*q;
            mom(i).S2(:,:,k) = mom(i).S2(:,:,k) + bsxfun(@times,q,x)'*x;
        end
    end
end