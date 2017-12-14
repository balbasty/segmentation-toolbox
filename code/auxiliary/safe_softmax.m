function [Q,ll] = safe_softmax(Q)
maxQ   = nanmax(Q,[],2);
Q      = exp(bsxfun(@minus,Q,maxQ));
sQ     = nansum(Q,2);
if nargout==2
    ll = nansum(log(sQ) + maxQ);
end
Q      = bsxfun(@rdivide,Q,sQ);
%=======================================================================