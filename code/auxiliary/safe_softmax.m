function [Q,ll] = safe_softmax(Q)
maxQ   = max(Q,[],2);
Q      = exp(bsxfun(@minus,Q,maxQ));
sQ     = sum(Q,2);
if nargout==2
    ll = sum(log(sQ) + maxQ);
end
Q      = bsxfun(@rdivide,Q,sQ);
%=======================================================================