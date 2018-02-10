function [Q,ll] = safe_softmax(Q)
maxQ = max(Q,[],2);
Q    = exp(bsxfun(@minus,Q,maxQ));
sQ   = nansum(Q,2);
ll   = nansum(log(sQ)+maxQ);
Q    = bsxfun(@rdivide,Q,sQ);
%=======================================================================