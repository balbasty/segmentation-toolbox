function Q = softmax(Q,safe)
if nargin<2, safe = false; end
if safe
    maxQ = max(Q,[],2);
    Q    = bsxfun(@minus,Q,maxQ);
end
Q  = exp(Q);
sQ = sum(Q,2);
Q  = bsxfun(@rdivide,Q,sQ);
%==========================================================================