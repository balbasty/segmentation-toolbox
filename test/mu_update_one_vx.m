clear;

K = 7;
N = 20;

rng('default');
rng(1);

a0 = exp(randn(K,1)); % log of tissue priors at some voxel (log-normal distribution)
a0 = a0/sum(a0); 
a0 = log(a0);

rng('default');
rng(2);

b = exp(randn(K,N)); % likelihoods 

rng('default');
rng(3);

c = exp(randn(K,N));  
c = bsxfun(@rdivide,c,sum(c)); % log of mixing weights
c = log(c);

%%
subplot(2,1,1);
ll = [];
a  = a0;
for i=1:100
    ea = exp(a);
    
    eab = bsxfun(@times,ea,b);
    r   = bsxfun(@rdivide,eab,sum(eab,1));
    r   = sum(r,2);
    
    a = log(1/N*sum(ea)*r + eps);

    ll = [ll; sum(log(sum(bsxfun(@times,ea,b),1)))]; plot(ll); drawnow;
end
disp(ll(end))

%%
subplot(2,1,2);
ll = [];
a  = a0;
for i=1:100
    ac  = bsxfun(@plus,a,c);    
    eac = exp(ac);
    
    eacb = bsxfun(@times,eac,b);
    
    r = bsxfun(@rdivide,eacb,sum(eacb,1));
    r = sum(r,2);
    
    p = 1./sum(eac,1);
    p = bsxfun(@times,p,exp(c)); % [1xN] x [KxN]
    p = sum(p,2);
    
    a = log(r./p + eps);

    ll = [ll; sum(log(sum(bsxfun(@times,eac,b),1)))]; plot(ll); drawnow;
end
disp(ll(end))