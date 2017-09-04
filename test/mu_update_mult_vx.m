clear;

K = 7;
N = 20;
I = 1;

rng('default');
rng(1);

a0 = exp(randn(K,I)); % log of tissue priors at some voxel (log-normal distribution)
a0 = bsxfun(@rdivide,a0,sum(a0)); 
a0 = log(a0);

rng('default');
rng(2);

b = exp(randn(K,I,N)); % likelihoods 

rng('default');
rng(3);

c0 = exp(randn(K,N));  
c0 = bsxfun(@rdivide,c0,sum(c0)); 
c0 = log(c0);
c  = [];
for i=1:I
   c(:,i,:) = c0;
end
% c  = bsxfun(@rdivide,c,sum(c)); % log of mixing weights

%%
subplot(2,1,1);
ll = [];
a  = a0;
for iter=1:100
    
    ea = exp(a);
    
    tmp = 0;
    for n=1:N        
        eab = bsxfun(@times,ea,b(:,:,n));
        r   = bsxfun(@rdivide,eab,sum(eab,1));
        tmp = tmp + r;    
    end
    
    a = log(1/N*bsxfun(@times,sum(ea),tmp) + eps);

    ll = [ll; sum(sum(log(sum(bsxfun(@times,ea,b),1)),2),3)]; plot(ll); drawnow;
end
disp(ll(end))

%%
subplot(2,1,2);
ll = [];
a  = a0;
for iter=1:100

    ac   = bsxfun(@plus,a,c);    
    eac0 = exp(ac);
        
    anum = 0;
    aden = 0;
    for n=1:N
   
        ac  = bsxfun(@plus,a,c(:,:,n));    
        eac = exp(ac);
    
        eacb = bsxfun(@times,eac,b(:,:,n));

        r    = bsxfun(@rdivide,eacb,sum(eacb,1));
        anum = anum + r;
        
        pn = 1./sum(eac,1);
        cn = exp(c(:,:,n));
        
        p = bsxfun(@times,pn,cn);
        aden = aden + p;
    end
    
    a = log(anum./aden + eps);

    ll = [ll; sum(sum(log(sum(bsxfun(@times,eac0,b),1)),2),3)]; plot(ll); drawnow;
end
disp(ll(end))