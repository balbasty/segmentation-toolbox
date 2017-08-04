clc;

rng('default');
rng(1);

K=7;
N=20;

lik=exp(randn(K,N));

% Rescaling of tissue priors
w=exp(randn(K,N)*3);
w=bsxfun(@rdivide,w,sum(w));

% Naive way to do things
mu=ones(K,1)/K;

subplot(2,1,1);
ll=[];
for i=1:100
    wmu  = bsxfun(@times,mu,w);
    mnom = bsxfun(@rdivide,wmu,sum(wmu,1));
    
    ml  = bsxfun(@times,mnom,lik);
    tmp = bsxfun(@rdivide,ml,sum(ml,1));
    
    mu = sum(tmp,2);
    mu = mu/sum(mu);
    
    ll = [ll; sum(log(sum(bsxfun(@times,mnom,lik),1)),2)]; %plot(ll); drawnow;
end
disp(ll(end))

% Better way (I hope)
mu=ones(K,1)/K;

subplot(2,1,2);
ll=[];
for i=1:100
    wmu  = bsxfun(@times,mu,w);
    mnom = bsxfun(@rdivide,wmu,sum(wmu,1));
    
    ml    = bsxfun(@times,mnom,lik);
    munum = bsxfun(@rdivide,ml,sum(ml,1));
    
    p     = 1./sum(wmu,1);
    muden = bsxfun(@times,p,w);
    
    mu = sum(munum,2)./sum(muden,2);
    mu = mu/sum(mu);
    
    ll = [ll; sum(log(sum(bsxfun(@times,mnom,lik),1)),2)]; %plot(ll); drawnow;
end
disp(ll(end))