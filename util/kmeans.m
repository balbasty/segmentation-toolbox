function [w,mn,vr] = kmeans(xv,nC,init)   
if nargin<3,     
    init{1} = 'rand';
    init{2} = [];
end

xv(xv==0) = NaN;

label = NaN(size(xv,1),nC);

N  = size(xv,1);
D  = size(xv,2);
M  = 2^D;   

in = bi2de(double(isfinite(xv)))+1;
%initialize means------------------------------------------------------

if strcmp(init{1},'eqspace')
    if isempty(init{2})
    st = [1/nC:1/nC:1]'; 
    else
        o  = init{2};
        st = [o:(1-o)/(nC-1):1]'; 
    end
    mu = permute(0.8*st*max(xv),[3,2,1]);

elseif strcmp(init{1},'rand')
     if ~isempty(init{2})
         xv(sum(xv,2)<=init{2}) = NaN;
         in = bi2de(double(isfinite(xv))) + 1;
     end
    mu  = permute(xv(randsample(find(in==M),nC),:),[3,2,1]);

elseif strcmp(init{1},'custom')
    if isempty(init{2})
        error('Have to specify means initialization')
    else
        st = init{2};
        mu = permute(st'*max(xv),[3,2,1]);
    end
end

for m=2:M

    conf_log = logical(de2bi(m-1,D));

    xvsub  = xv(repmat(m,[N,1])==in,conf_log);
    musub  = mu(1,conf_log,:);        
    N1     = size(xvsub,1); 
    D1     = size(xvsub,2);        
    dis    = zeros(N1,nC);

    %start iterations--------------------------------------------------               
    for it=1:20

        %labeling------------------------------------------------------
        for k=1:nC
        dis(:,k) = squeeze(sqrt(sum(bsxfun(@minus,xvsub,musub(:,:,k)).^2,2))); 
        end
        labelsub = bsxfun(@minus,dis,min(dis,[],2))==0;

        %updating means------------------------------------------------
        labelw = permute(repmat(labelsub,[1,1,D1]),[1,3,2]);
        m0     = sum(labelw,1);
        m1     = sum(labelw.*repmat(xvsub,[1,1,nC]),1);
        clear labelw
        musub  = m1./m0;
    end
    clear dis

    label(repmat(m,[N,1])==in,:) = labelsub;
    clear labelsub
end
clear xvsub

mom       = SuffStats(xv,label);
[w,mn,vr] = GaussiansFromSuffStats(mom);
%========================================================================== 

%==========================================================================
function mom = SuffStats(X,Q,mom)
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
%==========================================================================

%==========================================================================
function [w,mn,vr] = GaussiansFromSuffStats(mom)
K  = size(mom(1).s0,2);
N  = numel(mom(1).ind);
mg = zeros(1,K);
mn = zeros(N,K);
vr = zeros(N,N,K);
ll = -Inf;
for k=1:K,
    C  = eye(N);
    mu = zeros(N,1);
    s0 = 0;
    for i=1:numel(mom), s0 = s0 + mom(i).s0(k); end
    for iter=1:1024,
        s1     = zeros(N,1);
        S2     = zeros(N);
        P      = inv(C);
        s1i    = zeros(N,1);
        S2i    = zeros(N,N);
        old_ll = ll;
        ll     = 0;
        for i=1:numel(mom),
            if mom(i).s0(k),
                ind            = mom(i).ind;
                mux            = mu( ind,:);
                muy            = mu(~ind,:);
                Pyy            = P(~ind,~ind);
                Pyx            = P(~ind, ind);
                R              = Pyy\Pyx;
                Ex             = mom(i).s1(:,k)/mom(i).s0(k);
                Exx            = mom(i).S2(:,:,k)/mom(i).s0(k);
                tmp            = R*(mux-Ex)*muy';
                s1i( ind)      = mom(i).s1(:,k);
                S2i( ind, ind) = mom(i).S2(:,:,k);
                s1i(~ind)      = mom(i).s0(k)*(R*(mux-Ex) + muy);
                S2i(~ind,~ind) = mom(i).s0(k)*(R*(mux*mux'-mux*Ex'-Ex*mux'+Exx)*R' + tmp + tmp' + muy*muy' + inv(Pyy));
                S2i( ind,~ind) = mom(i).s0(k)*(R*(mux*Ex'-Exx) + muy*Ex');
                S2i(~ind, ind) = S2i( ind,~ind)';
                s1             = s1 + s1i;
                S2             = S2 + S2i;

                % Compute objective function
                S              = mom(i).s0(k)*(mux*mux') + mom(i).S2(:,:,k) - mux*mom(i).s1(:,k).' - mom(i).s1(:,k)*mux.';
                ll             = ll - mom(i).s0(k)*sum(log(diag(chol(C(ind,ind)*2*pi)))) - 0.5*trace(S/C(ind,ind));
            end
        end
        mu = s1/s0;
        C  = (S2 - s1*s1'/s0)/(s0+eps);
        %fprintf('%d\t%g\n', iter,ll);
        if ll-old_ll < 1e-12, break; end
    end
    mg(1,k)   = s0;
    mn(:,k)   = mu;
    vr(:,:,k) = C;
end
w = mg/sum(mg);
% fprintf('ML (simple)\t%g\n', ll);
% disp(mu)
% disp(C)
%==========================================================================