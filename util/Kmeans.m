function labels = Kmeans(xv,nC,init)   
if nargin<3,     
    init{1} = 'rand';
    init{2} = [];
end

labels = NaN(size(xv,1),nC);

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
    rng('default');
    rng(1);
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

    labels(repmat(m,[N,1])==in,:) = labelsub;
    clear labelsub
end
%========================================================================== 