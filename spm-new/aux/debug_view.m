function debug_view(type,fig,lkp,varargin)
if ~isempty(fig)
    buf = varargin{1};
    nz  = numel(buf);
    dz  = floor(nz/2) + 1;
    d   = [size(buf(1).msk) nz];        
    Kb  = max(lkp);
    K   = numel(lkp);  
    
    slice = NaN(d(1:2));
    
    prob_colmap = [zeros(128,1), linspace(0,1,128)', linspace(1,0,128)';
                   linspace(0,1,128)', linspace(1,0,128)', zeros(128,1)];
    clims_CT    = [0 100];
               
    set(0,'CurrentFigure',fig);
    
    % Responsibilities
    if strcmp(type,'responsibilities')
        gmm = varargin{2}; 
        mg  = varargin{3}; 
        wp  = varargin{4}; 
        
        Q = zeros([d K],'single');
        for z=1:nz
            qt               = zeros([prod(d(1:2)) K],'single');
            qt(buf(z).msk,:) = latent(buf(z).f,buf(z).bf,mg,gmm,buf(z).dat,lkp,wp,buf(z).code);            
            Q(:,:,z,:)       = reshape(single(qt),[d(1:2) 1 K]);
        end        
                 
        for k=1:K              
            subplot(3,K,k);
            slice = Q(:,:,floor(d(3)/2) + 1,k);
            imagesc(slice'); axis image xy off; title(['q, k=' num2str(lkp(k))]); colormap(pink);               
            
            subplot(3,K,K + k);
            slice = permute(Q(:,floor(d(2)/2) + 1,:,k),[3 1 2]);
            imagesc(slice); axis image xy off; title(['q, k=' num2str(lkp(k))]); colormap(pink);   
            
            subplot(3,K,2*K + k);
            slice = permute(Q(floor(d(1)/2) + 1,:,:,k),[2 3 1]);
            imagesc(slice'); axis image xy off; title(['q, k=' num2str(lkp(k))]); colormap(pink);   
        end 
    end
    
    % Bias field
    if strcmp(type,'bf')
        modality = varargin{2};
        
        N = numel(buf(dz).f);
        for n=1:N
            subplot(N,3,n); 
            slice(buf(dz).msk) = buf(dz).f{n};   
            if strcmp(modality,'CT'), imagesc(slice',clims_CT); axis image xy off; title(['X, n=' num2str(n)]); colormap(gray);    
            else                      imagesc(slice'); axis image xy off; title(['X, n=' num2str(n)]); colormap(gray);    
            end
            
            subplot(N,3,N + n);            
            slice(buf(dz).msk) = buf(dz).bf{n};
            imagesc(slice'); axis image xy off; colormap(gray); colorbar;
            title(['B, n=' num2str(n)]);
            
            subplot(N,3,2*N + n);            
            slice(buf(dz).msk) = buf(dz).bf{n}.*buf(dz).f{n};
            if strcmp(modality,'CT'), imagesc(slice',clims_CT); axis image xy off; title(['X, n=' num2str(n)]); colormap(gray);    
            else                      imagesc(slice'); axis image xy off; title(['X, n=' num2str(n)]); colormap(gray);    
            end
            title(['BX, n=' num2str(n)]);
        end  
    end
        
    % Template
    if strcmp(type,'template')    
        wp    = varargin{2};
        b     = buf(dz).dat;
        alpha = bsxfun(@times,b,wp);
        alpha = bsxfun(@rdivide,alpha,sum(alpha,2));
        
        K1 = floor(sqrt(Kb));
        K2 = ceil(Kb/K1);                                               
        for k=1:Kb            
            subplot(K1,K2,k);            
            slice(buf(dz).msk) = alpha(:,k);
            imagesc(slice'); axis image xy off; colormap(pink);
            title(['k=' num2str(k), ', wp=' num2str(round(wp(k),2))]);
        end 
    end

    % Convergence
    if strcmp(type,'convergence')    
        L  = varargin{2};
        F  = numel(L);
        F1 = floor(sqrt(F));
        F2 = ceil(F/F1);                                               
        for f=1:F            
            subplot(F1,F2,f);            
            plot(0:numel(L{f}) - 1, L{f},'-');
            if     f==1, title('ll');
            elseif f==2, title('llrb');
            elseif f==3, title('llr');
            end
        end
    end
    
    drawnow
end
%==========================================================================   