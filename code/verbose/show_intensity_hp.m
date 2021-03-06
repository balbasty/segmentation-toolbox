function show_intensity_hp(fig,obj)
set(0,'CurrentFigure',fig); clf(fig);    
 
M = numel(obj);
 
all_ct = true;
for m=1:M
    if strcmp(obj{m}{1}.modality,'MRI')
        all_ct = false;
    end
end
 
if all_ct
    M = 1;    
end
 
N0 = 0;
for m=1:M
    pr = obj{m}{1}.segment.gmm.pr;   
    N1 = size(pr.W,1);
    
    N0 = N0 + max(N1 - 1,1);
end
 
i = 0;
for m=1:M
    pr  = obj{m}{1}.segment.gmm.pr;
    lkp = obj{m}{1}.segment.lkp;    
    K   = numel(pr.b);
    N1  = size(pr.W, 1);
        
    for n=1:N1      
        
        if N1>1
            % ND case
            if n == 1                
                continue;
            end
            ind  = [1 n];
        else
            % 1D case
            ind  = 1;
        end
        
        i      = i + 1;
        k1     = 1;
        CM     = jet(max(lkp.part));
        labels = cell(1,K);
        for k=1:K
            if k1~=lkp.part(k)
                k1 = k1 + 1;  
            end
            
            subplot(1,N0,i)
            h = plot_gaussian(pr.m(ind,k), pr.n(k)*pr.W(ind,ind,k), CM(k1,:), K > 1);
            h.LineWidth = 1;
            labels{k} = ['k' num2str(k1)];
        end
        legend(labels);        
    end
    
    if all_ct
        break;
    end
end                                  
drawnow
%==========================================================================
 
%==========================================================================
function h = plot_gaussian(mu, Lambda, colour, holdax)
 
    if nargin < 4
        holdax = false;
        if nargin < 3
            colour = 'b';
        end
    end
    if holdax
        hold on
    end
    
    if numel(mu) > 1
    
        % 2D plot
        
        mu = mu(1:2);
        Sigma2 = spm_matcomp('Inv', Lambda(1:2,1:2));
        Sigma = sqrt(Sigma2);
        [x1,x2] = meshgrid(linspace(mu(1)-3*Sigma(1,1),mu(1)+3*Sigma(1,1),25)', ...
                           linspace(mu(2)-3*Sigma(2,2),mu(2)+3*Sigma(2,2),25)');
        y = mvnpdf([x1(:) x2(:)],mu',Sigma2);
        y = reshape(y, [length(x2) length(x1)]);
        [~,h] = contour(x1,x2,y,1, 'color',colour);                
    else
        
        % 1D plot
        
        Sigma2 = 1/Lambda;
        Sigma  = sqrt(Sigma2);
        x = linspace(mu-3*Sigma,mu+3*Sigma,25)';
        y = mvnpdf(x,mu,Sigma2);
        h = plot(x, y, 'color', colour);
        set(gca,'YTickLabel',[]);
    end
    
    if holdax
        hold off
    end
%==========================================================================