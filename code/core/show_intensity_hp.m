function show_intensity_hp(fig,obj)
set(0,'CurrentFigure',fig); clf(fig);    

M = numel(obj);
for m=1:M
    pr = obj{m}{1}.segment.gmm.pr;
    
    K = numel(pr.b);
    N = size(pr.W, 1);
    
    i = 0;
    for n=1:N

        if N > 1
            % ND case
            if n == 1
                continue;
            end
            ind = [1 n];
            ncol = N-1;
        else
            % 1D case
            ind = 1;
            ncol = 1;
        end

        i = i+1;
        for k=1:K
            subplot(1,ncol,i)
            h = plot_gaussian(pr.m(ind,k), pr.n(k)*pr.W(ind,ind,k), 'red', K > 1);
            h.LineWidth = 2;
        end
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
        [~,h] = contour(x1,x2,y,1, 'color', colour);
        
    else
        
        % 1D plot
        
        Sigma2 = 1/Lambda;
        Sigma  = sqrt(Sigma2);
        x = linspace(mu-3*Sigma,mu+3*Sigma,25)';
        y = mvnpdf(x,mu,Sigma2);
        h = plot(x, y, 'color', colour);
        
    end
    
    if holdax
        hold off
    end
%==========================================================================