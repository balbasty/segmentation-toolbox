function [Xhat,j,h] = ADMM_img_solver(pm,Y,lam,opts)
% FORMAT [Xhat,j,h] = ADMM_img_solver(pm,Y,lam,opts)
%__________________________________________________________________________  
  r(1) = Inf;
  for b=1:opts.bregmitr % Bregman iterations
    if b > 1
      opts.X0 = Xhat;
    end
    
%     if opts.verbose > 0
%       fprintf('Deconvolution started... '); 
%       tin = tic;
%     end
    [Xhat,j,h] = run_algorithm(pm,Y,lam,opts);
%     if opts.verbose > 0
%       t = toc(tin);
%       fprintf('Completed in %5.1f s. and %i iters.\n',t,j)
%     end
    
%     if ~isempty(Xref)
%       for m=1:numel(Xref)
%         fprintf('b=%i, PSNR(Xhat,Xref)=%3.2f, SSIM(Xhat,Xref)=%3.3f\n',b,calc_psnr(Xhat{m},opts.Xref{m}),calc_ssim(Xhat{m},opts.Xref{m}))
%       end
%     end

    if opts.bregmitr > 1     
        Y0 = Y; 
        for m=1:numel(Y)
            R    =  Y{m} - pm.A(Xhat{m});
            Y{m} = Y0{m} + R;
            e(b) = sqrt(sum(sum(sum((Y0{m} - pm.A(Xhat{m})).^2))));
            if b > 1
                r(b) = abs((e(b - 1)*(1 + 10*eps) - e(b))/e(b));
            end      
        end

        if r(b) < 0.5*1e-2
            break
        end
    end
          
    if opts.verbose > 0 && opts.bregmitr > 1 && numel(Y) == 1 &&  size(Y{1},3) == 1
      figure(2003)
      subplot(221); 
      imagesc(R); colormap(gray); axis off; axis image; title('R')
      subplot(222)
      imagesc(Y{1}); colormap(gray); axis off; axis image; title('Y0 + R')
      subplot(223)
      plot(0:b - 1,e,'-')
      subplot(224)
      plot(0:b - 1,r,'-')
      drawnow
    end
    
  end

  if ~isempty(Xref)
    for m=1:numel(Xref)
      fprintf('b=%i, PSNR(Xhat,Xref)=%3.2f, SSIM(Xhat,Xref)=%3.3f\n',b,calc_psnr(Xhat{m},opts.Xref{m}),calc_ssim(Xhat{m},opts.Xref{m}))
    end
  end
    
  %========================================================================
  %-function [X,j,history] = run_algorithm(pm,Y,lam,opts)
  %========================================================================
  function [X,j,history] = run_algorithm(pm,Y,lam,opts) 
    if nargin < 4, opts = []; end

    % Set algorithm options
    lam2    = set_opts(opts,'lam2',   0,0,1e6);
    lstol   = set_opts(opts,'lstol',  0.5*1e-4,1e-8,1e-2);
    maxitr  = set_opts(opts,'maxitr', 100,1,1e5);
    nonneg  = set_opts(opts,'nonneg', true);
    reltol  = set_opts(opts,'reltol', 1e-4,1e-8,1);
    rho     = set_opts(opts,'rho',    1e1,1e-6,1e5);
    runpar  = set_opts(opts,'runpar', true);
    usegpu  = set_opts(opts,'usegpu', false);
    vx      = set_opts(opts,'vx',     [1 1 1]);
    verbose = set_opts(opts,'verbose',0);
    X0      = set_opts(opts,'X0',     {});
    Xref    = set_opts(opts,'Xref',   {});    

    lsopts.reltol = lstol;

    c   = (1 + lam2)^(-0.5);
    lam = c*lam;        
    
    % -----------------------------------------------------------------------
    % Initialise algorithm 
    M = numel(Y); % Number of modalities
    
%     if ~isa(Y{1},'single')
%       % Cast to single precision
%       for m=1:M
%         Y{m} = single(Y{m});
%       end
%     end

    if M == 1 || usegpu
      % Disable parallel execution
      runpar = false;
    end

    dim  = size(Y{1});
    ndim = 3; 
    if numel(dim) == 2 
      dim(3) = 1;
      ndim   = 2;
    end

    % Define Laplacian filter (in Fourier space)    
    L = zeros(dim,'single');  
    if dim(1)>=2
      tmp        = 1/(vx(1)^2);
      L(  1,1,1) = L(  1,1,1) + tmp*2;
      L(  2,1,1) = L(  2,1,1) - tmp;
      L(end,1,1) = L(end,1,1) - tmp;
    end
    if dim(2)>=2
      tmp        = 1/(vx(2)^2);
      L(1,  1,1) = L(1,  1,1) + tmp*2;
      L(1,  2,1) = L(1,  2,1) - tmp;
      L(1,end,1) = L(1,end,1) - tmp;
    end
    if dim(3)>=2
      tmp        = 1/(vx(3)^2);
      L(1,1,  1) = L(1,1,  1) + tmp*2;
      L(1,1,  2) = L(1,1,  2) - tmp;
      L(1,1,end) = L(1,1,end) - tmp;
    end
    L = fftn(L);

    % Initialize algorithm variables
    X = cell(M,1);
    W = cell(M,ndim);
    for m=1:M
      if isempty(X0)
        X{m} = zeros(dim,'single');   
      else
        X{m} = X0{m};   
      end
      for d=1:ndim   
        W{m,d} = zeros(dim,'single');
      end
    end
    U   = W; 
    RHS = cell(M,1); 
    clear X0
    
    if usegpu && hasGPU
      % Send data to the GPU
      L = gpuArray(L);
      for m=1:M
        X{m} = gpuArray(X{m});
        Y{m} = gpuArray(Y{m});
        for d=1:ndim 
          U{m,d} = gpuArray(U{m,d});
          W{m,d} = gpuArray(W{m,d});
        end
        lam(m) = gpuArray(lam(m));    
      end
      rho    = gpuArray(rho);

      verbose = 0;
    end

    % Verbose stuff
    history = struct();
    if ~isempty(Xref)
      Hpsnr = gobjects(2,1); 
    end
    if verbose > 1
      F(2) = figure(2002); 
    end
    if verbose > 2 
      set(gcf,'units','normalized','outerposition',[0 0 1 1]); 
    end
    if verbose > 0
      F(1) = figure(2001);
      H    = gobjects(7,1);
    end
      
    % -----------------------------------------------------------------------
    % Run algorithm
    Dx = cell(M,1); AtY = Dx;
    for m=1:M
      [Dx{m,1:ndim}] = grad(X{m},vx,'single'); 
      AtY{m}         = pm.At(Y{m});
    end

    for j=1:maxitr  
      if verbose > 0
        oU = U;
      end

      % ---------------------------------------------------------------------
      % Step 1: U-subproblem
      Unorm = eps;
      for m=1:M
        for d=1:ndim
          U{m,d} = lam(m)*Dx{m,d}+(1/rho)*W{m,d};
          Unorm  = Unorm + U{m,d}.^2;
        end
      end
      Unorm = sqrt(Unorm);
      scale = max(Unorm - 1/rho,0)./Unorm;
      clear Unorm
      scale(~isfinite(scale)) = 0;
      for m=1:M
        for d=1:ndim
          U{m,d} = U{m,d}.*scale;
        end
      end
      clear scale
      
      for m=1:M
        for d=1:ndim
          U{m,d} = sqrt(1 + lam2)*U{m,d};
        end
      end

      % ---------------------------------------------------------------------
      % Step 2: X-subproblem
      for m=1:M
        % Compute divergence
        if ndim==3
          DtU = lam(m)*dive(U{m,:},vx);
          DtW = lam(m)*dive(W{m,:},vx);
        else
          DtU = lam(m)*dive(U{m,:},[],vx);
          DtW = lam(m)*dive(W{m,:},[],vx);
        end

        RHS{m} = DtU - (1/rho)*DtW + (1/rho)*AtY{m};
        LHS{m} = @(X) lam(m)^2*laplacian(X,vx) + 1/rho*pm.AtA(X);
      end
      clear DtU DtW
        
      parfor (m=1:M,runpar*M)
        %X{m} = real(ifftn(fftn(RHS{m})./(L*(lambda(m)^2) + ((tau{m})/rho))));            
        X{m} = CG_img_solver(LHS{m},RHS{m},X{m},lsopts);

        if nonneg
          X{m}(X{m} < 0) = 0;
        end
      end
      clear RHS LHS
      
      % ---------------------------------------------------------------------
      % Step 3: W-subproblem
      for m=1:M
        [Dx{m,1:ndim}] = grad(X{m},vx,'single');
        for d=1:ndim
          W{m,d} = W{m,d} + rho*(lam(m)*Dx{m,d} - U{m,d});
        end
      end

      % ---------------------------------------------------------------------
      % Do miscellenaous

      % Compute objective value and relative change
      [history.objval(j+1),history.f(j+1),history.tv(j+1),history.l2(j+1)] = objval;
      history.relchg(j)   = abs((history.objval(j)*(1+10*eps) - history.objval(j+1))/history.objval(j+1));

      % Compute primal and dual residuals
      if verbose > 0
        ss1 = 0;
        ss2 = 0;
        for m=1:M
          tmp1 = cell(1,ndim);
          for d=1:ndim
            ss1     = ss1 + sum(sum(sum((lam(m)*Dx{m,d} - U{m,d}).^2)));
            tmp1{d} = U{m,d} - oU{m,d};
          end
          if ndim==3
            tmp = lam(m)*dive(tmp1{:},vx);
          else
            tmp = lam(m)*dive(tmp1{:},[],vx);
          end
          clear tmp1
          
          ss2 = ss2 + sum(sum(sum(tmp.^2)));
          clear tmp
        end
        history.rnorm(j) = sqrt(ss1);
        history.snorm(j) = sqrt(ss2)*rho;
      end       

      % Verbose stuff
      if verbose > 1, show_fig; end
      if verbose > 0, plot_convergence; end

      if abs(history.relchg(j)) < reltol && j > 10
        % Break if converged
        break; 
      end
    end    
    
    %========================================================================
    %-function [val,f,tv,l2] = objval
    %========================================================================
    function [val,f,tv,l2] = objval
      f  = 0;
      tv = 0;
      l2 = 0;
      for m=1:M
        f = f + sum(sum(sum((pm.A(X{m}) - Y{m}).^2)));
        for d=1:ndim
          tv = tv + lam(m)^2*(Dx{m,d}.^2);
        end
        if lam2 > 0
          lap = real(ifftn(fftn(lam(m)^2*X{m}).*L)); 
          l2  = l2 + X{m}(:)'*lap(:);
        end
      end
      tv  = sum(sum(sum(sqrt(tv))));
      val = 0.5*f + tv + lam2*l2;
    end

    %========================================================================
    %-function plot_convergence
    %========================================================================
    function plot_convergence
      if ~isempty(Xref)
        for m=1:M
          history.psnr(j,m) = calc_psnr(X{m},Xref{m});
        end
      end
      if j == 1
        set(0, 'CurrentFigure', F(1));
        subplot(231); H(1)=semilogy(0, 0, '-'); title('relchg'); grid on; axis tight;
        subplot(232); H(2)=semilogy(0:numel(history.rnorm)-1, history.rnorm', '-');  title('rnorm'); grid on; axis tight;
        subplot(233); H(3)=semilogy(0:numel(history.snorm)-1, history.snorm', '-'); title('snorm'); grid on; axis tight;
        subplot(234); H(4)=semilogy(0:numel(history.objval)-1, history.objval', '-'); title('objval'); grid on; axis tight;
        subplot(235); H(5:7)=semilogy(repmat(0:numel(history.objval)-1,3,1)', [history.f',history.tv',history.l2'], '-'); title('f (b), tv (r), l2 (y)'); grid on; axis tight;
        if ~isempty(Xref)
          subplot(236); Hpsnr=plot(0:size(history.psnr,1)-1, history.psnr, '-'); title('PSNR'); grid on; axis tight;
        end
      else
        set(H(1),'XData',0:numel(history.relchg)-2,'YData',history.relchg(2:end));
        set(H(2),'XData',0:numel(history.rnorm)-1,'YData',history.rnorm);
        set(H(3),'XData',0:numel(history.snorm)-1,'YData',history.snorm);
        for h=4:7
          if h == 4
            set(H(4),'XData',0:numel(history.objval)-1,'YData',history.objval);
          elseif h == 5
            set(H(5),'XData',0:numel(history.f)-1,'YData',history.f);
          elseif h == 6
            set(H(6),'XData',0:numel(history.tv)-1,'YData',history.tv);
          elseif h == 7
            set(H(7),'XData',0:numel(history.l2)-1,'YData',history.l2);          
          end
        end
        if ~isempty(Xref)
          for m=1:M
            set(Hpsnr(m),'XData',0:size(history.psnr,1)-1,'YData',history.psnr(:,m));
          end
        end
      end
      drawnow
    end

    %========================================================================
    %-function show_fig
    %========================================================================
    function show_fig
      set(0, 'CurrentFigure', F(2));
      if dim(3)>1
        for m=1:M
          img  = X{m};
          img1 = img(:,:,floor(size(img,3)/2));
          img2 = squeeze(img(floor(size(img,1)/2),:,:));
          img3 = squeeze(img(:,floor(size(img,2)/2),:));
          subplot(M,3,3*(m-1)+1); imagesc(img1); colormap(gray); axis off; axis image;
          subplot(M,3,3*(m-1)+2); imagesc(img2); colormap(gray); axis off; axis image;
          title(['j=' num2str(j)])
          subplot(M,3,3*(m-1)+3); imagesc(img3); colormap(gray); axis off; axis image;
        end
      else
        for m=1:M
          img = X{m};
          subplot(1,M,m); imagesc(img); colormap(gray); axis off; axis image;
          title(['j=' num2str(j)])
        end
        pause(0.1)
      end
      drawnow
    end
  end
end
%==========================================================================

%==========================================================================
function Div = dive(U,V,W,vd)  

    if nargin < 4, vd = ones(3,1); end

    if isa(U,'gpuArray')
        vd = gpuArray(vd);
    end

    if size(U,3) == 1
        Du = [-U(:,1), -diff(U(:,1:end-1),1,2), U(:,end-1)];
        Dv = [-V(1,:); -diff(V(1:end-1,:),1,1); V(end-1,:)];
        Div = Du./vd(2) + Dv./vd(1);
    else
        Du = cat(2, -U(:,1,:), -diff(U(:,1:end-1,:),1,2), U(:,end-1,:)); 
        Dv = cat(1, -V(1,:,:), -diff(V(1:end-1,:,:),1,1), V(end-1,:,:));
        Dw = cat(3, -W(:,:,1), -diff(W(:,:,1:end-1),1,3), W(:,:,end-1));
        Div = Du./vd(2) + Dv./vd(1) + Dw./vd(3);
    end
end
%==========================================================================

%==========================================================================
function [Dx,Dy,Dz] = grad(X, vd, precision) 

    if nargin < 2, vd = ones(3,1); end
    if nargin < 3, precision = 'single'; end

    if isa(X,'gpuArray')
    vd = gpuArray(vd);
    end

    if size(X,3) == 1
    if isa(X,'gpuArray')
      Dx = [diff(X,1,2), gpuArray(zeros(size(X,1),1,precision))]./vd(2);
      Dy = [diff(X,1,1); gpuArray(zeros(1,size(X,2),precision))]./vd(1);
      Dz = 0;
    else
      Dx = [diff(X,1,2), zeros(size(X,1),1,precision)]./vd(2);
      Dy = [diff(X,1,1); zeros(1,size(X,2),precision)]./vd(1);
      Dz = 0;
    end    
    else
    if isa(X,'gpuArray')
      Dx = cat(2, diff(X,1,2), gpuArray(zeros(size(X,1),1,size(X,3),precision)))./vd(2);
      Dy = cat(1, diff(X,1,1), gpuArray(zeros(1,size(X,2),size(X,3),precision)))./vd(1);
      Dz = cat(3, diff(X,1,3), gpuArray(zeros(size(X,1),size(X,2),1,precision)))./vd(3);  
    else
      Dx = cat(2, diff(X,1,2), zeros(size(X,1),1,size(X,3),precision))./vd(2);
      Dy = cat(1, diff(X,1,1), zeros(1,size(X,2),size(X,3),precision))./vd(1);
      Dz = cat(3, diff(X,1,3), zeros(size(X,1),size(X,2),1,precision))./vd(3);  
    end
    end
end
%==========================================================================

%==========================================================================
function L = laplacian(X,vd)

    if nargin < 2, vd = ones(3,1); end

    [dx,dy,dz] = grad(X,vd);
    L          = dive(dx,dy,dz,vd);
end
%==========================================================================