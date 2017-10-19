function [X,j] = CG_img_solver(LHS,RHS,X0,opts)
% Conjugate gradient solver
% FORMAT [X,j] = CG_img_solver(LHS,RHS,X0,opts)
%
% LHS - Function handle containing left-hand side of linear system
%
% RHS - 3D matrix containing right hand side of linear systen
%
% X0 - 3D matrix containing initial guess [default: zeros(size(RHS))]
%
% opts - a structure containing various options.  The fields are:
%
%        maxiter - maximum number of iterations of algorithm [default: 100]
%
%        reltol - convergence tolerance [default: 1e-4]
%                 The algorithm has converged when the relative difference in 
%                 objective function, between two subsequent iterations, is
%                 below this value.
%
%        verbose - verbose level [default: 0]
%                  0 - No verbose
%                  1 - Show progress bar
%                  2 - (1) and plot convergence properties at each
%                  iteration
%                  3 - (1,2) and show reconstruction at each iteration
%                  4 - (1,2,3) and show reconstruction in fullscreen
%
%        Xref - reference solution [default: {}]
%               If a reference solution is available, for example when
%               using simulated data, the SNR can be computed between the
%               known reference and the reconstructed image. This SNR is
%               then plotted if verbose > 0.
%
% Xrec - cell array of reconstructed/recovered images
%
% j - total number of iterations required
%__________________________________________________________________________

  if nargin < 3,  X0 = zeros(size(RHS),'single'); end
  if isempty(X0), X0 = zeros(size(RHS),'single'); end
  if nargin < 4,  opts = []; end
  
  % -----------------------------------------------------------------------
  % Options
  maxiter = set_opts(opts,'maxiter', 200,1,1000);  
  objfun  = set_opts(opts,'objfun',    []); 
  reltol  = set_opts(opts,'reltol',  1e-6,1e-8,1);
  verbose = set_opts(opts,'verbose', 0);  
  Xref    = set_opts(opts,'Xref',    []);   
  
  % -----------------------------------------------------------------------
  % Initilisation  
  X = X0;
  clear X0
  
  normRHS = sqrt(sum(RHS(:).*RHS(:))); % Norm of RHS
  R       = RHS - LHS(X);              % Residual RHS - LHS(x)
  normR   = sum(R(:).*R(:));           % R'R
  P       = R;                         % Initial conjugate directions P
  beta    = 0;                         % Initial search direction for new P
  if isa(RHS,'gpuArray')    
    beta = gpuArray(beta);
    if verbose > 0, verbose = 0; end
  end
  
  % Verbose stuff
  if verbose > 1, F(2) = figure(2002); end
  if verbose > 2, set(gcf,'units','normalized','outerposition',[0 0 1 1]); end 
  if verbose > 0,
    history = struct();
    F(1) = figure(2001); 
    H = gobjects(3,1);
  end
  
  objval = [];
  
  j = 1;
  while sqrt(normR) > reltol*normRHS,
    % ---------------------------------------------------------------------
    % Calculate conjugate directions P which defines the direction of descent
    P = R + beta*P;

    % ---------------------------------------------------------------------
    % Finds the step size of the conj. gradient descent
    AtAP  = LHS(P);
    alpha = normR / sum(P(:).*AtAP(:));

    % ---------------------------------------------------------------------
    % Perform conj. gradient descent, obtaining updated X and R, using the calculated
    % P and alpha
    X = X + alpha *P; 
    R  = R - alpha*AtAP;

    % ---------------------------------------------------------------------
    % Finds the step size for updating P
    RtRp  = normR;
    normR = sum(R(:).*R(:));   
    beta  = normR / RtRp;            
    
    % Verbose stuff          
    if verbose > 1, show_fig; end
    if verbose > 0
      if ~isempty(objfun)
        objval(j) = objfun(X);
      end
      plot_convergence; 
    end     
    
    % Check if converged
    if j >= maxiter, break; end;                   
	
	j = j + 1; 
  end;
  
  %========================================================================
  %-function plot_convergence
  %========================================================================
  function plot_convergence  
    history.rnorm(j)  = sqrt(normR);  
    if ~isempty(Xref), history.snr(j) = calc_snr(X,Xref); end

    if j == 1
      set(0, 'CurrentFigure', F(1));        
      subplot(131); H(1)=semilogy(0:numel(history.rnorm)-1, history.rnorm', '-');  title('rnorm'); grid on; axis tight;                  
      if ~isempty(Xref)
        subplot(132); H(2)=plot(0:size(history.snr,1)-1, history.snr, '-'); title('SNR'); grid on; axis tight;          
      end
      if ~isempty(objfun)
        subplot(133); H(3)=plot(0:size(objval,1)-1, objval, '-'); title('objval'); grid on; axis tight;          
      end
    else
      set(H(1),'XData',0:numel(history.rnorm)-1,'YData',history.rnorm);
      if ~isempty(Xref)
        set(H(2),'XData',0:numel(history.snr)-1,'YData',history.snr);
      end
      if ~isempty(objfun)
        set(H(3),'XData',0:numel(objval)-1,'YData',objval);
      end      
    end
    drawnow
  end
  
  %========================================================================
  %-function show_fig
  %========================================================================
  function show_fig   
    set(0, 'CurrentFigure', F(2));
    if size((RHS),3) > 1        
      img1 = X(:,:,floor(size(X,3)/2));
      img2 = squeeze(X(floor(size(X,1)/2),:,:));
      img3 = squeeze(X(:,floor(size(X,2)/2),:));          
      subplot(131); imagesc(img1); colormap(gray); axis off; axis image;
      subplot(132); imagesc(img2); colormap(gray); axis off; axis image;
      title(['j=' num2str(j)])
      subplot(133); imagesc(img3); colormap(gray); axis off; axis image;
    else
      imagesc(X); colormap(gray); axis off; axis image;
      title(['j=' num2str(j)])
      pause(0.1)
    end      
    drawnow
  end  
end