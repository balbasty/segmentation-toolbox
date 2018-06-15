function varargout = update_gmm(varargin)
% Update VB-GMM parameters (gmm.po) and mixture weights (wp and mg)
% FORMAT varargout = update_gmm(varargin)
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Read function input
%--------------------------------------------------------------------------
ll       = varargin{1};
llr      = varargin{2};
llrb     = varargin{3};
buf      = varargin{4};
mg       = varargin{5};
gmm      = varargin{6};
wp       = varargin{7};
lkp      = varargin{8};
wp_reg   = varargin{9};
iter     = varargin{10};
tol1     = varargin{11};
nm       = varargin{12};
nitgmm   = varargin{13};
do_wp    = varargin{14};
fig      = varargin{15};
L        = varargin{16};
print_ll = varargin{17};
wp_l     = varargin{18};
do_mg    = varargin{19};
resp     = varargin{20};
mrf      = varargin{21};
save_resp = varargin{22};

for subit=1:nitgmm
    oll = ll;
    ll  = llrb + llr;
    
    % Compute responsibilities and moments
    [mom,dll,mrf,mgm,resp] = compute_moments(buf,lkp,mg,gmm,wp,wp_l,resp,mrf);        
    ll                     = ll + dll;     
    
    % Add up 0:th moment
    s0 = 0;
    for i=2:numel(mom), s0 = s0 + mom(i).s0; end
    
    % Get total number of voxels
    nvox = 0;
    for z=1:numel(buf)
        nvox = nvox + buf(z).Nm;
    end
    
    if do_wp
        % Update tissue weights
        wp = update_wp(lkp,s0,mgm,nvox,wp_reg);
    end
    
    if do_mg
        % Update Gaussian weights
        mg = update_mg(lkp,s0);
    end
    
    % Update GMM posteriors
    [llvb,gmm] = spm_VBGaussiansFromSuffStats(mom,gmm);
    ll         = ll + sum(sum(llvb));  

    my_fprintf('MOG:\t%g\t%g\t%g\n',ll,llr,llrb,print_ll);
    L{1}(end + 1) = ll;

    if subit>1 || iter>1
        debug_view('convergence',fig{4},lkp,buf,L);
    end
    
    if subit>1 && ll-oll<tol1*nm
        % Improvement is small, so go to next step
        break;
    end
end

% % Save final responsibilities
ll                   = llrb + llr;
[mom,dll,mrf,~,resp] = compute_moments(buf,lkp,mg,gmm,wp,wp_l,resp,mrf,save_resp);        
ll                   = ll + dll;     
   
llvb = spm_VBGaussiansFromSuffStats(mom,gmm);
ll   = ll + sum(sum(llvb));  
    
my_fprintf('MOG:\t%g\t%g\t%g\n',ll,llr,llrb,print_ll);
L{1}(end + 1) = ll;
debug_view('convergence',fig{4},lkp,buf,L);

% Write function output
%--------------------------------------------------------------------------
varargout{1} = ll;
varargout{2} = mg;
varargout{3} = gmm;
varargout{4} = wp;
varargout{5} = L;
varargout{6} = mrf;
varargout{7} = mom;
varargout{8} = resp;
%==========================================================================