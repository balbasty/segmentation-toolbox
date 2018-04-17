function varargout = update_gmm(varargin)
% Update VB-GMM parameters (gmm.po) and mixture weights (wp and mg)
% FORMAT varargout = update_gmm(varargin)
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Read function input
%--------------------------------------------------------------------------
ll         = varargin{1};
llr        = varargin{2};
llrb       = varargin{3};
buf        = varargin{4};
mg         = varargin{5};
gmm        = varargin{6};
wp         = varargin{7};
lkp        = varargin{8};
wp_reg     = varargin{9};
mix_wp_reg = varargin{10};
iter       = varargin{11};
tol1       = varargin{12};
nm         = varargin{13};
nitgmm     = varargin{14};
do_wp      = varargin{15};
fig        = varargin{16};
L          = varargin{17};
print_ll   = varargin{18};
wp_lab     = varargin{19};

for subit=1:nitgmm
    oll = ll;
    ll  = llrb + llr;
    
    % Compute responsibilities and moments
    [mom,dll,mgm] = compute_moments(buf,lkp,mg,gmm,wp,wp_lab);        
    ll            = ll + dll;     
    
    % Add up 0:th moment
    s0 = 0;
    for i=2:numel(mom), s0 = s0 + mom(i).s0; end
    
    nvox = 0;
    for z=1:numel(buf)
        nvox = nvox + buf(z).Nm;
    end
    
    if do_wp
        % Update tissue weights
        wp = update_wp(lkp,s0,mgm,nvox,mix_wp_reg,wp_reg);
    end
    
    % Update mixing proportions
    mg = update_mg(lkp,s0);
        
    % Update means and variances
    [dll,gmm] = spm_VBGaussiansFromSuffStats(mom,gmm);
    ll        = ll + sum(sum(dll));  

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

% Write function output
%--------------------------------------------------------------------------
varargout{1} = ll;
varargout{2} = mg;
varargout{3} = gmm;
varargout{4} = wp;
varargout{5} = L;
%==========================================================================