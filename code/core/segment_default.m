function pars = segment_default(pars,test_level)
if nargin<2, test_level = 0; end

if isstring(pars) || ischar(pars)
    % Parameters are in a JSON-file
    pars = parse_json(pars);
end

% General parameters
%--------------------------------------------------------------------------
if ~isfield(pars,'dir_output')
    pars.dir_output = './temp-data';
end
if test_level
    pars.dir_output = [pars.dir_output '_test'];
end
if exist(pars.dir_output,'dir'), rmdir(pars.dir_output,'s'); end; mkdir(pars.dir_output);

if ~isfield(pars,'dir_local')
    pars.dir_local = '';
end
if ~isfield(pars,'address_local')
    pars.address_local = '';
end
if ~isfield(pars,'dir_template')
    pars.dir_template = fullfile(pars.dir_output,'model');
end   
if ~isfield(pars,'use_2d_data')
    pars.use_2d_data = false;
end   

% Template specific parameters
%--------------------------------------------------------------------------
if ~isfield(pars,'do_template') 
    pars.do_template = true; 
end 
if pars.use_2d_data && ~pars.do_template
    if ~isfield(pars,'pth_template_2d')
        error('~isfield(pars,''pth_template_2d'')')
    else
        if isempty(pars.pth_template_2d)
            error('isempty(pars.pth_template_2d)')
        end
        pars.pth_template = pars.pth_template_2d;
    end
end
if ~isfield(pars,'pth_template')
    pars.pth_template = '';
end
if isfield(pars,'pth_template')
    if ~isempty(pars.pth_template)
        V1     = spm_vol(pars.pth_template); 
        pars.K = numel(V1);
    end
end    

if ~isfield(pars,'K')
    pars.K = 6;          
end
has_ct = inspect_ct_data(pars);
if has_ct && isempty(pars.pth_template)
    pars.K = 10;
end

if ~isfield(pars,'niter')
    pars.niter = 30;
end
if ~isfield(pars,'tol')
    pars.tol = 1e-4;
end
if ~isfield(pars,'verbose')
    pars.verbose = 4;
end
if ~isfield(pars,'vx_tpm')
    pars.vx_tpm = 1.5;
end
if ~isfield(pars,'sparam0')
    pars.sparam0 = 60;
end
if ~isfield(pars,'sparam1')
    pars.sparam1 = 20;
end
if ~isfield(pars,'uniform')
    pars.uniform = true;
end
if ~isfield(pars,'dt')
    pars.dt = [spm_type('float32') spm_platform('bigend')];
end
if ~isfield(pars,'sum_temp_der')
    pars.sum_temp_der = false;
end
if ~isfield(pars,'crop_template')
    pars.crop_template = 0;
end

% Data-set specific parameters (m=1,...,M)
%--------------------------------------------------------------------------
if ~isfield(pars,'dat')
    error('pars.dat needs to be defined!'); 
end

M  = numel(pars.dat);
S0 = 0;
for m=1:M    

    if ~pars.use_2d_data
        if ~isfield(pars.dat{m},'dir_data')
            error('~isfield(pars.dat,''dir_data'')')
        else
            if isempty(pars.dat{m}.dir_data)
                error('isempty(pars.dat.dir_data)')
            end
        end
    end
    
    if pars.use_2d_data
        if ~isfield(pars.dat{m},'dir_data_2d')
            error('~isfield(pars.dat,''dir_data_2d'')')
        else
            if isempty(pars.dat{m}.dir_data_2d)
                error('isempty(pars.dat.dir_data_2d)')
            end
            pars.dat{m}.dir_data = pars.dat{m}.dir_data_2d;
        end
    end

    % General parameters
    %----------------------------------------------------------------------    
    if ~isfield(pars.dat{m},'S')
        pars.dat{m}.S = Inf;        
    end    
    if     test_level==3, pars.dat{m}.S = min(16,pars.dat{m}.S);
    elseif test_level==2, pars.dat{m}.S = min(8,pars.dat{m}.S);   
    elseif test_level==1, pars.dat{m}.S = min(4,pars.dat{m}.S);   
    end 
    S0 = S0 + pars.dat{m}.S;
    
    if ~isfield(pars.dat{m},'modality')
        pars.dat{m}.modality = 'MRI';
    end
    if ~isfield(pars.dat{m},'healthy')
        pars.dat{m}.healthy = true;
    end    
    if ~isfield(pars.dat{m},'fwhm')
        pars.dat{m}.fwhm = 1;
    end
    if ~isfield(pars.dat{m},'est_fwhm')
        pars.dat{m}.est_fwhm = false;
    end        
    if ~isfield(pars.dat{m},'print_subj_info')
        if test_level==1
            pars.dat{m}.print_subj_info = true;
        else
            pars.dat{m}.print_subj_info = false;
        end    
    end   
    if ~isfield(pars.dat{m},'Affine')
        pars.dat{m}.Affine = eye(4);
    end
    
    % MI affine parameters
    %----------------------------------------------------------------------
    if ~isfield(pars.dat{m},'maff')
        pars.dat{m}.maff = struct;
    end
    if ~isfield(pars.dat{m}.maff,'do_maff')
        pars.dat{m}.maff.do_maff = true;
    end  
    if ~isfield(pars.dat{m}.maff,'affreg')
        pars.dat{m}.maff.affreg = 'mni';
    end
    if ~isfield(pars.dat{m}.maff,'maff_done')
        pars.dat{m}.maff.maff_done = false;
    end
    
    % Segmentation parameters
    %----------------------------------------------------------------------
    if ~isfield(pars.dat{m},'segment')
        pars.dat{m}.segment = struct;
    end
    if ~isfield(pars.dat{m}.segment,'gmm')
        pars.dat{m}.segment.gmm = struct;
    end        
    if ~isfield(pars.dat{m}.segment,'lkp')
        pars.dat{m}.segment.lkp = struct;
    end
    
    if isfield(pars.dat{m}.segment,'pth_prior')
        tmp                        = load(pars.dat{m}.segment.pth_prior);
        pars.dat{m}.segment.gmm.pr = tmp.pr;
        
        if isfield(pars.dat{m}.segment.gmm.pr,'part')
            pars.dat{m}.segment.lkp.part = pars.dat{m}.segment.gmm.pr.part;
        end                
    end        
        
    pars.dat{m}.segment.lkp.keep = 1:pars.K;    
    
    if ~isfield(pars.dat{m}.segment.lkp,'rem')
        pars.dat{m}.segment.lkp.rem = [];
    else
        msk                          = ismember(pars.dat{m}.segment.lkp.keep,pars.dat{m}.segment.lkp.rem);
        pars.dat{m}.segment.lkp.keep = pars.dat{m}.segment.lkp.keep(~msk);
    end    
    if ~isfield(pars.dat{m}.segment,'nlkp')
        pars.dat{m}.segment.nlkp = 1;
    end
    if ~isfield(pars.dat{m}.segment.lkp,'part')
        pars.dat{m}.segment.lkp.part = reshape(repelem(1:pars.K,1,pars.dat{m}.segment.nlkp),1,[]);
    end   
    if ~isfield(pars.dat{m}.segment.lkp,'lab')
        pars.dat{m}.segment.lkp.lab = [];
        pars.dat{m}.segment.wp_l    = 0;
    end    
    if ~isfield(pars.dat{m}.segment,'wp_l')
        pars.dat{m}.segment.wp_l = 0.8;
    end 
    if isempty(pars.dat{m}.segment.lkp.lab)
        pars.dat{m}.segment.lkp.lab = [];
    end
    if ~isfield(pars.dat{m}.segment,'samp')
        pars.dat{m}.segment.samp = 1;
    end    
    if pars.use_2d_data
        pars.dat{m}.segment.samp = 1;
    end
    if ~isfield(pars.dat{m}.segment,'do_bf')
        pars.dat{m}.segment.do_bf = true;
    end      
    if ~isfield(pars.dat{m}.segment,'do_def')
        pars.dat{m}.segment.do_def = true;
    end   
    if ~isfield(pars.dat{m}.segment,'do_wp')
        pars.dat{m}.segment.do_wp = true;
    end 
    if ~isfield(pars.dat{m}.segment,'do_mg') 
        pars.dat{m}.segment.do_mg = true; 
    end     
    if ~isfield(pars.dat{m}.segment,'do_missing_data')
        pars.dat{m}.segment.do_missing_data = false;
    end       
    if ~isfield(pars.dat{m}.segment,'kmeans_dist')
        pars.dat{m}.segment.kmeans_dist = 'cityblock';
    end   
    if ~isfield(pars.dat{m}.segment,'class_ix')
        pars.dat{m}.segment.class_ix = [];
    end   
    if ~isfield(pars.dat{m}.segment,'init_clust')
        pars.dat{m}.segment.init_clust = 'magnitude';
    end    
    if ~isfield(pars.dat{m}.segment,'print_ll')
        if test_level==1
            pars.dat{m}.segment.print_ll = true;
        else
            pars.dat{m}.segment.print_ll = false;
        end
    end   
    if ~isfield(pars.dat{m}.segment,'verbose')
        if test_level==1
            pars.dat{m}.segment.verbose = true;
        else
            pars.dat{m}.segment.verbose = false;
        end      
    end 
    if ~isfield(pars.dat{m}.segment,'reg')
        pars.dat{m}.segment.reg = [0 0.001 0.5 0.05 0.2]*0.1;
    end
    if ~isfield(pars.dat{m}.segment,'biasreg')
        pars.dat{m}.segment.biasreg = 1e-3;
    end
    if ~isfield(pars.dat{m}.segment,'biasfwhm')
        pars.dat{m}.segment.biasfwhm = 100;
    end    
    if ~isfield(pars.dat{m}.segment,'nsubit_bf')
        pars.dat{m}.segment.nsubit_bf = 6;
    end    
    if ~isfield(pars.dat{m}.segment,'nsubit_def')
        pars.dat{m}.segment.nsubit_def = 3;
    end        
    if ~isfield(pars.dat{m}.segment,'nitgmm')
        pars.dat{m}.segment.nitgmm = 20;
    end    
    if ~isfield(pars.dat{m}.segment,'niter')
        pars.dat{m}.segment.niter = 30;
    end    
    if ~isfield(pars.dat{m}.segment,'tol1')
        pars.dat{m}.segment.tol1 = 1e-4;
    end    
    if ~isfield(pars.dat{m}.segment,'wp_reg')
        pars.dat{m}.segment.wp_reg = 'single-subject';
    end        
    if ~isfield(pars.dat{m}.segment,'mg')
        pars.dat{m}.segment.mg = lkppart_to_mg(pars.dat{m}.segment.lkp.part);
    end    
    if ~isfield(pars.dat{m}.segment,'constr_inthp')
        pars.dat{m}.segment.constr_inthp = false;
    end    
    
    % MRF specific
    if ~isfield(pars.dat{m}.segment,'mrf')
        pars.dat{m}.segment.mrf = struct;
    end        
    if ~isfield(pars.dat{m}.segment.mrf,'do_mrf')
        pars.dat{m}.segment.mrf.do_mrf = false;
    end 
    if ~isfield(pars.dat{m}.segment.mrf,'ml')
        pars.dat{m}.segment.mrf.ml = true;
    end  
    if ~isfield(pars.dat{m}.segment.mrf,'alpha')        
        pars.dat{m}.segment.mrf.alpha = 1e5;
    end  
    if ~isfield(pars.dat{m}.segment.mrf,'val_diag')        
        pars.dat{m}.segment.mrf.val_diag = 0.8;
    end  
    if ~isfield(pars.dat{m}.segment.mrf,'Upsilon')        
        pars.dat{m}.segment.mrf.Upsilon = [];
    end      
    if ~isfield(pars.dat{m}.segment.mrf,'ElnUpsilon')        
        pars.dat{m}.segment.mrf.ElnUpsilon = [];
    end  
    if ~isfield(pars.dat{m}.segment.mrf,'update_Upsilon')        
        pars.dat{m}.segment.mrf.update_Upsilon = false;
    end  
    if ~isfield(pars.dat{m}.segment.mrf,'KK')        
        pars.dat{m}.segment.mrf.KK = [];
    end 
    
    % Push resps parameters
    %----------------------------------------------------------------------
    if ~isfield(pars.dat{m},'push_resp')
        pars.dat{m}.push_resp = struct;
    end
    if ~isfield(pars.dat{m}.push_resp,'do_push_resp')
        pars.dat{m}.push_resp.do_push_resp = false;
    end   
    
    % Write results parameters
    %----------------------------------------------------------------------
    if ~isfield(pars.dat{m},'write_res')
        pars.dat{m}.write_res = struct;
    end
    if ~isfield(pars.dat{m}.write_res,'do_write_res')
        pars.dat{m}.write_res.do_write_res = true;
    end       
    if ~isfield(pars.dat{m}.write_res,'bb')
        pars.dat{m}.write_res.bb = NaN(2,3);
    end    
    if ~isfield(pars.dat{m}.write_res,'vox')
        pars.dat{m}.write_res.vox = NaN;
    end
    if ~isfield(pars.dat{m}.write_res,'cleanup_gwc')
        pars.dat{m}.write_res.cleanup_gwc = false;
    end
    if ~isfield(pars.dat{m}.write_res,'cleanup_lesion')
        pars.dat{m}.write_res.cleanup_lesion = false;
    end    
    if ~isfield(pars.dat{m}.write_res,'mrf')
        pars.dat{m}.write_res.mrf = 2;
    end
    if pars.dat{m}.segment.mrf.do_mrf
        pars.dat{m}.write_res.mrf = 0;
    end    
    if ~isfield(pars.dat{m}.write_res,'write_tc')
        pars.dat{m}.write_res.write_tc = true(pars.K,4);                
    end    
    if ~isfield(pars.dat{m}.write_res,'write_bf')
        pars.dat{m}.write_res.write_bf = true(1,2);
    end    
    if ~isfield(pars.dat{m}.write_res,'write_df')
        pars.dat{m}.write_res.write_df = true(1,2);
    end           
    if isfield(pars.dat{m}.write_res,'G')
        G                       = pars.dat{m}.write_res.G;
        pars.dat{m}.write_res.G = single(G);                
    end           
    
    pars.dat{m}.segment.do_bf0  = pars.dat{m}.segment.do_bf; 
    pars.dat{m}.segment.do_def0 = pars.dat{m}.segment.do_def; 
    pars.dat{m}.segment.do_wp0  = pars.dat{m}.segment.do_wp; 
    pars.dat{m}.segment.do_mg0  = pars.dat{m}.segment.do_mg;            
    pars.dat{m}.segment.reg0    = pars.dat{m}.segment.reg;                
    pars.dat{m}.segment.mrf.update_Upsilon0 = pars.dat{m}.segment.mrf.update_Upsilon;    
end

if S0==1
    pars.do_template = false;
end

if ~pars.do_template 
    pars.niter = 1;
end

if pars.do_template 
    % Building template 
    for m=1:M               
           pars.dat{m}.push_resp.do_push_resp = true; 
           pars.dat{m}.write_res.do_write_res = false;            

           pars.dat{m}.segment.wp_reg = 'build-template'; 
           pars.dat{m}.segment.do_def = false; 
           pars.dat{m}.segment.do_bf  = false;        
           pars.dat{m}.segment.do_wp  = false; 
           pars.dat{m}.segment.do_mg  = false; 
           pars.dat{m}.segment.niter  = 1; 
           pars.dat{m}.segment.nitgmm = 1;
           pars.dat{m}.segment.nsubit_def = 1;
           pars.dat{m}.segment.nsubit_bf  = 1;
           
           pars.dat{m}.segment.tol1 = 1e-5;
           
           pars.dat{m}.segment.mrf.update_Upsilon = false;
    
           pars.dat{m}.write_res.mrf = 0;
    end
end 
%==========================================================================

%==========================================================================
function pars = append_dir(pars)
f                 = pars.dir_output;
pars.dir_output   = fullfile(f,['tmp-' pars.name]);
pars.dir_template = fullfile(f,['res-' pars.name]);
%==========================================================================

%==========================================================================
function pars = parse_json(pth)
pars = spm_jsonread(pth);
M    = numel(pars.dat);

if isstruct(pars.dat)
    % Convert dat from struct-array to cell-array    
    opars    = pars;
    pars     = rmfield(pars,'dat');    
    pars.dat = cell(1,M);
    for m=1:M
        pars.dat{m} = opars.dat(m);
    end
end

% Some cleaning up
for m=1:M
    if isfield(pars.dat{m},'S')
        if isempty(pars.dat{m}.S)
            pars.dat{m}.S = Inf;
        end
        
        if ischar(pars.dat{m}.S)
            if strcmpi(pars.dat{m}.S,'Inf')
                pars.dat{m}.S = Inf;
            end
        end
    end
end
%==========================================================================

%==========================================================================
function [has_ct,flag_healthy_ct] = inspect_ct_data(pars)
has_ct      = false;
healthy = [];
M1      = 0;
for m=1:numel(pars.dat)
    if strcmp(pars.dat{m}.modality,'CT')
        has_ct  = true;
        if isfield(pars.dat{m},'healthy')
            healthy = [healthy pars.dat{m}.healthy];
        else
            healthy = [healthy true];
        end
        M1      = M1 + 1;
    end
end
if sum(healthy)==M1
    % All healthy
    flag_healthy_ct = 0;
elseif sum(healthy)==0
    % All lesion
    flag_healthy_ct = 1;
else
    % Mix of healthy and lesion        
    flag_healthy_ct = 2;
end
%==========================================================================

%==========================================================================
function mg = lkppart_to_mg(lkppart)
Kb = max(lkppart);
K  = numel(lkppart);
mg = ones(1,K);
for k1=1:Kb
    k     = find(lkppart==k1);
    mg(k) = mg(k)/numel(k);
end
%==========================================================================