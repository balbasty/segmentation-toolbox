function pars = segment_default(pars,test_level)
if nargin<2, test_level = 0; end

if isstring(pars) || ischar(pars)
    % Parameters are in a JSON-file
    pars = parse_json(pars);
end

% For CT data, initialises the GMM parameters by fitting a GMM to an
% accumulated histogram of image intensities.
pars = init_ct(pars,test_level);

% General parameters
%--------------------------------------------------------------------------
if ~isfield(pars,'name')
    pars.name = 'segmentation-toolbox';    
end
if test_level==1 || test_level==2, pars.name = 'test-local';
elseif test_level==3               pars.name = 'test-holly';   
end     
if ~isfield(pars,'dir_output')
    pars.dir_output = './temp-data/';
end
if ~isfield(pars,'dir_template')
    pars.dir_template = pars.dir_output;
end   
    
pars = append_dir(pars);

% Template specific parameters
%--------------------------------------------------------------------------
if ~isfield(pars,'K')
    pars.K = 6;
    if isfield(pars,'pth_template')
        if ~isempty(pars.pth_template)
            V1     = spm_vol(pars.pth_template); 
            pars.K = numel(V1);
        end
    end            
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
if ~isfield(pars,'pth_template')
    pars.pth_template = '';
end
if ~isfield(pars,'vx_tpm')
    pars.vx_tpm = 1.5;
end
if ~isfield(pars,'sparam')
    pars.sparam = [0 30 30];
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
if ~isfield(pars,'mrf')
    pars.mrf = false;
end

% Data-set specific parameters (m=1,...,M)
%--------------------------------------------------------------------------
if ~isfield(pars,'dat')
    error('pars.dat needs to be defined!'); 
end

M = numel(pars.dat);
for m=1:M
    if ~isfield(pars.dat{m},'dir_data')
        error('pars.dat.dir_data needs to be defined!'); 
    end

    % General parameters
    %----------------------------------------------------------------------    
    if ~isfield(pars.dat{m},'S')
        pars.dat{m}.S = Inf;        
    end    
    if     test_level==2 || test_level==3, pars.dat{m}.S = min(8,pars.dat{m}.S);
    elseif test_level==1,                  pars.dat{m}.S = 1;   
    end 
    if ~isfield(pars.dat{m},'modality')
        pars.dat{m}.modality = 'MRI';
    end
    if ~isfield(pars.dat{m},'trunc_ct')
        pars.dat{m}.trunc_ct = [];
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
    if ~isfield(pars.dat{m}.segment,'lkp')
        pars.dat{m}.segment.lkp = struct;
    end
    if ~isfield(pars.dat{m}.segment.lkp,'keep')
        pars.dat{m}.segment.lkp.keep = 1:pars.K;
    end    
    if ~isfield(pars.dat{m}.segment.lkp,'rem')
        pars.dat{m}.segment.lkp.rem = [];
    end    
    if ~isfield(pars.dat{m}.segment.lkp,'part')
        if isfield(pars.dat{m}.segment,'nlkp')
            pars.dat{m}.segment.lkp.part = reshape(repelem(1:pars.K,1,pars.dat{m}.segment.nlkp),1,[]);
        else
            pars.dat{m}.segment.lkp.part = 1:pars.K;
        end                       
    end     
    if ~isfield(pars.dat{m}.segment.lkp,'lab')
        pars.dat{m}.segment.lkp.lab = [];
    end
    if ~isfield(pars.dat{m}.segment,'pth_prior')
        pars.dat{m}.segment.pth_prior = '';
    end
    if ~isfield(pars.dat{m}.segment,'wp_lab')
        pars.dat{m}.segment.wp_lab = 0.5;
    end 
    if isempty(pars.dat{m}.segment.lkp.lab)
        pars.dat{m}.segment.wp_lab = 0;
    end
    if ~isfield(pars.dat{m}.segment,'samp')
        pars.dat{m}.segment.samp = 3;
    end    
    if ~isfield(pars.dat{m}.segment,'do_bf')
        pars.dat{m}.segment.do_bf = true;
    end      
    if ~isfield(pars.dat{m}.segment,'do_def')
     pars.dat{m}.segment.tol1 = 1e-4;   pars.dat{m}.segment.do_def = true;
    end   
    if ~isfield(pars.dat{m}.segment,'do_wp')
        pars.dat{m}.segment.do_wp = true;
    end 
    if ~isfield(pars.dat{m}.segment,'do_missing_data')
        pars.dat{m}.segment.do_missing_data = false;
    end       
    if ~isfield(pars.dat{m}.segment,'kmeans_dist')
        pars.dat{m}.segment.kmeans_dist = 'cityblock';
    end   
    if ~isfield(pars.dat{m}.segment,'kmeans_ix')
        pars.dat{m}.segment.kmeans_ix = [];
    end   
    if ~isfield(pars.dat{m}.segment,'init_clust')
        pars.dat{m}.segment.init_clust = '';
    end    
    if ~isempty(pars.dat{m}.segment.kmeans_ix)
        pars.dat{m}.segment.init_clust = 'user';
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
        pars.dat{m}.segment.biasreg = 1e-4;
    end
    if ~isfield(pars.dat{m}.segment,'biasfwhm')
        pars.dat{m}.segment.biasfwhm = 60;
    end    
    if ~isfield(pars.dat{m}.segment,'nsubit')
        pars.dat{m}.segment.nsubit = 8;
    end    
    if ~isfield(pars.dat{m}.segment,'nitgmm')
        pars.dat{m}.segment.nitgmm = 20;
    end    
    if ~isfield(pars.dat{m}.segment,'niter')
        pars.dat{m}.segment.niter = 30;
    end    
    if ~isfield(pars.dat{m}.segment,'tol1')
        pars.dat{m}.segment.tol1 = 0.5*1e-4;
    end    
    if ~isfield(pars.dat{m}.segment,'mix_wp_reg')
        pars.dat{m}.segment.mix_wp_reg = 0.9;
    end        
    if ~isfield(pars.dat{m}.segment,'mg')
        pars.dat{m}.segment.mg = lkppart_to_mg(pars.dat{m}.segment.lkp.part);
    end    
    if ~isfield(pars.dat{m}.segment,'constr_inthp')
        pars.dat{m}.segment.constr_inthp = false;
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
        pars.dat{m}.write_res.do_write_res = false;
    end       
    if ~isfield(pars.dat{m}.write_res,'bb')
        pars.dat{m}.write_res.bb = NaN(2,3);
    end    
    if ~isfield(pars.dat{m}.write_res,'vox')
        pars.dat{m}.write_res.vox = NaN;
    end
    if ~isfield(pars.dat{m}.write_res,'cleanup')
        pars.dat{m}.write_res.cleanup = false;
    end
    if ~isfield(pars.dat{m}.write_res,'mrf')
        pars.dat{m}.write_res.mrf = 1;
    end
    if ~isfield(pars.dat{m}.write_res,'write_tc')
        pars.dat{m}.write_res.write_tc      = false(pars.K,4);        
        pars.dat{m}.write_res.write_tc(:,1) = true;
    end    
    if ~isfield(pars.dat{m}.write_res,'write_bf')
        pars.dat{m}.write_res.write_bf = false(1,2);
    end    
    if ~isfield(pars.dat{m}.write_res,'write_df')
        pars.dat{m}.write_res.write_df = false(1,2);
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
    end
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