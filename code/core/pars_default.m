function pars = pars_default(pars,test_level)
if nargin<2, test_level = 0; end

% General parameters
%--------------------------------------------------------------------------
if ~isfield(pars,'name')
    pars.name = 'build-template';    
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
if ~isfield(pars,'do_segment')
    pars.do_segment = true;
end
    
pars = append_dir(pars);

% Template specific parameters
%--------------------------------------------------------------------------
if ~isfield(pars,'K')
    pars.K = 6;
end
if ~isfield(pars,'niter')
    pars.niter = 50;
end
if ~isfield(pars,'tol')
    pars.tol = 1e-4;
end
if ~isfield(pars,'verbose')
    pars.verbose = 3;
end
if ~isfield(pars,'pth_template')
    pars.pth_template = '';
end
if ~isfield(pars,'vx_tpm')
    pars.vx_tpm = 1.5;
end
if ~isfield(pars,'sparam')
    pars.sparam = [0.001 2 0];
end
if ~isfield(pars,'uniform')
    pars.uniform = true;
end
if ~isfield(pars,'dt')
    pars.dt = [spm_type('float32') spm_platform('bigend')];
end

% Data-set specific parameters
%--------------------------------------------------------------------------
if ~isfield(pars,'dat'), 
    error('pars.dat needs to be defined!'); 
end

M = numel(pars.dat);
for m=1:M
    if ~isfield(pars.dat{m},'dir_data'), 
        error('pars.dat.dir_data needs to be defined!'); 
    end

    if ~isfield(pars.dat{m},'do_preproc')
        pars.dat{m}.do_preproc = false;
    end
    if ~isfield(pars.dat{m},'S')
        pars.dat{m}.S = 8;
    end
    if test_level==2 || test_level==3, pars.dat{m}.S = 8;
    elseif test_level==1               pars.dat{m}.S = 1;   
    end 
    if ~isfield(pars.dat{m},'modality')
        pars.dat{m}.modality = 'MRI';
    end    
    
    % Pre-processing parameters
    %----------------------------------------------------------------------
    if ~isfield(pars.dat{m},'preproc')
        pars.dat{m}.preproc = struct;
    end
    if ~isfield(pars.dat{m}.preproc,'do_rem_corrupted')
        pars.dat{m}.preproc.do_rem_corrupted = false;
    end
    if ~isfield(pars.dat{m}.preproc,'tol_dist')
        pars.dat{m}.preproc.tol_dist = 4;
    end
    if ~isfield(pars.dat{m}.preproc,'tol_vx')
        pars.dat{m}.preproc.tol_vx = 5;
    end
    if ~isfield(pars.dat{m}.preproc,'verbose')
        pars.dat{m}.preproc.verbose = false;
    end
    if ~isfield(pars.dat{m}.preproc,'coreg_and_reslice')
        pars.dat{m}.preproc.coreg_and_reslice = true;
    end
    if ~isfield(pars.dat{m}.preproc,'do_reslice')
        pars.dat{m}.preproc.do_reslice = true;
    end
    if ~isfield(pars.dat{m}.preproc,'realign2mni')
        pars.dat{m}.preproc.realign2mni = true;
    end
    if ~isfield(pars.dat{m}.preproc,'crop')
        pars.dat{m}.preproc.crop = true;
    end
    if ~isfield(pars.dat{m}.preproc,'rem_neck')
        pars.dat{m}.preproc.rem_neck = true;
    end
    if ~isfield(pars.dat{m}.preproc,'skull_strip')
        pars.dat{m}.preproc.skull_strip = false;
    end
    
    % Segmentation parameters
    %----------------------------------------------------------------------
    if ~isfield(pars.dat{m},'segment')
        pars.dat{m}.segment = struct;
    end
    if ~isfield(pars.dat{m}.segment,'pth_prior')
        pars.dat{m}.segment.pth_prior = '';
    end
    if ~isfield(pars.dat{m}.segment,'dir_lab')
        pars.dat{m}.segment.dir_lab = '';
    end    
    if ~isfield(pars.dat{m}.segment,'K_rem')
        pars.dat{m}.segment.K_rem = [];
    end
    if ~isfield(pars.dat{m}.segment,'samp')
        pars.dat{m}.segment.samp = 3;
    end
    if ~isfield(pars.dat{m}.segment,'init_clust')
        if M>1
            pars.dat{m}.segment.init_clust = 'total';
        else
            pars.dat{m}.segment.init_clust = '';
        end
    end
    if ~isfield(pars.dat{m}.segment,'do_ml')
        pars.dat{m}.segment.do_ml = false;
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
    if ~isfield(pars.dat{m}.segment,'do_write_res')
        pars.dat{m}.segment.do_write_res = false;
    end   
    if ~isfield(pars.dat{m}.segment,'do_push_resp')
        pars.dat{m}.segment.do_push_resp = false;
    end   
    if ~isfield(pars.dat{m}.segment,'do_missing_data')
        pars.dat{m}.segment.do_missing_data = false;
    end       
    if ~isfield(pars.dat{m}.segment,'do_old_segment')
        pars.dat{m}.segment.do_old_segment = false;
    end   
    if ~isfield(pars.dat{m}.segment,'kmeans_dist')
        pars.dat{m}.segment.kmeans_dist = 'cityblock';
    end   
    if ~isfield(pars.dat{m}.segment,'print_ll')
        if test_level==1
            pars.dat{m}.segment.print_ll = true;
        else
            pars.dat{m}.segment.print_ll = false;
        end
    end   
    if ~isfield(pars.dat{m}.segment,'print_seg')
        if test_level==1
            pars.dat{m}.segment.print_seg = true;
        else
            pars.dat{m}.segment.print_seg = false;
        end        
    end   
    if ~isfield(pars.dat{m}.segment,'verbose')
        if test_level==1
            pars.dat{m}.segment.verbose = true;
        else
            pars.dat{m}.segment.verbose = false;
        end        
    end   
    if ~isfield(pars.dat{m}.segment,'trunc_ct')
        pars.dat{m}.segment.trunc_ct = [-Inf Inf];
    end   
    if ~isfield(pars.dat{m}.segment,'nlkp')
        pars.dat{m}.segment.nlkp = 1;
    end   
    if ~isfield(pars.dat{m}.segment,'lkp')
        pars.dat{m}.segment.lkp = reshape(repmat(1:pars.K,1,1),1,[]);
    end       
    if ~isfield(pars.dat{m}.segment,'bb')
        pars.dat{m}.segment.bb = NaN(2,3);
    end    
    if ~isfield(pars.dat{m}.segment,'vox')
        pars.dat{m}.segment.vox = NaN;
    end
    if ~isfield(pars.dat{m}.segment,'cleanup')
        pars.dat{m}.segment.cleanup = false;
    end
    if ~isfield(pars.dat{m}.segment,'mrf')
        pars.dat{m}.segment.mrf = 1;
    end
    if ~isfield(pars.dat{m}.segment,'affreg')
        pars.dat{m}.segment.affreg = 'mni';
    end
    if ~isfield(pars.dat{m}.segment,'reg')
        pars.dat{m}.segment.reg = [0 0.001 0.5 0.05 0.2]*0.1;
    end
    if ~isfield(pars.dat{m}.segment,'fwhm')
        pars.dat{m}.segment.fwhm = 1;
    end
    if ~isfield(pars.dat{m}.segment,'Affine')
        pars.dat{m}.segment.Affine = eye(4);
    end
    if ~isfield(pars.dat{m}.segment,'biasreg')
        pars.dat{m}.segment.biasreg = 1e-3*(1/5);
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
        pars.dat{m}.segment.tol1 = 1e-4;
    end    
    if ~isfield(pars.dat{m}.segment,'write_tc')
        pars.dat{m}.segment.write_tc        = true(pars.K,4);
        pars.dat{m}.segment.write_tc(:,2:4) = false;
    end    
    if ~isfield(pars.dat{m}.segment,'write_bf')
        pars.dat{m}.segment.write_bf = false(1,2);
    end    
    if ~isfield(pars.dat{m}.segment,'write_df')
        pars.dat{m}.segment.write_df = false(1,2);
    end    
end
%==========================================================================

%==========================================================================
function pars = append_dir(pars)
pth               = fileparts(pars.dir_output);
pars.dir_output   = fullfile(pth,['tmp-' pars.name]);
pth               = fileparts(pars.dir_template);
pars.dir_template = fullfile(pth,['res-' pars.name]);
%==========================================================================