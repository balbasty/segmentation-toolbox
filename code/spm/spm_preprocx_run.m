function spm_preprocx_run(obj,im,K)

%==========================================================================
% Make directories
%==========================================================================

if ~exist(obj.dir_data,'dir'), mkdir(obj.dir_data); end
obj.dir_res = fullfile(obj.dir_data,obj.dir_res); 
if exist(obj.dir_res,'dir'), rmdir(obj.dir_res,'s'); end; mkdir(obj.dir_res);

obj.dir_Twarp = fullfile(obj.dir_data,'Twarp');
if exist(obj.dir_Twarp,'dir'), rmdir(obj.dir_Twarp,'s'); end; mkdir(obj.dir_Twarp);

%==========================================================================
% Load and process image data
%==========================================================================

[V,labels] = load_and_process_images(obj,im);

%==========================================================================
% Initialise template
%==========================================================================

obj = init_TPM(obj,V,K);

%==========================================================================
% Initialise some debugging output
%==========================================================================

[obj.fig,fig_TPM,fig_L,obj.verbose] = create_fig(obj);

%==========================================================================
% Initialise algorithm i/o 
%==========================================================================

obj = get_obj(obj,V,K,im,labels);

%==========================================================================
% Run algorithm
%==========================================================================

fprintf('==============================================\n')   
fprintf('Algorithm started (')
fprintf(datestr(now,'mmmm dd, yyyy HH:MM:SS'))
fprintf(')\n')
fprintf('==============================================\n')   

nitermain = obj{1}{1}.nitermain;
tolmain   = obj{1}{1}.tolmain;
dopr      = obj{1}{1}.dopr;
use_mog   = ~isempty(obj{1}{1}.lkp);

L       = -Inf; % Lower bound of complete model
for iter=1:nitermain
    fprintf('iter=%d======================\n',iter);  
           
    obj = set_iter(obj,iter);       
    
    % Update subject specific parameters (clusters, bias, affine, deformations, template derivatives)
    [obj,L,munum,muden,Nm] = update_subjects(obj,L); 
     
    plot_objval(L,fig_L); 
        
    if obj{1}{1}.dotpm
        % Update the template 
        update_global_TPM(obj,munum,muden);
              
        show_TPM(obj,fig_TPM);
    end  
    
    if dopr && use_mog && obj{1}{1}.vb && iter>=3                                      
        % Update subject specific parameters (clusters, bias, deformations)
        [obj,L] = update_subjects(obj,L); 
              
        % Update intensity prior
        obj = update_global_prior(obj);
        
        plot_objval(L,fig_L);
    end
        
    % Check convergence----------------------------------------------------
    if ~(abs(L(end)-L(end - 1))>2*tolmain*Nm)        
        fprintf('==============================================\n')                        
        fprintf('Algorithm converged in %d iterations.\n',iter)                        
        fprintf('==============================================\n')                        
        
        save_mog_parameters(obj)        
        break;
    end
end
%==========================================================================

%==========================================================================
function [obj,munum,muden] = segment(obj)

% Load template
logtpm = spm_load_logpriors8(obj.pth_logTPM,obj.tiny,obj.deg,obj.uniform);

if obj.doaff && ((obj.use_tpm && obj.iter==1) || (~obj.use_tpm && obj.iter==3))
    % Affinely register template to 1st image (n=1)
    obj.Affine = spm_logmaff8(obj.image(1),obj.samp,(obj.fwhm+1)*16,obj.descrip,logtpm,obj.Affine,'mni');            
    obj.Affine = spm_logmaff8(obj.image(1),obj.samp, obj.fwhm,      obj.descrip,logtpm,obj.Affine,'mni');                                    
    obj.doaff  = 0;
end                 

if ~obj.dodef0
    % Do not estimate deformations
    obj.dodef = 0;
elseif obj.dodef0 && obj.use_tpm && obj.iter==1
    % If a template is used, start estimating deformations from 1st
    % iteration
    obj.dodef = 1;
elseif ~obj.use_tpm && ~obj.dotpm        
    % If no template, and no template update, do not estimate deformations
    obj.dodef = 0;
end

if ~obj.use_tpm && obj.iter==1 && obj.dotpm   
    % 1st iteration
    %----------------------------------------------------------------------    
    obj.dodef = 0;              % No deformation update for 1st iteration
elseif ~obj.use_tpm && obj.iter==2 && obj.dotpm                      
    % 2nd iteration    
    %----------------------------------------------------------------------
    if obj.clear_pars
        % Re-estimate cluster and bias field parameters based on template constructed after 1st iteration
        obj = clear_pars(obj);
    end
elseif ~obj.use_tpm && obj.iter==3 && obj.dotpm                      
    % 3rd iteration
    %----------------------------------------------------------------------    
    if obj.clear_pars
        % Re-estimate cluster and bias field parameters based on template constructed after 2nd iteration
        obj = clear_pars(obj);    
        
        obj.clear_pars = 0;
    end    
elseif ~obj.use_tpm && obj.iter==4 && obj.dotpm       
    % 4th iteration: start estimating deformations
    %----------------------------------------------------------------------
    obj.dodef = obj.dodef0; 
end         

if isfield(obj,'pth_Twarp')
    % To save memory, load deformation from file
    load(obj.pth_Twarp,'Twarp')
    obj.Twarp = Twarp;
end

% Run segmentation algorithm (spm_preprocx)
munum  = 0;
muden  = 0;
obj.status = 0; % All OK   
try       
    if nargout>=2
        [obj,munum,muden] = spm_preprocx(obj,logtpm);
    else
        obj               = spm_preprocx(obj,logtpm);
    end
catch ME
    % Error!
    obj.status = 1;
    str_output = ['Error in: ' obj.image(1).fname '\n' ME.message '\n'];
    fprintf(str_output)
end

if isfield(obj,'pth_Twarp')
    % To save memory, delete deformation from result structure
    Twarp = obj.Twarp;
    save(obj.pth_Twarp,'Twarp')
    obj   = rmfield(obj,'Twarp');
end
%==========================================================================

%==========================================================================
function save_mog_parameters(obj)
vb      = obj{1}{1}.vb;
dir_res = obj{1}{1}.dir_res;
iter    = obj{1}{1}.iter;

M       = numel(obj);
mog_par = cell(1,M);
for m=1:M
    S          = numel(obj{m});
    mog_par{m} = cell(1,S);
    for s=1:S        
        par.iter          = iter;
        par.wp            = obj{m}{s}.wp;
        par.mg            = obj{m}{s}.mg;
        if vb
            par.pr        = obj{m}{s}.pr;
            par.po        = obj{m}{s}.po;
            mog_par{m}{s} = par;
        else
            par.mn        = obj{m}{s}.mn;
            par.vr        = obj{m}{s}.vr;
            mog_par{m}{s} = par;
        end
   end
end
fname = fullfile(dir_res,'mog_par.mat');
save(fname,'mog_par')
%==========================================================================

%==========================================================================
function obj = update_global_prior(obj)
dir_res = obj{1}{1}.dir_res;

for m=1:numel(obj)
    pr  = update_intensity_prior(obj{m});
    S_m = numel(obj{m});
    for s=1:S_m
        obj{m}{s}.pr = pr;
    end      
    save(fullfile(dir_res,['pr_m' num2str(m) '.mat']),'pr'); 
end
%==========================================================================        

%==========================================================================
function update_global_TPM(obj,munum,muden)
pth_logTPM = obj{1}{1}.pth_logTPM;
fwhm_TPM   = obj{1}{1}.fwhm_TPM;
tiny       = obj{1}{1}.tiny;
mrf        = obj{1}{1}.mrf;
deg        = obj{1}{1}.deg;
dir_res    = obj{1}{1}.dir_res;

if obj{1}{1}.uniform
    % Template will no longer be uniform after update    
    M = numel(obj);
    for m=1:M
        S = numel(obj{m});
        for s=1:S
            obj{m}{s}.uniform = 0;
        end
    end
end

% Load template
Nii = nifti(pth_logTPM);
d   = size(Nii.dat(:,:,:,:));
dm  = d(1:3);
K   = d(4);

% Smooth template
[munum,muden] = smooth_template(munum,muden,dm,fwhm_TPM);

% Update template
logmu            = log(munum./muden + tiny);
logmu            = reshape(logmu',[dm K]);                                          
Nii.dat(:,:,:,:) = logmu;        
clear munum muden

if mrf
    % Use a MRF cleanup procedure
    logtpm = spm_load_logpriors8(pth_logTPM,tiny,deg,uniform);

    phi = double(identity(logtpm.d));
    b   = spm_sample_logpriors8(logtpm,phi(:,:,:,1),phi(:,:,:,2),phi(:,:,:,3));
    clear phi

    Q = zeros([logtpm.d(1:3),K],'single');
    for k=1:K
       Q(:,:,:,k) = b{k};
    end
    clear b

    P        = zeros([logtpm.d(1:3),K],'uint8');
    nmrf_its = 1;           
    T        = 1;
    G        = ones([K,1],'single')*T;
    vx2      = 1./single(sqrt(sum(logtpm.M(1:3,1:3).^2)));
    for i=1:nmrf_its
        spm_mrf(P,Q,G,vx2);
    end             
    clear Q logtpm

    P                = double(P)/255;
    Nii.dat(:,:,:,:) = log(P + eps*eps);
    clear P
end   
clear Nii logmu

% Convert log template to probability template and save to disk
logTPM2TPM(pth_logTPM,dir_res);
%==========================================================================

%==========================================================================
function obj = set_iter(obj,iter)
M = numel(obj);
for m=1:M
    S = numel(obj{m});
    for s=1:S
        obj{m}{s}.iter = iter;
    end
end
%==========================================================================

%==========================================================================
function show_TPM(obj,fig_TPM)
if ~isempty(fig_TPM) 
    pth_logTPM = obj{1}{1}.pth_logTPM;
    tiny       = obj{1}{1}.tiny;
    deg        = obj{1}{1}.deg;
    uniform    = obj{1}{1}.uniform;

    logtpm = spm_load_logpriors8(pth_logTPM,tiny,deg,uniform);

    phi = double(identity(logtpm.d));
    b   = spm_sample_logpriors8(logtpm,phi(:,:,:,1),phi(:,:,:,2),phi(:,:,:,3));
    clear phi

    K  = numel(b);
    K1 = floor(sqrt(K)); 
    K2 = ceil(K/K1); 
    set(0,'CurrentFigure',fig_TPM);                                        
    for i=1:K    
        subplot(K1,K2,i);
        imagesc(b{i}(:,:,floor(logtpm.d(3)/2) + 1)'); axis image xy off; colormap(gray);
    end 
    clear b
    drawnow
end
%==========================================================================

%==========================================================================
function obj = clear_pars(obj)
% Re-estimate cluster parameters based on template constructed after 1st iteration
if isfield(obj,'mg'), obj = rmfield(obj,'mg'); end
if isfield(obj,'wp'), obj = rmfield(obj,'wp'); end
if isfield(obj,'mn'), obj = rmfield(obj,'mn'); end
if isfield(obj,'vr'), obj = rmfield(obj,'vr'); end
if isfield(obj,'po'), obj = rmfield(obj,'po'); end
if isfield(obj,'pr'), obj = rmfield(obj,'pr'); end                        

% Re-estimate bias field as well
if isfield(obj,'Tbias'), obj = rmfield(obj,'Tbias'); end  
%==========================================================================     

%==========================================================================
function plot_objval(L,fig_L)
if ~isempty(fig_L)  
    set(0,'CurrentFigure',fig_L);                
    plot(0:numel(L) - 1,L,'b-','LineWidth',1);   hold on            
    plot(0:numel(L) - 1,L,'b.','markersize',10); hold off  
end
%==========================================================================

%==========================================================================
function [obj,L,munum,muden,Nm] = update_subjects(obj,L)
num_workers  = obj{1}{1}.num_workers;
run_on_holly = obj{1}{1}.run_on_holly;
iter         = obj{1}{1}.iter;

dotpm = nargout>=3;

munum = single(0); 
muden = single(0);  
ll    = 0;
Nm    = 0;
for m=1:numel(obj)
    if run_on_holly
        
    else
        obj_m = obj{m};
        S     = numel(obj_m);
        ix    = [];    
%         for s=1:S
        parfor (s=1:S,num_workers)
            if dotpm
                fprintf('iter=%d, m=%d, s=%d (TPM)\n',iter,m,s);  
                [obj_m{s},num,den] = segment(obj_m{s});
            else
                fprintf('iter=%d, m=%d, s=%d (prior)\n',iter,m,s);  
                obj_m{s} = segment(obj_m{s});
                num      = 0;
                den      = 0;
            end
            if obj_m{s}.status
                % spm_preprocx errored...remove subject
                ix    = [ix,s];
            else
                munum = munum + num;
                muden = muden + den;            
                ll    = ll + obj_m{s}.ll;
                Nm    = Nm + obj_m{s}.nm;
            end
        end
        obj_m(ix) = [];
        obj{m}    = obj_m;
    end
end

obj_empty = 1;
for m=1:numel(obj)
    if ~isempty(obj{m})
        obj_empty = 0;
    end
end
if obj_empty
    error('obj_empty')
end

L = [L,ll];
%==========================================================================

%==========================================================================
function obj = get_obj(obj0,V,K,im,labels)
M   = numel(V);
obj = cell(1,M);
for m=1:M    
    S      = numel(V{m});
    obj{m} = cell(1,S);    
    for s=1:S
        obj{m}{s}.image    = V{m}{s};
        N                  = numel(V{m}{s});
        obj{m}{s}.biasfwhm = obj0.biasfwhm*ones(1,N);
        obj{m}{s}.biasreg  = obj0.biasreg*ones(1,N);       

        obj{m}{s}.use_tpm = obj0.use_tpm;
        
        obj{m}{s}.vb         = obj0.vb;        
        obj{m}{s}.wp_reg     = obj0.wp_reg;
        obj{m}{s}.clear_pars = 1;
        if obj0.dotpm
            % One Gaussian per tissue for template construction
            obj{m}{s}.lkp = 1:K; 
        else
            if isscalar(obj0.lkp)
                obj{m}{s}.lkp = repelem(1:K,obj0.lkp);
            else
                obj{m}{s}.lkp = obj0.lkp;
            end
        end        

        obj{m}{s}.Affine  = eye(4);
        obj{m}{s}.reg     = obj0.rparam;
        obj{m}{s}.samp    = obj0.samp;
        obj{m}{s}.fwhm    = 0;
        obj{m}{s}.verbose = obj0.verbose;

        obj{m}{s}.dobias     = obj0.dobias;
        obj{m}{s}.dodef0     = obj0.dodef; 
        obj{m}{s}.dotpm      = obj0.dotpm;         
        obj{m}{s}.doaff      = obj0.doaff;
        obj{m}{s}.dopr       = obj0.dopr;
                
        obj{m}{s}.nitermain  = obj0.nitermain;
        obj{m}{s}.tolmain    = obj0.tolmain;       
        obj{m}{s}.tolseg     = obj0.tolseg;
        obj{m}{s}.niter      = obj0.niter;
        obj{m}{s}.niter1     = obj0.niter1;
        obj{m}{s}.nsubitmog  = obj0.nsubitmog;
        obj{m}{s}.nsubitbf   = obj0.nsubitbf;
        obj{m}{s}.nitdef     = obj0.nitdef;

        obj{m}{s}.descrip = im{m}{3};
        obj{m}{s}.healthy = im{m}{4}; 
        obj{m}{s}.labels  = labels{m}{s}; 
        
        obj{m}{s}.fig = obj0.fig; 
              
        obj{m}{s}.pth_logTPM = obj0.pth_logTPM;
        obj{m}{s}.fwhm_TPM   = obj0.fwhm_TPM;
        obj{m}{s}.tiny       = obj0.tiny;
        obj{m}{s}.mrf        = obj0.mrf;
        obj{m}{s}.deg        = obj0.deg;
        obj{m}{s}.dir_res    = obj0.dir_res;
        obj{m}{s}.uniform    = obj0.uniform;
        obj{m}{s}.iter       = 0;
        
        obj{m}{s}.num_workers  = obj0.num_workers;
        obj{m}{s}.run_on_holly = obj0.run_on_holly;
        
        if S>1
            % Allocate initial deformation field (only for multiple subjects) 
            [obj{m}{s}.pth_Twarp] = allocate_deformation(obj{m}{s},m,s,obj0.dir_Twarp);
        end
    end
end
%==========================================================================

%==========================================================================
function [fig,fig_TPM,fig_L,verbose] = create_fig(obj)
verbose     = obj.verbose;
num_workers = obj.num_workers;
figix       = obj.figix;
dotpm       = obj.dotpm;

fig = cell(4,1);

if verbose==2 && ~num_workers
    spm_figure('Create','Interactive');
    figure(figix)

    for i=1:size(fig,1)              
        fig{i} = figure(figix + i);
        clf(fig{i})
    end 
end  

if verbose && dotpm
    fig_L   = figure(figix + 5); 
    fig_TPM = figure(figix + 6);    
    if num_workers
        verbose = 0;
    end
else
    fig_TPM = [];
    fig_L   = [];
end

try
    distFig; 
catch
    warning('distFig not available')
end
drawnow
%==========================================================================

%==========================================================================
function logTPM2TPM(pth_logTPM,dir_write)
if nargin<2, dir_write = './data/results'; end

pth_TPM = 'TPM.nii';
pth_TPM = fullfile(dir_write,pth_TPM);

Nii = nifti(pth_logTPM);
img = Nii.dat(:,:,:,:);
d   = size(img);
dm  = d(1:3);
K   = d(4);

img = reshape(img,[prod(dm) K]);
img = safe_softmax(img);
img = reshape(img,d);

[dir_write,nam,ext] = fileparts(pth_TPM);

vols = cell(K,1);
for k=1:K    
    vols{k}       = fullfile(dir_write,[nam num2str(k) ext]);
    create_nii(vols{k},img(:,:,:,k),Nii.mat,16,nam);
end

matlabbatch{1}.spm.util.cat.vols  = vols;
matlabbatch{1}.spm.util.cat.name  = pth_TPM;
matlabbatch{1}.spm.util.cat.dtype = 16;
spm_jobman('run',matlabbatch);

delete(fullfile(dir_write,[nam '.mat']));

for k=1:K
    delete(vols{k});
end
%==========================================================================

%==========================================================================
function [munum,muden] = smooth_template(munum,muden,d,fwhm)
if nargin<4, fwhm = 0.5; end

K     = size(munum,1);
munum = reshape(munum',[d K]);
muden = reshape(muden',[d K]);
for k=1:K
    img = munum(:,:,:,k);            
    img = smooth_img(img,fwhm);
    munum(:,:,:,k) = img;

    img = muden(:,:,:,k);            
    img = smooth_img(img,fwhm);
    muden(:,:,:,k) = img;
end
munum = reshape(munum,[prod(d) K])';
muden = reshape(muden,[prod(d) K])';
%==========================================================================

%==========================================================================
function pth_Twarp = allocate_deformation(obj,m,s,dir_Twarp)
d0 = obj.image(1).dim;
vx = sqrt(sum(obj.image(1).mat(1:3,1:3).^2));
sk = max([1 1 1],round(obj.samp*[1 1 1]./vx));
x0 = ndgrid(1:sk(1):d0(1),1:sk(2):d0(2),1);
z0 = 1:sk(3):d0(3);
d  = [size(x0) length(z0)];

Twarp     = zeros([d 3],'single');    
pth_Twarp = fullfile(dir_Twarp,['Twarp_m' num2str(m) '_s' num2str(s) '.mat']);     
save(pth_Twarp,'Twarp');
%==========================================================================