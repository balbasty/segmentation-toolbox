function spm_preprocx_run(obj,im,K)

%==========================================================================
% Prepare or disable parallel processing
%==========================================================================

manage_parpool(obj.num_workers);

%==========================================================================
% Make some directories
%==========================================================================

if obj.run_on_holly
    [~,obj.dir_data,~] = read_directory_details('directory_details.txt',obj.holly.jnam);    
end

if ~exist(obj.dir_data,'dir'), mkdir(obj.dir_data); end
obj.dir_res = fullfile(obj.dir_data,obj.dir_res); 
if exist(obj.dir_res,'dir'), rmdir(obj.dir_res,'s'); end; mkdir(obj.dir_res);

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

nitermain    = obj.nitermain;
tolmain      = obj.tolmain;
dopr         = obj.dopr;
dotpm        = obj.dotpm;
vb           = obj.vb;
use_mog      = ~isempty(obj.lkp);
pth_logTPM   = obj.pth_logTPM;
tiny         = obj.tiny;
deg          = obj.deg;
fwhm_TPM     = obj.fwhm_TPM;
mrf          = obj.mrf;
dir_data     = obj.dir_data;
dir_res      = obj.dir_res;
num_workers  = obj.num_workers;
run_on_holly = obj.run_on_holly;
holly        = obj.holly; 

pth_obj = get_pth_obj(obj,V,K,im,labels,num_workers,run_on_holly);
clear V obj

%==========================================================================
% Prepare for running the algorithm on the Holly cluster
% OBS! For this to work, this code needs to be located on the Holly server
%==========================================================================
   
holly = init_holly(pth_obj,dir_data,run_on_holly,holly);

%==========================================================================
% Run algorithm
%==========================================================================

print_algorithm_started;

L = -Inf; % Lower bound of complete model
for iter=1:nitermain
    fprintf('iter=%d======================\n',iter);  
                   
    % Update subject specific parameters (clusters, bias, affine, deformations, template derivatives)
    [L,munum,muden,Nm] = update_subjects(pth_obj,L,num_workers,run_on_holly,holly); 
     
    plot_objval(L,fig_L); 
        
    if dotpm
        % Update the template 
        update_global_TPM(munum,muden,pth_logTPM,tiny,deg,fwhm_TPM,mrf,dir_res);                
              
        show_TPM(fig_TPM,pth_logTPM,tiny,deg);
    end  
    
    if dopr && use_mog && vb && iter>=3                                      
        % Update subject specific parameters (clusters, bias, deformations)
        L = update_subjects(pth_obj,L,num_workers,run_on_holly,holly); 
              
        % Update intensity prior
        update_global_prior(pth_obj,dir_res,num_workers);
        
        plot_objval(L,fig_L);
    end
        
    % Check convergence----------------------------------------------------
    fprintf('L=%d\n',L(end));
    
    if ~(abs(L(end)-L(end - 1))>2*tolmain*Nm)        
        fprintf('==============================================\n')                        
        fprintf('Algorithm converged in %d iterations.\n',iter)                        
        fprintf('==============================================\n')                        
          
        break;
    end
end
%==========================================================================

%==========================================================================
function [L,munum,muden,Nm] = update_subjects(pth_obj,L,num_workers,run_on_holly,holly)

% Set flag deciding if template calculations should be performed in
% spm_preproc or not
%--------------------------------------------------------------------------
M     = numel(pth_obj);
dotpm = nargout>=2;
for m=1:M
    S         = numel(pth_obj{m});
    pth_obj_m = pth_obj{m};  
    parfor (s=1:S,num_workers)
        obj          = matfile(pth_obj_m{s},'Writable',true);
        obj.diff_TPM = dotpm;  
    end
end


if run_on_holly
    % Run subject specific jobs on the FIL cluster (Holly)
    %----------------------------------------------------------------------            
    [ll,munum,muden,Nm] = parfor_holly(pth_obj,holly,num_workers);
else 
    % Run subject specific jobs using MATLAB parfor
    %----------------------------------------------------------------------
    [ll,munum,muden,Nm] =  parfor_matlab(pth_obj,num_workers);
end
L = [L,ll];

% Display how many segmentations have errored
%--------------------------------------------------------------------------
tot_status = 0;
for m=1:M
    S         = numel(pth_obj{m});
    pth_obj_m = pth_obj{m};  
    parfor (s=1:S,num_workers)
        tmp        = load(pth_obj_m{s},'-mat','status');       
        tot_status = tot_status + tmp.status;
    end
end
if tot_status
    fprintf('%d number of job(s) failed.\n',tot_status);
end
%==========================================================================

%==========================================================================
function [ll,munum,muden,Nm] =  parfor_matlab(pth_obj,num_workers)
M     = numel(pth_obj);
munum = 0; muden = 0; ll = 0; Nm = 0;
for m=1:M                    
    pth_obj_m = pth_obj{m};
    S         = numel(pth_obj_m);
%         for s=1:S                                       
    parfor (s=1:S,num_workers)    
        % Run segmentation routine                               
        obj = load(pth_obj_m{s});            
        obj = segment(obj);                                       

        if obj.status==0
            munum = munum + double(obj.munum);
            muden = muden + double(obj.muden);            
            ll    = ll + obj.ll;
            Nm    = Nm + obj.nm;
        end

        save_in_parfor(pth_obj_m{s},obj,'-struct');
    end        
end  
%==========================================================================

%==========================================================================
function update_global_prior(pth_obj,dir_res,num_workers)
M   = numel(pth_obj);
obj = cell(1,M);

% Load posteriors and priors from all subjects and store in struct
for m=1:M
    S         = numel(pth_obj{m});
    pth_obj_m = pth_obj{m};
    obj{m}    = cell(1,S);
    obj_m     = obj{m};
%     for s=1:S
    parfor (s=1:S,num_workers)
        tmp         = load(pth_obj_m{s},'-mat','po','pr');        
        obj_m{s}.po = tmp.po;
        obj_m{s}.pr = tmp.pr;
    end
    obj{m} = obj_m;
end

% Update prior based on posteriors and previous priors, then save new prior
for m=1:M
    pr  = update_intensity_prior(obj{m});
    save(fullfile(dir_res,['pr_m' num2str(m) '.mat']),'pr'); 
    
    S         = numel(obj{m});
    pth_obj_m = pth_obj{m};
    parfor (s=1:S,num_workers)
        obj1    = matfile(pth_obj_m{s},'Writable',true);
        obj1.pr = pr;   
    end          
end
%==========================================================================        

%==========================================================================
function update_global_TPM(munum,muden,pth_logTPM,tiny,deg,fwhm_TPM,mrf,dir_res)

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
    logtpm = spm_load_logpriors8(pth_logTPM,tiny,deg,0);

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
function show_TPM(fig_TPM,pth_logTPM,tiny,deg)
if ~isempty(fig_TPM) 
    uniform = 0;
    logtpm  = spm_load_logpriors8(pth_logTPM,tiny,deg,uniform);

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
function plot_objval(L,fig_L)
if ~isempty(fig_L) && numel(L)>=3
    set(0,'CurrentFigure',fig_L);                
    plot(0:numel(L(3:end)) - 1,L(3:end),'b-','LineWidth',1);   hold on            
    plot(0:numel(L(3:end)) - 1,L(3:end),'b.','markersize',10); hold off  
end
%==========================================================================

%==========================================================================
function pth_obj = get_pth_obj(obj0,V,K,im,labels,num_workers,run_on_holly)
dir_data = obj0.dir_data;
dir_obj  = fullfile(dir_data,'obj');
if exist(dir_obj,'dir'), rmdir(dir_obj,'s'); end; mkdir(dir_obj);

M       = numel(V);
pth_obj = cell(1,M);
for m=1:M    
    S          = numel(V{m});
    pth_obj{m} = cell(1,S);  
    pth_obj_m  = pth_obj{m};
%     for s=1:S
    parfor (s=1:S,num_workers)
        obj = struct;
        
        obj.s      = s;
        obj.m      = m;
        obj.status = 0;
        
        N         = numel(V{m}{s});
        obj.image = V{m}{s};
        if run_on_holly
            for n=1:N
                fname              = obj.image(n).fname;
                nfname             = ['/' fname(7:end)];
                obj.image(n).fname = nfname;            
            end
        end        
        
        obj.biasfwhm = obj0.biasfwhm*ones(1,N);
        obj.biasreg  = obj0.biasreg*ones(1,N);       

        obj.use_tpm = obj0.use_tpm;
        
        obj.vb         = obj0.vb;        
        obj.wp_reg     = obj0.wp_reg;
        obj.clear_pars = 1;
        if obj0.dotpm
            % One Gaussian per tissue for template construction
            obj.lkp = 1:K; 
        else
            if isempty(obj0.lkp)
                obj.lkp = 0;
            elseif isscalar(obj0.lkp)
                obj.lkp = repelem(1:K,obj0.lkp);
            else
                obj.lkp = obj0.lkp;
            end
        end        

        obj.Affine  = eye(4);
        obj.reg     = obj0.rparam;
        obj.samp    = obj0.samp;
        obj.fwhm    = 0;
        obj.verbose = obj0.verbose;

        obj.dobias   = obj0.dobias;
        obj.dodef0   = obj0.dodef; 
        obj.dotpm    = obj0.dotpm;  
        obj.diff_TPM = obj0.dotpm;  
        obj.doaff    = obj0.doaff;
        obj.dopr     = obj0.dopr;
        obj.dowp     = obj0.dowp;
        obj.dowp0    = obj0.dowp;
                
        obj.nitermain = obj0.nitermain;
        obj.tolmain   = obj0.tolmain;       
        obj.tolseg    = obj0.tolseg;
        obj.niter     = obj0.niter;
        obj.niter1    = obj0.niter1;
        obj.nsubitmog = obj0.nsubitmog;
        obj.nsubitbf  = obj0.nsubitbf;
        obj.nitdef    = obj0.nitdef;

        obj.descrip = im{m}{3};
        obj.healthy = im{m}{4}; 
        obj.labels  = labels{m}{s}; 
        
        obj.munum = 0;
        obj.muden = 0;
        obj.ll    = 0;
        obj.nm    = 0;
                    
        obj.fig = obj0.fig; 
              
        obj.pth_logTPM = obj0.pth_logTPM;
        if run_on_holly
            fname          = obj.pth_logTPM;
            nfname         = ['/' fname(7:end)];
            obj.pth_logTPM = nfname;            
        end
        
        obj.fwhm_TPM = obj0.fwhm_TPM;
        obj.tiny     = obj0.tiny;
        obj.mrf      = obj0.mrf;
        obj.deg      = obj0.deg;
        obj.dir_res  = obj0.dir_res;
        obj.uniform  = obj0.uniform;
        obj.iter     = 0;
        
        obj.num_workers  = obj0.num_workers;
        obj.run_on_holly = obj0.run_on_holly;
        
        pth_obj_m{s} = fullfile(dir_obj,['obj-m' num2str(m) '-s' num2str(s) '.mat']);
        save_in_parfor(pth_obj_m{s},obj,'-struct');
    end
    pth_obj{m} = pth_obj_m;
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
function manage_parpool(num_workers)
poolobj = gcp('nocreate');
if ~isempty(poolobj)    
    if num_workers==0
        delete(poolobj);
    elseif poolobj.NumWorkers~=num_workers
        delete(poolobj);
        parpool('local',num_workers);
    end
end
%==========================================================================

%==========================================================================
function print_algorithm_started
fprintf('==============================================\n')  
fprintf('Algorithm started (')
fprintf(datestr(now,'mmmm dd, yyyy HH:MM:SS'))
fprintf(')\n')
fprintf('==============================================\n\n')   
%==========================================================================