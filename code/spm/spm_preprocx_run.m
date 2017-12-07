function spm_preprocx_run(obj,im,K)

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

[obj,K] = init_TPM(obj,V,K);

%==========================================================================
% Initialise some debugging output
%==========================================================================

[fig_seg,fig_TPM,fig_L,obj.verbose] = create_fig(obj,false);

%==========================================================================
% Initialise algorithm i/o 
%==========================================================================

nitermain    = obj.nitermain;
tolmain      = obj.tolmain;
dotpm        = obj.dotpm;
vb           = obj.vb;
use_mog      = ~isempty(obj.lkp);
dopr         = obj.dopr && use_mog && vb;
pth_logTPM   = obj.pth_logTPM;
tiny         = obj.tiny;
deg          = obj.deg;
fwhm_TPM     = obj.fwhm_TPM;
dir_data     = obj.dir_data;
dir_res      = obj.dir_res;
run_on_holly = obj.run_on_holly;
num_workers  = obj.num_workers;

pth_obj = get_pth_obj(obj,V,K,im,labels,run_on_holly);
clear V obj

%==========================================================================
% Prepare for running the algorithm on the Holly cluster
% OBS! For this to work, this code needs to be located on the Holly server
%==========================================================================
   
if run_on_holly
    holly = init_holly(pth_obj,dir_data,obj.holly_jnam);
else
    holly = [];
end

%==========================================================================
% Run algorithm
%==========================================================================

print_algorithm_progress('started');

L = -Inf; % Lower bound of complete model
for iter=1:nitermain
    fprintf('iter=%d======================\n',iter);  
                   
    % Update subject specific parameters (clusters, bias, affine, deformations, template derivatives)
    [L,munum,muden,Nm,holly] = update_subjects(pth_obj,L,num_workers,run_on_holly,holly,fig_seg,pth_logTPM); 
     
    % Plot objective value
    plot_objval(L,fig_L); 
        
    if dotpm
        % Update template 
        update_global_TPM(munum,muden,pth_logTPM,tiny,fwhm_TPM,pth_obj);                
              
        show_TPM(fig_TPM,pth_logTPM,tiny,deg);
    end  
    
    if dopr
        % Update intensity prior
        update_global_prior(pth_obj,im);
    end
        
    % Check convergence
    fprintf('L=%d\n',L(end));    
    if ~(abs(L(end)-L(end - 1))>2*tolmain*Nm) && iter>9
        fprintf('==============================================\n')                        
        fprintf('Algorithm converged in %d iterations.\n',iter)                        
        fprintf('==============================================\n')                        
          
        % Convert log template to probability template and save to disk
        logTPM2TPM(pth_logTPM,dir_res);
    
        break;
    end
end

print_algorithm_progress('finished');
%==========================================================================

%==========================================================================
function [L,munum,muden,Nm,holly] = update_subjects(pth_obj,L,num_workers,run_on_holly,holly,fig,pth_logTPM)
V      = spm_vol(pth_logTPM);
dm_tpm = V(1).dim;
K      = numel(V);
clear V

if run_on_holly
    % Run subject specific jobs on the FIL cluster (Holly)
    %----------------------------------------------------------------------            
    [ll,munum,muden,Nm,status,holly] = parfor_holly(pth_obj,holly,dm_tpm,K);
else 
    % Run subject specific jobs using MATLAB parfor
    %----------------------------------------------------------------------
    [ll,munum,muden,Nm,status] = parfor_matlab(pth_obj,num_workers,fig,dm_tpm,K);
end
L = [L,ll];

if status
    % Display how many segmentations that have errored
    fprintf(2,'==============================================\n');
    fprintf(2,'%d job(s) with status~=0\n',status);
    fprintf(2,'==============================================\n\n');
end
%==========================================================================

%==========================================================================
function [ll,munum,muden,Nm,tot_status] =  parfor_matlab(pth_obj,num_workers,fig,dm_tpm,K)
% Estimate on all subjects
M = numel(pth_obj);
manage_parpool(num_workers);
for m=1:M                    
    pth_obj_m = pth_obj{m};
    S         = numel(pth_obj_m);    
    parfor (s=1:S,num_workers)                                  
%     for s=1:S
        obj = load(pth_obj_m{s});            
        
        % Run segmentation routine         
        obj = segment(obj,fig);                                       

        save_in_parfor(pth_obj_m{s},obj,'-struct');
    end     
end  

% Read results from estimations
tic;
munum = zeros([dm_tpm K],'single'); muden = munum; ll = 0; Nm = 0; tot_status = 0;
for m=1:M                  
    S = numel(pth_obj{m});
    for s=1:S 
        fprintf('.');
        obj = load(pth_obj{m}{s}); 
        
        tot_status = tot_status + obj.status;
        if obj.status==0
            [~,ix_x,ix_y,ix_z] = bb_info(obj.bb);

            Nii    = nifti(obj.pth_munum);
            munum1 = Nii.dat(:,:,:,:); 
            munum(ix_x,ix_y,ix_z,:) = munum(ix_x,ix_y,ix_z,:) + munum1;

            Nii    = nifti(obj.pth_muden);
            muden1 = Nii.dat(:,:,:,:); 
            muden(ix_x,ix_y,ix_z,:) = muden(ix_x,ix_y,ix_z,:) + muden1; 

            ll = ll + obj.ll;
            Nm = Nm + obj.nm;
        end
    end
end
fprintf('\n');
fprintf('Elapsed time (parfor_matlab): %d s\n',round(toc))      
%==========================================================================

%==========================================================================
function update_global_prior(pth_obj,im)
tic
M = numel(im);

% Map modalities
img_mod = cell(1,M);
for m=1:M
    img_mod{m} = im{m}{3};
end
[~,M1,img_mod] = unique(img_mod);
M1             = numel(M1);
% img_mod = 1:M;
% M1      = M;

% Load posteriors and priors from all subjects and store in struct
obj = cell(1,M1);
for m=1:M
    m1     = img_mod(m);
    S      = numel(pth_obj{m});    
    cnt    = numel(obj{m1});
    for s=1:S
        tmp = load(pth_obj{m}{s},'-mat','po','pr','status');
        if tmp.status==0
            cnt             = cnt + 1;
            obj{m1}{cnt}.po = tmp.po;
            obj{m1}{cnt}.pr = tmp.pr;            
        end
    end
end

% Update prior based on posteriors and previous priors, then save new prior
pr = cell(1,M1);
for m=1:M1
    if ~isempty(obj{m})
        pr{m} = update_intensity_prior(obj{m});        
    end
end
clear obj

for m=1:M  
    m1        = img_mod(m);
    S         = numel(pth_obj{m});
    for s=1:S    
        obj1    = matfile(pth_obj{m}{s},'Writable',true);
        obj1.pr = pr{m1};   
    end              
end
toc
%==========================================================================        

%==========================================================================
function update_global_TPM(munum,muden,pth_logTPM,tiny,fwhm_TPM,pth_obj,softmax_TPM,dir_res)
if nargin<7, softmax_TPM = false; end

% Smooth template
[munum,muden] = smooth_template(munum,muden,fwhm_TPM);

% Update template
logmu = log(munum./muden + tiny);
clear munum muden

% Save updated template
Nii = nifti(pth_logTPM);
Nii.dat(:,:,:,:) = logmu;        
clear Nii logmu

if 0
    % Crop template according to subject bounding boxes
    M = numel(pth_obj);
    S = 0;
    for m=1:M, S = S + numel(pth_obj{m}); end

    bb  = [];
    cnt = 1;
    for m=1:M
        S = numel(pth_obj{m});    
        for s=1:S
            obj = load(pth_obj{m}{s},'-mat','bb','status');
            if obj.status==0   
                bb1  = obj.bb;
                is3D = numel(bb1)==6;
                if is3D
                    % 3D
                    bb1 = [bb1(3) bb1(4);bb1(1) bb1(2);bb1(5) bb1(6)];
                else
                    % 2D
                    bb1 = [bb1(3) bb1(4);bb1(1) bb1(2);1 1];
                end
                bb(:,:,cnt) = bb1;
                cnt         = cnt + 1;
            end
        end
    end

    if is3D
        mn_bb = min(bb,[],3);
        mx_bb = max(bb,[],3);
        nbb   = [mn_bb(:,1) mx_bb(:,2)];
        nbb   = floor(nbb);
        nbb(nbb==0) = 1;

        V  = spm_vol(pth_logTPM);
        od = V(1).dim;

        for k=1:numel(V)
            subvol(V(k),nbb','tmp');        
        end

        delete(pth_logTPM);
        [pth,nam,ext] = fileparts(V(1).fname);
        fname         = fullfile(pth,['tmp' nam ext]);
        movefile(fname,pth_logTPM);

        V  = spm_vol(pth_logTPM);
        nd = V(1).dim;

        fprintf('dim(TPM)=[%s], dim(nTPM)=[%s]\n\n',sprintf('%d ',od),sprintf('%d ',nd));
    end
end

if softmax_TPM
    % Convert log template to probability template and save to disk
    logTPM2TPM(pth_logTPM,dir_res);
end
%==========================================================================

%==========================================================================
function show_TPM(fig_TPM,pth_logTPM,tiny,deg)
if ~isempty(fig_TPM) 
    logtpm  = spm_load_logpriors8(pth_logTPM,tiny,deg);

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
        title(['mu, k=' num2str(i)]);
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
    title('Lower bound')
end
%==========================================================================

%==========================================================================
function pth_obj = get_pth_obj(obj0,V,K,im,labels,run_on_holly)
dir_data = obj0.dir_data;

dir_obj  = fullfile(dir_data,'obj');
if exist(dir_obj,'dir'), rmdir(dir_obj,'s'); end; mkdir(dir_obj);

dir_def  = fullfile(dir_data,'def');
if exist(dir_def,'dir'), rmdir(dir_def,'s'); end; mkdir(dir_def);

dir_munum  = fullfile(dir_data,'munum');
if exist(dir_munum,'dir'), rmdir(dir_munum,'s'); end; mkdir(dir_munum);

dir_muden  = fullfile(dir_data,'muden');
if exist(dir_muden,'dir'), rmdir(dir_muden,'s'); end; mkdir(dir_muden);

M       = numel(V);
pth_obj = cell(1,M);
for m=1:M       
    S          = numel(V{m});
    pth_obj{m} = cell(1,S);
    for s=1:S  
        obj = struct;
        
        obj.s      = s;
        obj.m      = m;
        obj.status = false;
        
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
        obj.dodef    = obj0.dodef; 
        obj.dodef0   = obj0.dodef; 
        obj.dotpm    = obj0.dotpm;  
        obj.doaff    = obj0.doaff;
        obj.dowp     = obj0.dowp;
        obj.dowp0    = obj0.dowp;
                   
        obj.tolseg    = obj0.tolseg;
        obj.niter     = obj0.niter;
        obj.niter1    = obj0.niter1;
        obj.nsubitmog = obj0.nsubitmog;
        obj.nsubitbf  = obj0.nsubitbf;
        obj.nitdef    = obj0.nitdef;

        obj.descrip = im{m}{3};
        obj.healthy = im{m}{4};         
        
        obj.def_done = false;
        
        obj.ll = 0;
        obj.nm = 0;
              
        obj.bb = [];
        
        obj.pth_logTPM = obj0.pth_logTPM;
        if run_on_holly
            fname          = obj.pth_logTPM;
            nfname         = ['/' fname(7:end)];
            obj.pth_logTPM = nfname;            
        end
        
        obj.tiny = obj0.tiny;        
        obj.deg  = obj0.deg;        
        obj.iter = 0;
               
        % Path to deformations
        pth_def     = fullfile(dir_def,['def-m' num2str(m) '-s' num2str(s) '.nii']);                       
        obj.pth_def = pth_def;        
        if run_on_holly        
            fname       = pth_def;
            nfname      = ['/' fname(7:end)];
            obj.pth_def = nfname;                        
        end    
        
        % Path to template updates
        pth_munum     = fullfile(dir_munum,['munum-m' num2str(m) '-s' num2str(s) '.nii']);  
        obj.pth_munum = pth_munum;                
        if run_on_holly        
            obj.pth_munum_l = pth_munum;
            fname           = pth_munum;
            nfname          = ['/' fname(7:end)];
            obj.pth_munum   = nfname;                        
        end   
        
        pth_muden     = fullfile(dir_muden,['muden-m' num2str(m) '-s' num2str(s) '.nii']); 
        obj.pth_muden = pth_muden;                
        if run_on_holly        
            obj.pth_muden_l = pth_muden;
            fname           = pth_muden;
            nfname          = ['/' fname(7:end)];
            obj.pth_muden   = nfname;                        
        end   
        
        % Allocate labels
%         obj.labels  = labels{m}{s}; 
        
        % Store path to obj        
        pth_obj{m}{s} = fullfile(dir_obj,['obj-m' num2str(m) '-s' num2str(s) '.mat']);
        
        save(pth_obj{m}{s},'-struct','obj')
    end
end
%==========================================================================

%==========================================================================
function [fig_seg,fig_TPM,fig_L,verbose] = create_fig(obj,distf)
if nargin<2, distf = true; end

verbose      = obj.verbose;
num_workers  = obj.num_workers;
figix        = obj.figix;
dotpm        = obj.dotpm;
run_on_holly = obj.run_on_holly;

fig_seg = cell(4,1);

if verbose==2 && ~num_workers && ~run_on_holly
    spm_figure('Create','Interactive');
    figure(figix)

    for i=1:size(fig_seg,1)              
        fig_seg{i} = figure(figix + i);
        clf(fig_seg{i})
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

if distf
    try
        distFig; 
    catch
        warning('distFig not available')
    end
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
function [munum,muden] = smooth_template(munum,muden,fwhm)
if nargin<3, fwhm = 0.5; end

K = size(munum,4);
for k=1:K
    img = munum(:,:,:,k);            
    img = smooth_img(img,fwhm);
    munum(:,:,:,k) = img;

    img = muden(:,:,:,k);            
    img = smooth_img(img,fwhm);
    muden(:,:,:,k) = img;
end
%==========================================================================

%==========================================================================
function print_algorithm_progress(status)
fprintf('==============================================\n')  
fprintf(['Algorithm ' status ' ('])
fprintf(datestr(now,'mmmm dd, yyyy HH:MM:SS'))
fprintf(')\n')
fprintf('==============================================\n\n')   
%==========================================================================