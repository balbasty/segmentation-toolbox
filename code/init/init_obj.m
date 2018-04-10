function [obj,pars] = init_obj(pars)
M = numel(pars.dat);

%--------------------------------------------------------------------------
dir_subjects = fullfile(pars.dir_output,'subjects');
if exist(dir_subjects,'dir'), rmdir(dir_subjects,'s'); end; mkdir(dir_subjects);   

%--------------------------------------------------------------------------
tot_S = 0;
for m=1:M
    S = numel(pars.dat{m}.V);    
    for s=1:S
        tot_S = tot_S + 1; 
    end
end

if ~(tot_S>1)   
    pars.niter = 1; 
end

%--------------------------------------------------------------------------
fig = cell(1,4);
if pars.dat{1}.segment.verbose
    for i=1:numel(fig), fig{i} = figure(i); clf(figure(i)); end    
end
if pars.niter>1
    if pars.verbose>3, fig{8} = figure(8); clf(figure(8)); end
    if pars.verbose>2, fig{7} = figure(7); clf(figure(7)); end
    if pars.verbose>1, fig{6} = figure(6); clf(figure(6)); end
    if pars.verbose>0, fig{5} = figure(5); clf(figure(5)); end        
end
drawnow

pars.fig = fig;

%--------------------------------------------------------------------------
if M==1, S0 = 8;
else     S0 = floor(8/M);
end

%--------------------------------------------------------------------------
V1   = spm_vol(pars.pth_template); 
d1   = [V1(1).dim numel(V1)];
d_gr = [d1(1:3),d1(4) - 1];
d_H  = [d1(1:3) round(((d1(4) - 1)*d1(4))/2)];

%--------------------------------------------------------------------------
obj        = cell(1,M);
rand_subjs = cell(1,M);
cnt        = 1;
for m=1:M           
    V             = pars.dat{m}.V;
    N             = numel(V{1});
    S             = numel(V);
    rand_subjs{m} = randperm(S,min(S0,S));   
    obj{m}        = cell(1,S);
    
    pars.dat{m}.segment.biasreg  = pars.dat{m}.segment.biasreg*ones(1,N);
    pars.dat{m}.segment.biasfwhm = pars.dat{m}.segment.biasfwhm*ones(1,N);     
    
    if V{1}(1).dim(3)==1
        % Data is 2D
        %------------------------------------------------------------------
        pars.dat{m}.maff.do_maff = false;
        pars.dat{m}.segment.samp = 1;
    end
    
    pars.dat{m}.write_res.write_bf = repmat(pars.dat{m}.write_res.write_bf,N,1);  
        
    for s=1:S                
        obj1 = struct;
             
        %------------------------------------------------------------------
        dir_s      = fullfile(dir_subjects,['s-' num2str(s) '_m-' num2str(m)]);
        obj1.dir_s = dir_s;
        mkdir(dir_s);         
        
        %------------------------------------------------------------------        
        obj1.image  = V{s};
        obj1.labels = pars.dat{m}.labels{s};
        
        %------------------------------------------------------------------
        obj1.s            = cnt;
        obj1.tot_S        = tot_S;
        
        obj1.modality        = pars.dat{m}.modality;
        obj1.pth_template    = pars.pth_template;                        
        obj1.dt              = pars.dt;        
        obj1.uniform         = pars.uniform;    
        obj1.dir_template    = pars.dir_template;    
        obj1.trunc_ct        = pars.dat{m}.trunc_ct;
        obj1.fwhm            = pars.dat{m}.fwhm;
        obj1.Affine          = pars.dat{m}.Affine;
        obj1.print_subj_info = pars.dat{m}.print_subj_info;
        obj1.ll_template     = 0;
        obj1.status          = 0;
        obj1.iter            = 1;
        
        %------------------------------------------------------------------
        obj1.maff = pars.dat{m}.maff;
                
        %------------------------------------------------------------------
        obj1.push_resp = pars.dat{m}.push_resp;
        
        %------------------------------------------------------------------                 
        obj1.write_res = pars.dat{m}.write_res;        
        
        %------------------------------------------------------------------                     
        obj1.segment = pars.dat{m}.segment;
        
        obj1.segment.reg0    = pars.dat{m}.segment.reg;
        obj1.segment.do_bf0  = pars.dat{m}.segment.do_bf;
        obj1.segment.do_def0 = pars.dat{m}.segment.do_def;
        obj1.segment.do_wp0  = pars.dat{m}.segment.do_wp;
           
        obj1.segment.bf_dc     = zeros(1,N);
        obj1.segment.avg_bf_dc = zeros(1,N);
        
        % lkp
        lkp = pars.dat{m}.segment.lkp;
        
        if ~isempty(lkp.lab) && ~isempty(lkp.rem)
           error('~isempty(lkp.lab) && ~isempty(lkp.rem)') 
        end
            
        if ~isempty(pars.dat{m}.labels{s}) && ~isempty(lkp.lab)
            labels  = pars.dat{m}.labels{s};
            labels  = uint8(labels.private.dat(:,:,:));
            if min(labels(:))==0
               labels = labels + 1; 
            end
            mn = min(labels(:));
            mx = max(labels(:));
            
            if isempty(lkp.lab)
                lkp.lab = mn:mx;
                lkp.lab = [lkp.lab zeros(1,pars.K - mx)];
            end
        end
        
        if ~isempty(lkp.rem)
            msk      = ismember(lkp.keep,lkp.rem);
            lkp.keep = lkp.keep(~msk);
        end    
        
        obj1.segment.lkp = lkp;
        
        % gmm
        obj1.segment.gmm = struct;
        pth_prior        = pars.dat{m}.segment.pth_prior;
        if ~isempty(pth_prior)
            tmp                 = load(pth_prior,'-mat');
            obj1.segment.gmm.pr = tmp.pr;
        end
          
        %------------------------------------------------------------------
        [~,nam] = fileparts(obj1.image(1).fname);

        dir_resp = fullfile(dir_s,'resp'); 
        mkdir(dir_resp); 

        pth_resp = cell(1,pars.K);
        for k=1:pars.K
            fname       = fullfile(dir_resp,['resp-', num2str(k), '-' nam, '.nii']);
            pth_resp{k} = fname;  
        end
        obj1.pth_resp = pth_resp;
    
        dir_der = fullfile(dir_s,'der');      
        mkdir(dir_der);

        pth_gr = cell(1,d_gr(end));
        for k=1:d_gr(end)
            fname     = fullfile(dir_der,['gr-', num2str(k), '-' nam, '.nii']);
            pth_gr{k} = fname;  
        end
        obj1.pth_gr = pth_gr;

        pth_H = cell(1,d_H(end));
        for k=1:d_H(end)
            fname    = fullfile(dir_der,['H-', num2str(k), '-' nam, '.nii']);
            pth_H{k} = fname;  
        end
        obj1.pth_H = pth_H;

        dir_vel = fullfile(dir_s,'vel');   
        mkdir(dir_vel);

        obj1.pth_vel = fullfile(dir_vel,['vel-' nam '.nii']);   

        dir_write = fullfile(dir_s,'write');   
        mkdir(dir_write);

        obj1.dir_write = dir_write;
        
        %------------------------------------------------------------------
        obj{m}{s} = obj1;
        
        cnt = cnt + 1;
    end
end

% A random sample of subjects (used for visualising algorithm progress)
pars.rand_subjs = rand_subjs;

% If pars.dat{m}.segment.kmeans_hist, fit GMM to accumulative histogram.
% For setting initial GMM parameters.
obj = fit_gmm_to_acc_hist(obj,pars);
%==========================================================================