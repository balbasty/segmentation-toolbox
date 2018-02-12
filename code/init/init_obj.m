function [obj,niter,fig,rand_subjs] = init_obj(V,im,pth_template,uniform,do_show_seg,do_show_results,pars)
M = numel(V);

%--------------------------------------------------------------------------
dir_subjects = fullfile(pars.dir_output,'subjects');
if exist(dir_subjects,'dir'), rmdir(dir_subjects,'s'); end; mkdir(dir_subjects);   

%--------------------------------------------------------------------------
tot_subj = 0;
for m=1:M
    S = numel(V{m});    
    for s=1:S, tot_subj = tot_subj + 1; end
end

%--------------------------------------------------------------------------
if tot_subj>1     
    niter = pars.niter;    
else    
    niter = 1;    
end

%--------------------------------------------------------------------------
fig = cell(1,4);
if do_show_seg
    for i=1:numel(fig), fig{i} = figure(i); clf(figure(i)); end
    drawnow
end
if niter>1
    if do_show_results>0, fig{5} = figure(5); clf(figure(5)); end
    if do_show_results>1, fig{6} = figure(6); clf(figure(6)); end
    if do_show_results>2, fig{7} = figure(7); clf(figure(7)); end
    drawnow
end

%--------------------------------------------------------------------------
if M==1, S0 = 8;
else,    S0 = floor(8/M);
end

%--------------------------------------------------------------------------
V1   = spm_vol(pth_template); 
d1   = [V1(1).dim numel(V1)];
d_gr = [d1(1:3),d1(4) - 1];
d_H  = [d1(1:3) round(((d1(4) - 1)*d1(4))/2)];

%--------------------------------------------------------------------------
obj        = cell(1,M);
rand_subjs = cell(1,M);
cnt        = 1;
for m=1:M
    S             = numel(V{m});       
    rand_subjs{m} = randperm(S,min(S0,S));   
    obj{m}        = cell(1,S);
    for s=1:S
        modality   = im{m}{3}; 
        samp       = im{m}{5};
        init_clust = im{m}{6};
        
        obj_s = struct;
             
        %------------------------------------------------------------------
        dir_s     = fullfile(dir_subjects,['s-' num2str(s) '_m-' num2str(m)]);
        mkdir(dir_s); 
        obj_s.dir_s = dir_s;
        
        %------------------------------------------------------------------
        N         = numel(V{m}{s});
        obj_s.image = V{m}{s};
        [~,nam]   = fileparts(obj_s.image(1).fname);
        
        %------------------------------------------------------------------
        obj_s.s            = cnt;
        obj_s.iter         = 0;
        obj_s.tot_S        = tot_subj;
        obj_s.bb           = NaN(2,3);
%         obj.bb       = [-90 -126 -72; 90 90 108];
        obj_s.vox          =  NaN;
        obj_s.cleanup      =  0;
        obj_s.mrf          =  1;
        obj_s.affreg       = 'mni';
        obj_s.reg          = [0 0.001 0.5 0.05 0.2]*0.1;
        obj_s.reg0         = obj_s.reg;
        obj_s.fwhm         = 1;
        obj_s.samp         = samp;
        obj_s.Affine       = eye(4);
        obj_s.biasreg      = 1e-3*(1/5)*ones(1,N);
        obj_s.biasfwhm     = 60*ones(1,N);  
        obj_s.modality     = modality;
        obj_s.ml           = pars.ml;
        obj_s.do_bf        = pars.do_bf;
        if strcmp(modality,'CT')
            obj_s.do_bf    = false;
        end
        obj_s.do_def       = pars.do_def;
        obj_s.do_wp        = pars.do_wp;
        obj_s.print_seg    = pars.print_seg;
        obj_s.print_ll     = pars.print_ll;
        obj_s.do_write_res = pars.do_write_res;
        obj_s.do_push_resp = false;
        obj_s.bf_dc        = zeros(1,N);
        obj_s.avg_bf_dc    = zeros(1,N);
        obj_s.ll_template  = 0;
        obj_s.pth_template = pth_template;        
        obj_s.status       = 0;
        obj_s.true         = true;
        obj_s.nsubit       = 8;
        obj_s.nitgmm       = 20;
        obj_s.niter        = 30;            
        obj_s.tol1         = 1e-4;
        obj_s.gmm          = struct;
        obj_s.kmeans_dist  = pars.kmeans_dist;
        obj_s.init_clust   = init_clust;
        obj_s.missing_data = pars.missing_data;

        %------------------------------------------------------------------
        obj_s.uniform = uniform;            
        Kb            = numel(V1);
        if isfield(pars,'lkp')
            obj_s.lkp = pars.lkp;
            if Kb~=max(obj_s.lkp)
               error('Kb~=max(obj.lkp)') 
            end
        else
            obj_s.lkp = 1:Kb;   
        end
        
        K1          = 1:Kb;
        K_rem       = im{m}{4};
        msk         = ismember(K1,K_rem);
        K_keep      = K1(~msk);
        obj_s.K_lab = {K_keep,K_rem};
        
        %------------------------------------------------------------------   
%         if ~obj_s.ml
%             pr.m   = zeros(N,Kb);
%             pr.b   = ones(1,Kb);
%             pr.n   = (N - .999)*ones(1,Kb);
%             pr.W   = repmat(eye(N,N),[1 1 Kb]);
% 
%             obj_s.gmm.pr = pr;
%         end
        
        %------------------------------------------------------------------        
        dir_resp = fullfile(dir_s,'resp'); 
        mkdir(dir_resp); 
        
        pth_resp = cell(1,Kb);
        for k=1:Kb
            fname       = fullfile(dir_resp,['resp-', num2str(k), '-' nam, '.nii']);
            pth_resp{k} = fname;  
        end
        obj_s.pth_resp = pth_resp;
  
        %------------------------------------------------------------------      
        dir_der = fullfile(dir_s,'der');      
        mkdir(dir_der);
        
        pth_gr = cell(1,d_gr(end));
        for k=1:d_gr(end)
            fname     = fullfile(dir_der,['gr-', num2str(k), '-' nam, '.nii']);
            pth_gr{k} = fname;  
        end
        obj_s.pth_gr = pth_gr;
        
        pth_H = cell(1,d_H(end));
        for k=1:d_H(end)
            fname    = fullfile(dir_der,['H-', num2str(k), '-' nam, '.nii']);
            pth_H{k} = fname;  
        end
        obj_s.pth_H = pth_H;
        
        %------------------------------------------------------------------
        dir_vel = fullfile(dir_s,'vel');   
        mkdir(dir_vel);
        
        obj_s.pth_vel = fullfile(dir_vel,['vel-' nam '.nii']);   
        
        %------------------------------------------------------------------
        dir_write = fullfile(dir_s,'write');   
        mkdir(dir_write);
        
        obj_s.dir_write = dir_write;            
        
        %------------------------------------------------------------------
        obj{m}{s} = obj_s;
        
        cnt = cnt + 1;
    end
end
%==========================================================================
