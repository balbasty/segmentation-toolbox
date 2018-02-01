function [pth_obj,niter,fig,rand_subjs,pth_H,pth_gr] = init_obj(V,im,pth_template,uniform,pars)
print_ll     = pars.print_ll;
ml           = pars.ml;
do_bf        = pars.do_bf;
do_def       = pars.do_def;
dir_output   = pars.dir_output;
show_fig     = pars.show_fig;
write        = pars.write;

M        = numel(V);
modality = im{1}{3}; 
samp     = im{1}{5};

%--------------------------------------------------------------------------
dir_subjects = fullfile(dir_output,'subjects');
if exist(dir_subjects,'dir'), rmdir(dir_subjects,'s'); end; mkdir(dir_subjects);   
        
%--------------------------------------------------------------------------
tot_S = 0;
for m=1:M
    S = numel(V{m});    
    for s=1:S, tot_S = tot_S + 1; end
end

%--------------------------------------------------------------------------
if tot_S>1 
    do_push_resp = true;
    niter        = 30;    
else
    do_push_resp = false;
    niter        = 1;    
end

%--------------------------------------------------------------------------
fig = cell(1,4);
if show_fig
    for i=1:numel(fig), fig{i} = figure(i); clf(figure(i)); end
end
if niter>1
    fig{5} = figure(5); clf(figure(5));
    fig{6} = figure(6); clf(figure(6));
    fig{7} = figure(7); clf(figure(7));
end
drawnow

%--------------------------------------------------------------------------
if M==1, S0 = 8;
else,    S0 = floor(8/M);
end

%--------------------------------------------------------------------------
pth_obj    = {};
rand_subjs = cell(1,M);
for m=1:M
    S          = numel(V{m});    
    rand_subjs{m} = randperm(S,min(S0,S));   
    for s=1:S
        obj = struct;
             
        %------------------------------------------------------------------
        dir_s     = fullfile(dir_subjects,['s-' num2str(s) '_m-' num2str(m)]);
        if exist(dir_s,'dir'), rmdir(dir_s,'s'); end; mkdir(dir_s);    
        obj.dir_s = dir_s;
        
        %------------------------------------------------------------------
        N         = numel(V{m}{s});
        obj.image = V{m}{s};
        [~,nam]   = fileparts(obj.image(1).fname);
        
        %------------------------------------------------------------------
        obj.s        = s;
        obj.iter     = 1;
        obj.tot_S    = tot_S;
        obj.bb       = NaN(2,3);
%         obj.bb       = [-90 -126 -72; 90 90 108];
        obj.vox      =  NaN;
        obj.cleanup  =  0;
        obj.mrf      =  1;
        obj.affreg   = 'mni';
        obj.reg      = [0 0.001 0.5 0.05 0.2]*0.1;
        obj.fwhm     = 1;
        obj.samp     = samp;
        obj.Affine   = eye(4);
%         if strcmp(modality,'CT')
%             obj.biasreg  = 1e-4*(1/5)*ones(1,N);
%             obj.biasfwhm = 30*ones(1,N);
%         elseif strcmp(modality,'MRI')
            obj.biasreg  = 1e-3*(1/5)*ones(1,N);
            obj.biasfwhm = 60*ones(1,N);    
%         end
        obj.modality     = modality;
        obj.ml           = ml;
        obj.do_bf        = do_bf;
        obj.do_def       = do_def;
        obj.print_ll     = print_ll;
        obj.write        = write;
        obj.do_push_resp = do_push_resp;
        obj.bf_dc        = zeros(1,N);
        obj.avg_bf_dc    = zeros(1,N);
        obj.ll_template  = 0;
        obj.pth_template = pth_template;
        
        %------------------------------------------------------------------
        if tot_S>1
            obj.nsubit = 12;
            obj.niter  = 1;
            obj.nitgmm = 20;      
            obj.tol1   = 1e-5;
            obj.do_def = false;
        else
            obj.nsubit = 12;
            obj.niter  = 30;
            obj.nitgmm = 20;
            obj.tol1   = 1e-4;
        end
        
        %------------------------------------------------------------------
        obj.uniform = uniform; 
        V1          = spm_vol(pth_template);        
        obj.lkp     = 1:numel(V1);   
        % obj.lkp      = [1 1 2 2 3 3 4 4 5 5 5 6 6];

        %------------------------------------------------------------------
        K      = max(obj.lkp);
        pr.m   = zeros(N,K);
        if strcmp(modality,'CT')
            pr.m = pr.m + 1000;
        end
        pr.b   = ones(1,K);
        pr.n   = (N - .999)*ones(1,K);
        pr.W   = repmat(eye(N,N),[1 1 K]);
        obj.pr = pr;

        %------------------------------------------------------------------        
        dir_resp = fullfile(dir_s,'resp');
        if exist(dir_resp,'dir'), rmdir(dir_resp,'s'); end; mkdir(dir_resp);   
        pth_resp = cell(1,K);
        for k=1:K
            fname       = fullfile(dir_resp,['resp-', num2str(k), '-' nam, '.nii']);
            pth_resp{k} = fname;  
        end
        obj.pth_resp = pth_resp;
  
        %------------------------------------------------------------------               
        pth_H      = fullfile(fileparts(pth_template),'H.nii');  
        obj.pth_H  = pth_H;
        pth_gr     = fullfile(fileparts(pth_template),'gr.nii');  
        obj.pth_gr = pth_gr;
        
        %------------------------------------------------------------------
        dir_vel     = fullfile(dir_s,'vel');
        if exist(dir_vel,'dir'), rmdir(dir_vel,'s'); end; mkdir(dir_vel);            
        obj.pth_vel = fullfile(dir_vel,['vel-' nam '.nii']);   
        
        %------------------------------------------------------------------
        dir_write     = fullfile(dir_s,'write');
        if exist(dir_write,'dir'), rmdir(dir_write,'s'); end; mkdir(dir_write);            
        obj.dir_write = dir_write;            
        
        %------------------------------------------------------------------
        dir_obj       = fullfile(dir_s,'obj');
        if exist(dir_obj,'dir'), rmdir(dir_obj,'s'); end; mkdir(dir_obj);                   
        pth_obj{m}{s} = fullfile(dir_obj,['obj-' nam '.mat']);        
        save(pth_obj{m}{s},'-struct','obj')
    end
end
%==========================================================================