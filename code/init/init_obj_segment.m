function [obj,pars] = init_obj_segment(pars)
M = numel(pars.dat);

%--------------------------------------------------------------------------
dir_subjects = fullfile(pars.dir_output,'subject-data');
if exist(dir_subjects,'dir'), rmdir(dir_subjects,'s'); end; mkdir(dir_subjects);   

dir_animations      = fullfile(pars.dir_output,'animations');
if exist(dir_animations,'dir'), rmdir(dir_animations,'s'); end; mkdir(dir_animations);   
pars.dir_animations = dir_animations;

dir_local = pars.dir_local;
if ~isempty(dir_local)
    if exist(dir_local,'dir'), rmdir(dir_local,'s'); end; mkdir(dir_local);   
end

%--------------------------------------------------------------------------
fig = cell(1,11);
if pars.dat{1}.segment.verbose
        name = 's-resp';
        fh   = findobj('Type','Figure','Name',name);
        if ~isempty(fh)
            fig{1} = fh;
        else
            fig{1} = figure('Name',name,'NumberTitle','off'); 
        end
        clf(fig{1});   
        name = 's-bf';
        fh   = findobj('Type','Figure','Name',name);
        if ~isempty(fh)
            fig{2} = fh;
        else
            fig{2} = figure('Name',name,'NumberTitle','off'); 
        end
        clf(fig{2});   
        name = 's-tpm';
        fh   = findobj('Type','Figure','Name',name);
        if ~isempty(fh)
            fig{3} = fh;
        else
            fig{3} = figure('Name',name,'NumberTitle','off'); 
        end
        clf(fig{3});           
        name = 's-ll';
        fh   = findobj('Type','Figure','Name',name);
        if ~isempty(fh)
            fig{4} = fh;
        else
            fig{4} = figure('Name',name,'NumberTitle','off'); 
        end
        clf(fig{4});                   
end
if pars.niter>1
    if pars.verbose>3, 
        name = 'bf';
        fh   = findobj('Type','Figure','Name',name);
        if ~isempty(fh)
            fig{11} = fh;
        else
            fig{11} = figure('Name',name,'NumberTitle','off'); 
        end
        clf(fig{11});         
    end
    if pars.verbose>2, 
        name = 'def';
        fh   = findobj('Type','Figure','Name',name);
        if ~isempty(fh)
            fig{10} = fh;
        else
            fig{10} = figure('Name',name,'NumberTitle','off'); 
        end
        clf(fig{10});        
    end
    if pars.verbose>2, 
        name = 'll';
        fh   = findobj('Type','Figure','Name',name);
        if ~isempty(fh)
            fig{9} = fh;
        else
            fig{9} = figure('Name',name,'NumberTitle','off'); 
        end
        clf(fig{9});        
    end
    if pars.verbose>2, 
        name = 'prior';
        fh   = findobj('Type','Figure','Name',name);
        if ~isempty(fh)
            fig{8} = fh;
        else
            fig{8} = figure('Name',name,'NumberTitle','off'); 
        end
        clf(fig{8});         
    end
    if pars.verbose>2, 
        name = 'resps';
        fh   = findobj('Type','Figure','Name',name);
        if ~isempty(fh)
            fig{7} = fh;
        else
            fig{7} = figure('Name',name,'NumberTitle','off'); 
        end
        clf(fig{7});         
    end
    if pars.verbose>1, 
        name = 'tpm';
        fh   = findobj('Type','Figure','Name',name);
        if ~isempty(fh)
            fig{6} = fh;
        else
            fig{6} = figure('Name',name,'NumberTitle','off'); 
        end
        clf(fig{6});         
    end
    if pars.verbose>0, 
        name = 'L';
        fh   = findobj('Type','Figure','Name',name);
        if ~isempty(fh)
            fig{5} = fh;
        else
            fig{5} = figure('Name',name,'NumberTitle','off'); 
        end
        clf(fig{5});         
    end        
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
    rand_subjs{m} = 1:min(S0,S);%randperm(S,min(S0,S));   
    obj{m}        = cell(1,S);
    
    pars.dat{m}.segment.biasreg  = pars.dat{m}.segment.biasreg*ones(1,N);
    pars.dat{m}.segment.biasfwhm = pars.dat{m}.segment.biasfwhm*ones(1,N);     
    
    if V{1}(1).dim(3)==1
        % Is 2D
        %------------------------------------------------------------------
        pars.dat{m}.maff.do_maff = false;                             
    end
    
    pars.dat{m}.write_res.write_bf = repmat(pars.dat{m}.write_res.write_bf,N,1);  
        
    for s=1:S                
        obj1 = struct;
             
        %------------------------------------------------------------------
        dir_s      = fullfile(dir_subjects,['s-' num2str(s) '_m-' num2str(m)]);
        obj1.dir_s = dir_s;
        mkdir(dir_s);         
        
        if ~isempty(dir_local)
            dir_local_s = fullfile(dir_local,['s-' num2str(s) '_m-' num2str(m)]);
            mkdir(dir_local_s); 
        end        
        
        %------------------------------------------------------------------        
        obj1.image  = V{s};
        obj1.labels = pars.dat{m}.labels{s};
        
        %------------------------------------------------------------------ 
        obj1.do_template     = pars.do_template;
        obj1.modality        = pars.dat{m}.modality;
        obj1.pth_template    = pars.pth_template;                        
        obj1.est_fwhm        = pars.dat{m}.est_fwhm;                        
        obj1.dt              = pars.dt;        
        obj1.uniform         = pars.uniform;    
        obj1.dir_template    = pars.dir_template;    
        obj1.fwhm            = pars.dat{m}.fwhm;
        obj1.Affine          = pars.dat{m}.Affine;
        obj1.print_subj_info = pars.dat{m}.print_subj_info;
        
        obj1.s            = cnt;        
        obj1.ll_template     = 0;
        obj1.status          = 0;
        obj1.iter            = 1;
        
        obj1.address_local = pars.address_local;        
        
        %------------------------------------------------------------------
        obj1.maff = pars.dat{m}.maff;
                
        %------------------------------------------------------------------
        obj1.push_resp = pars.dat{m}.push_resp;
        
        %------------------------------------------------------------------                 
        obj1.write_res = pars.dat{m}.write_res;        
        if V{1}(1).dim(3)==1
            % Is 2D
            obj1.write_res.write_tc(:,2) = false; 
            obj1.write_res.write_df(:,:) = false; 
        end
        
        %------------------------------------------------------------------                     
        obj1.segment           = pars.dat{m}.segment;            
        obj1.segment.ll_all    = [];        
        obj1.segment.bf.dc     = zeros(1,N);
        obj1.segment.bf.avg_dc = zeros(1,N);        
                  
        % Build directory structure
        %------------------------------------------------------------------
        [~,nam] = fileparts(obj1.image(1).fname);

        dir_resp = fullfile(dir_s,'pushed-responsibilities'); 
        mkdir(dir_resp); 

        pth_resp = cell(1,pars.K);
        for k=1:pars.K
            fname       = fullfile(dir_resp,['push_resp-', num2str(k), '-' nam, '.nii']);
            pth_resp{k} = fname;  
        end
        obj1.pth_resp = pth_resp;
    
        dir_der      = fullfile(dir_s,'template-derivatives');      
        obj1.dir_der = dir_der;
        mkdir(dir_der);        
        
        if ~isempty(dir_local)
            dir_der_local      = fullfile(dir_local_s,'template-derivatives');      
            obj1.dir_der_local = dir_der_local;
            mkdir(dir_der_local);  
        else
            obj1.dir_der_local = '';
        end
        
        pth_gr = cell(1,d_gr(end));
        for k=1:d_gr(end)
            fname     = fullfile(dir_der,['gr-', num2str(k), '-' nam, '.nii']);
            pth_gr{k} = fname;  
        end
        obj1.pth_gr = pth_gr;

        if ~isempty(dir_local)
            pth_gr = cell(1,d_gr(end));
            for k=1:d_gr(end)
                fname     = fullfile(dir_der_local,['gr-', num2str(k), '-' nam, '.nii']);
                pth_gr{k} = fname;  
            end
            obj1.pth_gr_local = pth_gr;
        end
                
        pth_H = cell(1,d_H(end));
        for k=1:d_H(end)
            fname    = fullfile(dir_der,['H-', num2str(k), '-' nam, '.nii']);
            pth_H{k} = fname;  
        end
        obj1.pth_H = pth_H;

        if ~isempty(dir_local)
            pth_H = cell(1,d_H(end));
            for k=1:d_H(end)
                fname    = fullfile(dir_der_local,['H-', num2str(k), '-' nam, '.nii']);
                pth_H{k} = fname;  
            end
            obj1.pth_H_local = pth_H;
        end
                
        dir_vel = fullfile(dir_s,'deformations');   
        mkdir(dir_vel);

        obj1.pth_vel = fullfile(dir_vel,['deformation-' nam '.nii']);   

        dir_write = fullfile(dir_s,'results');   
        mkdir(dir_write);

        obj1.dir_write = dir_write;
        
        dir_seg = fullfile(dir_s,'segmentations'); 
        mkdir(dir_seg);
        
        obj1.dir_seg = dir_seg;
        
        %------------------------------------------------------------------
        obj{m}{s} = obj1;
        
        cnt = cnt + 1;
    end
end

% A random sample of subjects (used for visualising algorithm progress)
pars.rand_subjs = rand_subjs;
%==========================================================================
