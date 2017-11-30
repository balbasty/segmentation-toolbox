function obj = segment(obj)

if obj.status~=0
    fprintf('s=%d --> obj.status~=0',obj.s)
    return
end                        

if obj.diff_TPM                
    obj.iter = obj.iter + 1;
    fprintf('iter=%d, m=%d, s=%d (TPM)\n',obj.iter,obj.m,obj.s);                  

    if obj.iter==2
        obj.uniform = 0;
    end
else
    fprintf('iter=%d, m=%d, s=%d (prior)\n',obj.iter,obj.m,obj.s);  
end          

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
    obj.dodef = 0; % No deformation update for iteration 1, 2 and 3
    obj.dowp  = 0; % No tissue mixing weight update for iteration 1 and 2
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
    
    obj.dowp = obj.dowp0;
elseif ~obj.use_tpm && obj.iter==4 && obj.dotpm       
    % 4th iteration: start estimating deformations
    %----------------------------------------------------------------------
    obj.dodef = obj.dodef0; 
end         

% Run segmentation algorithm (spm_preprocx)
obj = spm_preprocx(obj,logtpm);
% obj.munum  = 0;
% obj.muden  = 0;
% obj.status = 0; % All OK   
% try       
%     obj = spm_preprocx(obj,logtpm);
% catch ME
%     % Error!
%     obj.status = 1;
%     str_output = ['Error in: ' obj.image(1).fname '\n' ME.message '\n'];
%     fprintf(str_output)
% end
%==========================================================================