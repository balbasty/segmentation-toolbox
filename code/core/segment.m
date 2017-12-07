function obj = segment(obj,fig)
if nargin<2, fig = cell(4,1); end
                               
obj.iter = obj.iter + 1;
fprintf('iter=%d, m=%d, s=%d\n',obj.iter,obj.m,obj.s);                  

if obj.iter==2
    obj.uniform = false;
end         

% Load template
logtpm = spm_load_logpriors8(obj.pth_logTPM,obj.tiny,obj.deg,obj.uniform);

if obj.doaff && ((obj.use_tpm && obj.iter==1) || (~obj.use_tpm && obj.iter==3))
    % Affinely register template to 1st image
    %----------------------------------------------------------------------
    
    obj.Affine = spm_logmaff8(obj.image(1),obj.samp,(obj.fwhm+1)*16,obj.descrip,logtpm,obj.Affine,'mni');            
    obj.Affine = spm_logmaff8(obj.image(1),obj.samp, obj.fwhm,      obj.descrip,logtpm,obj.Affine,'mni');                                    
    obj.doaff  = false;
end                 

if ~obj.use_tpm && obj.dotpm     
    % When estimating TPMs
    %----------------------------------------------------------------------
    
    if obj.iter==1   
        % 1st iteration

        obj.dodef = false; % No deformation update for iteration 1, 2 and 3
        obj.dowp  = false; % No tissue mixing weight update for iteration 1 and 2
    elseif obj.iter==2                     
        % 2nd iteration    

        obj = clear_pars(obj); % Re-estimate cluster and bias field parameters

    elseif obj.iter==3                  
        % 3rd iteration
        
        obj = clear_pars(obj); % Re-estimate cluster and bias field parameters                

        obj.dowp  = obj.dowp0;
        obj.dodef = obj.dodef0; 
    end 
end         

% Run segmentation algorithm
if ~isempty(fig{1})
    obj = spm_preprocx(obj,logtpm,fig);
else 
    try       
        obj        = spm_preprocx(obj,logtpm,fig);
        obj.status = 0;
    catch ME
        % Error!
        obj.status = 1;
        str_output = ['Error in: ' obj.image(1).fname '\n' ME.message '\n'];
        fprintf(str_output)
    end
end
%==========================================================================