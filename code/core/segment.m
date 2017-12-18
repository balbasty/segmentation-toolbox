function obj = segment(obj,fig)
if nargin<2, fig = cell(4,1); end
                               
obj.iter = obj.iter + 1;
             
% Load template
logtpm = spm_load_logpriors8(obj.pth_logtpm,obj.deg);

if obj.doaff && ((obj.use_tpm && obj.iter==1) || (~obj.use_tpm && obj.iter==3))
    % Affinely register template to 1st image
    %----------------------------------------------------------------------
    
    obj.Affine = spm_logmaff8(obj.image(1),obj.samp,(obj.fwhm+1)*16,obj.descrip,logtpm,obj.Affine,'mni');            
    obj.Affine = spm_logmaff8(obj.image(1),obj.samp, obj.fwhm,      obj.descrip,logtpm,obj.Affine,'mni');  
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
                
        % Re-estimate cluster parameters        
        if isfield(obj.mog,'po'), obj.mog = rmfield(obj.mog,'po'); end 
        if isfield(obj.mog,'mn'), obj.mog = rmfield(obj.mog,'mn'); end 
        if isfield(obj.mog,'vr'), obj.mog = rmfield(obj.mog,'vr'); end 
    elseif obj.iter>=3                  
        % 3rd iteration

        obj.dowp  = obj.dowp0;
        obj.dodef = obj.dodef0; 
    end 
end         

% Run segmentation algorithm
try       
    obj        = spm_preprocx(obj,logtpm,fig);
    obj.status = 0;
catch ME
    % Error!        
    fprintf(['Error for image: ' obj.image(1).fname '\n'])
    for i=1:numel(ME.stack)
        disp([ME.stack(i).name ', line ' num2str(ME.stack(i).line)]);
    end
    disp(ME.message)
    obj.status = 1;
end
%==========================================================================