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