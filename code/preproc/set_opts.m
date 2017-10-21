function var = set_opts(opts,field,default,mn,mx)
% Set options with default, min and max defined
% FORMAT var = set_opts(opts,field,default,mn,mx)
%__________________________________________________________________________

  var = default;
  
  if isempty(opts)
    return
  end
  
  % has the option already been set?
  if ~isfield(opts,field) 
      % see if there is a capitalization problem:
      names = fieldnames(opts);
      for i = 1:length(names)
          if strcmpi(names{i},field)
              opts.(field) = opts.(names{i});
              opts=rmfield(opts,names{i});
              break;
          end
      end
  end

  if isfield(opts,field) && ~isempty(opts.(field))
      var = opts.(field);  % override the default
  end

  % perform error checking, if desired
  if nargin >= 4 && ~isempty(mn)
      if var < mn
          fprintf('Variable %s is %f, should be at least %f\n',...
              field,var,mn); error('variable out-of-bounds');
      end
  end
  if nargin >= 5 && ~isempty(mx)
      if var > mx
          fprintf('Variable %s is %f, should be less than %f\n',...
              field,var,mx+1); error('variable out-of-bounds');
      end
  end
  opts.(field) = var;