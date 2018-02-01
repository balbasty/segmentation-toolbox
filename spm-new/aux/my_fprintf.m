function count = my_fprintf(varargin)
verbose           = varargin{end};
if verbose, count = fprintf(varargin{1:end - 1});
else        count = 0;
end
%=======================================================================