function save_in_parfor(pth,var,is_struct)
if nargin<3, is_struct=''; end
if strcmp(is_struct,'-struct')
    save(pth,'-struct','var')
else
    save(pth,'var') 
end
    