function browse_subjects(folder,S)
if nargin<2, S = Inf; end

[V,S,N] = read_ims_from_dir(folder,S);

fname = cell(1,N*S);
cnt   = 1;
for s=1:S
    for n=1:N
        fname{cnt} = V{s}(n).fname;
        cnt        = cnt + 1;
    end
end

spm_check_registration(char(fname));
%==========================================================================