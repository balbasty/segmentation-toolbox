function V = copy2tmpdir(V,tmpdir)
S = numel(V);
D = numel(V{1});

if exist(tmpdir,'dir')
    rmdir(tmpdir,'s');
end
mkdir(tmpdir);

for s=1:S
    fname = V{s}.fname;
    pth   = fileparts(fname);
    
    splitstr = strsplit(pth,'/');
    
    f = fullfile(tmpdir,splitstr{end});
    if exist(f,'dir')
        rmdir(f,'s');
    end
    mkdir(f);

    copyfile(pth,f);
    
    for d=1:D
        [~,nam,ext] = fileparts(V{s}(d).fname);
        V{s}(d)     = spm_vol(fullfile(f,[nam ext]));
    end
end
%=======================================================================
