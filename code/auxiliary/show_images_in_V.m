function show_images_in_V(V)
S   = numel(V);
N   = numel(V{1});
nV  = {};
cnt = 1;
for s=1:S
    for n=1:N
        nV{cnt} = V{s}(n).fname;
        cnt     = cnt + 1;    
    end
end

S1 = min(N*S,20);
p  = randperm(N*S,S1);
nV = nV(p);

spm_check_registration(char(nV));
