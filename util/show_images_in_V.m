function show_images_in_V(V)

S   = numel(V);
nV  = {};
cnt = 1;
for s=1:S
    for n=1:numel(V{s})
        nV{cnt} = V{s}(n).fname;
        cnt     = cnt + 1;    
    end
end

S1 = min(S,20);
p  = randperm(S,S1);
nV = nV(p);

spm_check_registration(char(nV));
