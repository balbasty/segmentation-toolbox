function segment_holly(t)

pth0 = read_directory_details('directory_details.txt');    

pth     = fullfile(pth0,'pth_obj.mat');
tmp     = load(pth,'-mat');
pth_obj = tmp.pth_obj1;

pth   = fullfile(pth0,'split.mat');
tmp   = load(pth,'-mat');
split = tmp.split;

for t=1 + (t - 1)*split:t*split;
    if t>numel(pth_obj)
        return
    end
    
    obj = load(pth_obj{t},'-mat');   

    obj = segment(obj);

    save(pth_obj{t},'-struct','obj')
end
%==========================================================================
