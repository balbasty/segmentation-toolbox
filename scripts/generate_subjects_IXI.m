close all; clc;

addpath('../util/')
addpath('../core/')
addpath('../preproc/')

% Parameters
datadir  = '/home/mbrud/Data/IXI-clean/';
regdir   = '/home/mbrud/Data/IXI-subjects/';

reg   = false;
slice = true;

%% Get triplets of same subject images
d = dir(datadir);
d = d(3:end);
D = numel(d);

fnames = cell(D,1);
for i=1:D
    fnames{i} = dir(fullfile(datadir,d(i).name,'*.nii'));
end

%%
P = {};
for i1=1:numel(fnames{1})
    
    tag1 = fnames{1}(i1).name(1:6);

    Ptmp{1} = fullfile(datadir,d(1).name,fnames{1}(i1).name);
    exists = 0;
    for i2=2:D        
        for i3=1:numel(fnames{i2})

            tag2 = fnames{i2}(i3).name(1:6);

            if strcmp(tag1,tag2)
               exists = 1; 
               
               Ptmp{end + 1} = fullfile(datadir,d(i2).name,fnames{i2}(i3).name);
               
               break;
            end
        end
    end
    
    if exists
        P{end + 1} = Ptmp;
    end
    
    Ptmp = {};
end

%% Co-register triplets
if (exist(regdir,'dir') == 0)
    mkdir(regdir);
end

prefix = '';
N      = numel(P);
for i=1:N
    Ptmp = {};

    ref = P{i}(1);

    Ptmp{end + 1} = ref{1};

    for i2=1:numel(P{i}) - 1
        source = P{i}(i2 + 1);

        [pth,nam,ext] = fileparts(source{1});
        Ptmp{end + 1} = fullfile(pth,[prefix nam ext]); 
    end    

    f = fullfile(regdir,['S' num2str(i)]);

    if exist(f,'dir')
        rmdir(f,'s');
    end
    if (exist(f,'dir') == 0)
        mkdir(f);
    end 

    for i2=1:D           
        copyfile(Ptmp{i2},f);
    end    
end