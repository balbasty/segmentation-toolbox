function pars = read_images_segment(pars,verbose)
if nargin<2, verbose = true; end

M = numel(pars.dat);
for m=1:M % Loop over populations
    pth        = pars.dat{m}.dir_data;
    S          = pars.dat{m}.S;
    
    dir_subj   = dir(pth);
    dir_subj   = dir_subj(3:end);    
    dir_subj   = dir_subj([dir_subj.isdir]);
    S1         = numel(dir_subj);
    if S>S1, S = S1; end
    
    pars.dat{m}.S = S;  
    
    dir_channels = dir(fullfile(pth,dir_subj(1).name,'scans'));  
    dir_channels = dir_channels(3:end);            
    dir_channels = dir_channels([dir_channels.isdir]);
    N            = numel(dir_channels);

    labels = cell(1,S);
    scans  = cell(1,S);
    for s=1:S % Loop over subjects
        dir_channels = dir(fullfile(pth,dir_subj(s).name,'scans'));  
        dir_channels = dir_channels(3:end);            
        dir_channels = dir_channels([dir_channels.isdir]);

        for n=1:N % Loop over modalities
            dir_scans = dir(fullfile(pth,dir_subj(s).name,'scans',dir_channels(n).name,'*.nii'));                
            if isempty(dir_scans)
                dir_scans = dir(fullfile(pth,dir_subj(s).name,'scans',dir_channels(n).name,'*.img'));
            end

           scans{s}(n) = spm_vol(fullfile(pth,dir_subj(s).name,'scans',dir_channels(n).name,dir_scans(1).name));
        end

        dir_labels = fullfile(pth,dir_subj(s).name,'labels');
        labels1    = dir(fullfile(dir_labels,'*.nii'));
        if isempty(labels1)
            labels1 = dir(fullfile(dir_labels,'*.img'));
        end         

        if ~isempty(labels1)
            labels{s} = spm_vol(fullfile(dir_labels,labels1.name));
        end
    end
        
    pars.dat{m}.V      = scans;  
    pars.dat{m}.labels = labels;    
    
    if verbose
        fprintf('Loaded data from %d subject(s) having %d channel(s) each (%s)\n',S,N,pth); 
    end
end
%==========================================================================