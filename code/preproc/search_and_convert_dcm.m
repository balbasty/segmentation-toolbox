function P = search_and_convert_dcm(dicomdir,outdir)
% Searches folder structure for sets of DICOM files and converts these to NIFTI format.
% FORMAT P = search_and_convert_dcm(dicomdir,outdir)
%        dicomdir - Root directory of DICOM data
%        outdir - Where to write NIFTI files
%        P - cell array of NIFTI filenames

seen = {}; % cell array to store seen edges
dcm_DFS(dicomdir,[],seen,outdir);

P = dir(fullfile(outdir,'*.nii'));
P = {P(:).name};
P = cellfun(@(c) fullfile(outdir,c),P,'uni',false);
P = P';
%==========================================================================

%==========================================================================
function seen = dcm_DFS(R,v,seen,niidir)
% Recursive depth-first search (https://en.wikipedia.org/wiki/Depth-first_search)

if isempty(v)
    v = R;
end

% Look for DICOMs
P = [dir(fullfile(v,'*.dcm')); dir(fullfile(v,'*.DCM'))];
if ~isempty(P)
    % There are DICOM files in the folder
    P = {P(:).name};
    P = cellfun(@(c) fullfile(v,c),P,'uni',false);

    % Read DICOMs
    H   = spm_dicom_headers(char(P));
    spm_dicom_convert_new(H,'all','flat','nii',niidir); % new version, which adjusts for variable slice thickness
end 

seen{end + 1,1} = v;                        % label v as discovered

w = dirs_in_dir(v);                         % get adjacents edges to v

w(ismember(w,v)) = [];                      % remove edge v

for i=1:numel(w)
    if ~ismember(seen,w{i})                 % if w not labeled as discovered
        seen = dcm_DFS(R,w{i},seen,niidir); % call dir_DFS recursively
    end
end
%==========================================================================

%==========================================================================
function dirs = dirs_in_dir(P)
% Get all subfolders of folder

d    = dir(P);
isub = [d(:).isdir];
dirs = {d(isub).name}';

dirs(ismember(dirs,{'.','..'})) = [];

for i=1:numel(dirs)
   dirs{i} = fullfile(P,dirs{i});
end
%==========================================================================
