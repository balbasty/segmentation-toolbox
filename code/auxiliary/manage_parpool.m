function manage_parpool(num_workers)

nw = nbr_parfor_workers;
if num_workers>nw
    num_workers = nw;
end
    
poolobj = gcp('nocreate');
if ~isempty(poolobj) && num_workers==0
    delete(poolobj);
elseif ~isempty(poolobj) && poolobj.NumWorkers~=num_workers
    delete(poolobj);
    parpool('local',num_workers);
elseif isempty(poolobj) && num_workers
    parpool('local',num_workers);
end
%==========================================================================