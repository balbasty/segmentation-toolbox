%%
mex -v COPTIMFLAGS="-O3"  spm_vbmrf.c

%%
mex -v COPTIMFLAGS="-O3"  spm_vbmrf_lowerbound.c

%%
mex -v COPTIMFLAGS="-O3"  spm_vbmrf_update_Upsilon.c

%%
mex -v COPTIMFLAGS="-O3"  spm_vbmrf_all.c