function nw = get_nbr_of_cores
c  = parcluster('local');
nw = c.NumWorkers;