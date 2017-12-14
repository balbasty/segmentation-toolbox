function mom = mom_struct(K,N)
tiny = eps*eps;
mom  = struct('ind',[],'s0',tiny,'s1',[],'S2',[]);
for n=1:2^N       
    mom(n).ind = dec2bin(n-1,N)=='1'; % Inh/sum(h)/mddices
    Ni         = sum(mom(n).ind);
    mom(n).s0  = zeros(1,K) + tiny;   % Zeroeth moments
    mom(n).s1  = zeros(Ni,K);         % First moments
    mom(n).S2  = zeros(Ni,Ni,K);      % Second moments
end
%=======================================================================