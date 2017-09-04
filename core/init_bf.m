function chan = init_bf(N,C,V,d,lam,fwhm)
lam  = lam*ones(1,C);
fwhm = fwhm*ones(1,C);

[x0,y0,~] = ndgrid(1:d(1),1:d(2),1);
z0        = 1:d(3);
        
cl   = cell(C,1);
args = {'C',cl,'B1',cl,'B2',cl,'B3',cl,'T',cl};
chan = cell(N,1);
for n=1:N
    chan{n} = struct(args{:});
       
    for c=1:C
        vx = sqrt(sum(V{n,c}.mat(1:3,1:3).^2));
        d0 = V{n,c}.dim;
    
        % GAUSSIAN REGULARISATION for bias correction
        sd    = vx(1)*d0(1)/fwhm(c); d3(1) = ceil(sd*2); krn_x = exp(-(0:(d3(1)-1)).^2/sd.^2)/sqrt(vx(1));
        sd    = vx(2)*d0(2)/fwhm(c); d3(2) = ceil(sd*2); krn_y = exp(-(0:(d3(2)-1)).^2/sd.^2)/sqrt(vx(2));
        sd    = vx(3)*d0(3)/fwhm(c); d3(3) = ceil(sd*2); krn_z = exp(-(0:(d3(3)-1)).^2/sd.^2)/sqrt(vx(3));
        Cbias = kron(krn_z,kron(krn_y,krn_x)).^(-2)*lam(c);
        
        chan{n}(c).C = sparse(1:length(Cbias),1:length(Cbias),Cbias,length(Cbias),length(Cbias));    

        % Basis functions for bias correction
        chan{n}(c).B3 = spm_dctmtx(d(3),d3(3),z0);
        chan{n}(c).B2 = spm_dctmtx(d(2),d3(2),y0(1,:)');
        chan{n}(c).B1 = spm_dctmtx(d(1),d3(1),x0(:,1));

        % Initial parameterisation of bias field
        chan{n}(c).T = zeros(d3);
    end
end