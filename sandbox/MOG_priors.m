
%Set up some random Gaussian-Wishart hyperparameters, for each of N subjects
N=10;
D=3;
W=zeros(D,D,N);
m=zeros(D,N);
b=zeros(N,1);
n=zeros(N,1);
tmp=randn(D); tmp=tmp'*tmp;
for i=1:N,
    m(:,i)   = randn(D,1)*2+[1:D]';
    b(i)     = exp(randn(1)+2)+1;
    n(i)     = exp(randn(1)+2)+D;
    X        = randn(D);
    W(:,:,i) = n(i)*(X'*X*0.1+tmp);
end




%______________________________________________________________________________
%
% Starting estimates
m0 = zeros(D,1);
b0 = 1;
n0 = D-.999;
W0 = eye(D,D);
%______________________________________________________________________________


%______________________________________________________________________________
%
% Compute m_0

g = zeros(D,1);
H = zeros(D,D);
for i=1:N
    g = g + b0*n(i)*W(:,:,i)*(m(:,i)-m0);
    H = H + b0*n(i)*W(:,:,i);
end
m0 = m0 + H\g;
%______________________________________________________________________________


%______________________________________________________________________________
%
% Compute \beta_0

g_const = 0;
for i=1:N
    g_const = g_const - 0.5*(D/b(i) + n(i)*(m(:,i)-m0)'*W(:,:,i)*(m(:,i)-m0));
end
for subit=1:100,
    % Diff w.r.t. b0
    g  = 0.5*N*D/b0 + g_const; % Gradient
    H  = 0.5*N*D/b0^2;         % Hessian
    b0 = max(b0 + H\g,1e-5);
    if norm(g)==0, break; end
end
%______________________________________________________________________________


%______________________________________________________________________________
%
% Set up some constants

nW_const = zeros(D);
for i=1:N
    nW_const = nW_const + n(i)*W(:,:,i);
end

ElogLam = N*D*log(2);
for i=1:N
    ElogLam = ElogLam + 2*sum(log(diag(chol(W(:,:,i)))));
    for j=1:D
        ElogLam = ElogLam + psi((n(i)+1-j)/2);
    end
end

convergence = [];
E = -realmax;

for it=1:1000
    %______________________________________________________________________________
    %
    % Compute objective function (Equation 10.74 of Bishop)

    oE = E;
    logB = -n0*sum(log(diag(chol(W0)))) - n0*D/2*log(2) - D*(D-1)/4*log(pi);
    for j=1:D, logB = logB - gammaln((n0+1-j)/2); end
    E = (0.5*D*log(b0/(2*pi)) + logB)*N + 0.5*(n0-D-1)*ElogLam;
    for i=1:N
        e = 0.5*(-D*b0/b(i) - b0*n(i)*(m(:,i)-m0)'*W(:,:,i)*(m(:,i)-m0))...
          - 0.5*n(i)*trace(W0\W(:,:,i));
        E = E + e;
    end
    %if E-oE<abs(E)*eps*D^2, break; end
    if E-oE==0, break; end

    convergence = [convergence E];
    plot(convergence,'.-'); drawnow;


    %______________________________________________________________________________
    %
    % Compute \nu_0

    % Objective function terms containing n0:
    % NlogB = -n0*N*(sum(log(diag(chol(W0)))) + D/2*log(2));
    % for j=1:D, NlogB = NlogB - N*gammaln((n0+1-j)/2); end
    % E = NlogB + n0*0.5*ElogLam

    g = (sum(log(diag(chol(W0)))) + D/2*log(2))*N - 0.5*ElogLam;
    H = 0;
    for j=1:D
        g = g + N*psi(  (n0+1-j)/2)/2;
        H = H + N*psi(1,(n0+1-j)/2)/4;
    end
    n0 = max(n0 - H\g,D-0.99999);
    %______________________________________________________________________________



    %______________________________________________________________________________
    %
    % Compute W_0

    % Objective function terms containing W0:
    % E = -n0*N*sum(log(diag(chol(W0))));
    % for i=1:N
    %    E = E - 0.5*n(i)*trace(W0\W(:,:,i));
    % end

    C = inv(chol(W0));

    % Objective function terms containing W0, after
    % re-expressing using C = inv(chol(W0)):
    % E = n0*N*sum(log(diag(C)));
    % for i=1:N
    %    E = E - 0.5*n(i)*trace(C'*W(:,:,i)*C);
    % end

    G  = -n0*N*diag(1./diag(C)) + nW_const*C;
    for d=1:D,
        c        = C(1:d,d);
        g        = G(1:d,d);
        H        = nW_const(1:d,1:d);
        H(d,d)   = H(d,d) + n0*N/c(d)^2;
        C(1:d,d) = c - H\g;
    end
    C  = inv(C);
    W0 = C'*C;

end

% Show answers:
b0
m0
n0
W0




if false
    % For computing 1st and 2nd derivatives for updating the cholesky decomposition of W_0
    syms nn
    syms c11 c12 c13 c21 c22 c23 c31 c32 c33
    syms f11 f12 f13 f21 f22 f23 f31 f32 f33
    C=[c11 c12 c13; 0 c22 c23; 0 0 c33];
    F=[f11 f12 f13; f12 f22 f23; f13 f23 f33];

    f = find([1 0 0  1 1 0  1 1 1]);
    H=sym(zeros(6));
    g=sym(zeros(6,1));
    for i=1:6
        g(i) = simplify(diff(-nn*sum(log(diag(C))) + 0.5*trace(C.'*F*C),C(f(i))));
        for j=1:6
            H(i,j) = simplify(diff(diff(-nn*sum(log(diag(C))) + 0.5*trace(C.'*F*C),C(f(i))),C(f(j))));
        end
    end

end


