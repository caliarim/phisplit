function [PV,phi_cache]=phisplit(tau,A,V,p,phi_cache)
%PHISPLIT phi-function of a Kronecker sum applied to a tensor
%   with direction splitting.
%   PHISPLIT(TAU,A,V,P) returns an ND-array PV which is an approximation of
%   order TAU^2 of
%
%   PHI_P(TAU*K)*V(:),
%
%   that is PV(:) = PHI_P(TAU*K)*V(:) + O(TAU^2). Here TAU is a scalar,
%   A is a cell of length D containing the square full matrices
%   A{1},..., A{D}, V is an ND-array, and P is a nonnegative integer.
%   K is a matrix corresponding to the Kronecker sum of A{D}, A{D-1},..., A{1},
%   but it is not assembled. If P is equal to 0 then the approximation is in
%   fact exact.
%
%   [___] = PHISPLIT([], [], V, P, PHI_CACHE) uses the information contained in
%   the structure PHI_CACHE to form the desired approximation.
%   If P == 0, PHISPLIT uses the field exp of the structure, which is
%   a 1xD cell array containing the matrices PHI_0(tau*A{MU}), that is
%
%   PHI_CACHE.exp{MU} = PHI_0(tau*A{MU}), MU = 1,..., D.
%
%   If P > 0, PHISPLIT uses the field phi of the structure, which is
%   a PxD cell array containing the matrices PHI_P(tau*A{MU}), that is
%
%   PHI_CACHE.phi{P,MU} = PHI_P(tau*A{MU}), MU = 1,..., D.
%
%   [PV, PHI_CACHE] = PHISPLIT(___)  also returns the structure PHI_CACHE as
%   given in input or as computed by PHISPLIT itself.
%
%   Notice that PHISPLIT requires the function TUCKER from the package
%   KronPACK which can be added to the path by
%
%   addpath('extern/KronPACK/src')
%
%   Example
%         tau = 1/4; A{1} = randn(3); A{2} = randn(4); V = randn(3,4); p = 1;
%         K = kron(eye(4),A{1})+kron(A{2},eye(3));
%         pv = PHISPLIT(tau,A,V,p);
%         ref = phiquad(tau*K,p)*V(:);
%         norm(pv(:)-ref,inf)
%         tau = 1/8;
%         pv = PHISPLIT(tau,A,V,p);
%         ref = phiquad(tau*K,p)*V(:);
%         norm(pv(:)-ref,inf)
%
%   See also PHIQUAD.
%
%   [CC23] M. Caliari and F. Cassini.
%          Direction splitting of phi-functions in exponential integrators
%          for d-dimensional problems in Kronecker form, Submitted 2023

  d = length(size(V));

  if ((nargin < 5) || isempty(phi_cache) || isempty(fieldnames(phi_cache)))
    % Compute needed phi functions
    if (p==0)
      phi_cache.exp = cell(1,d);
      for mu = 1:d
        phi_cache.exp{mu} = phiquad(tau*A{mu},0);
      end
    else
      phi_cache.phi = cell(p,d);
      for mu = 1:d
        phi_cache.phi{p,mu} = phiquad(tau*A{mu},p);
      end
    end
  end
  % Assemble approximation
  if (p==0)
    PV = tucker(V,phi_cache.exp);
  else
    PV = tucker((factorial(p)^(d-1))*V,phi_cache.phi{p,:});
  end
end
%!demo
%! taurange = [1e-1,1e-2,1e-3,1e-4];
%! prange = 1:4;
%! d = 3;
%! n = [2,3,4];
%! for mu = 1:d
%!   A{mu} = randn(n(mu));
%! end
%! V = randn(n);
%! v = V(:);
%! AA = full(kronsum(A));
%! pn = prod(n);
%! counterp = 0;
%! for p = prange
%!   counterp = counterp + 1;
%!   countertau = 0;
%!   for tau = taurange
%!     countertau = countertau + 1;
%!     PV = phisplit(tau,A,V,p);
%!     M = [tau*AA,v,zeros(pn,p-1);zeros(p,pn),diag(ones(1,p-1),1)];
%!     Pvref=expm(M)*[zeros(pn+p-1,1);1];
%!     Pvref = Pvref(1:pn);
%!     err(counterp,countertau) = norm(PV(:)-Pvref,inf)/norm(Pvref,inf);
%!   end
%! end
%! figure;
%! for i = 1:length(prange)
%!   subplot(2,2,i)
%!   loglog(taurange,err(i,:),'x',taurange,err(i,end)*(taurange/taurange(end)).^(2),'-')
%!   legend('Error','Reference rate 2')
%!   title(sprintf('p=%i',prange(i)))
%! end
%!test
%! tau = rand;
%! d = 3;
%! n = [2,3,4];
%! for mu = 1:d
%!   A{mu} = randn(n(mu));
%! end
%! V = randn(n);
%! p = 0;
%! PV = phisplit(tau,A,V,p);
%! v = V(:);
%! AA = full(kronsum(A));
%! PVref = expm(tau*AA)*v;
%! assert(PV(:),PVref,1e-10)
%!test
%! tau = 0.01;
%! d = 2;
%! n = [2,3];
%! for mu = 1:d
%!   A{mu} = randn(n(mu));
%! end
%! V = randn(n);
%! p = 1;
%! PV = phisplit(tau,A,V,p);
%! v = V(:);
%! AA = full(kronsum(A));
%! pn = prod(n);
%! M = [tau*AA,v,zeros(pn,p-1);zeros(p,pn),diag(ones(1,p-1),1)];
%! Pvref=expm(M)*[zeros(pn+p-1,1);1];
%! Pvref = Pvref(1:pn);
%! assert(PV(:),Pvref,10*tau^2)
%!test
%! tau = 0.001;
%! d = 4;
%! n = [2,3,4,5];
%! for mu = 1:d
%!   A{mu} = randn(n(mu));
%! end
%! V = randn(n);
%! p = 2;
%! PV = phisplit(tau,A,V,p);
%! v = V(:);
%! AA = full(kronsum(A));
%! pn = prod(n);
%! M = [tau*AA,v,zeros(pn,p-1);zeros(p,pn),diag(ones(1,p-1),1)];
%! Pvref=expm(M)*[zeros(pn+p-1,1);1];
%! Pvref = Pvref(1:pn);
%! assert(PV(:),Pvref,10*tau^2)
%!test
%! tau = 0.001;
%! d = 3;
%! n = [4,2,3];
%! for mu = 1:d
%!   A{mu} = randn(n(mu));
%! end
%! V = randn(n);
%! p = 3;
%! PV = phisplit(tau,A,V,p);
%! v = V(:);
%! AA = full(kronsum(A));
%! pn = prod(n);
%! M = [tau*AA,v,zeros(pn,p-1);zeros(p,pn),diag(ones(1,p-1),1)];
%! Pvref=expm(M)*[zeros(pn+p-1,1);1];
%! Pvref = Pvref(1:pn);
%! assert(PV(:),Pvref,10*tau^2)
%!test
%! tau = 0.01;
%! d = 3;
%! n = [2,3,4];
%! p = 3;
%! for mu = 1:d
%!   A{mu} = randn(n(mu));
%!   phi_cache.exp{mu} = expm(tau*A{mu});
%!   [phi_cache.phi{1:3,mu}] = phiquad(tau*A{mu},p);
%! end
%! V = randn(n);
%! v = V(:);
%! AA = full(kronsum(A));
%! pn = prod(n);
%! PV = phisplit(tau,A,V,0,phi_cache);
%! Pvref=expm(tau*AA)*v;
%! assert(PV(:),Pvref,1e-10)
%! for ell = 1:p
%!   PV = phisplit(tau,A,V,ell,phi_cache);
%!   M = [tau*AA,v,zeros(pn,ell-1);zeros(ell,pn),diag(ones(1,ell-1),1)];
%!   Pvref=expm(M)*[zeros(pn+ell-1,1);1];
%!   Pvref = Pvref(1:pn);
%!   assert(PV(:),Pvref,10*tau^2)
%! end
%!test
%! tau = 0.01;
%! d = 3;
%! n = [2,3,4];
%! p = 1;
%! for mu = 1:d
%!   A{mu} = randn(n(mu));
%! end
%! V = randn(n);
%! v = V(:);
%! AA = full(kronsum(A));
%! pn = prod(n);
%! phi_cache = [];
%! [PV,phi_cache] = phisplit(tau,A,V,1,phi_cache);
%! PV2 = phisplit(tau,A,V,1,phi_cache);
%! M = [tau*AA,v;zeros(1,pn+1)];
%! Pvref=expm(M)*[zeros(pn,1);1];
%! Pvref = Pvref(1:pn);
%! assert(PV2(:),Pvref,10*tau^2)
%!test
%! tau = 0.01;
%! d = 3;
%! n = [2,3,4];
%! p = 1;
%! for mu = 1:d
%!   A{mu} = randn(n(mu));
%! end
%! V = randn(n);
%! v = V(:);
%! AA = full(kronsum(A));
%! pn = prod(n);
%! [PV,phi_cache] = phisplit(tau,A,V,1);
%! PV2 = phisplit(tau,A,V,1,phi_cache);
%! M = [tau*AA,v;zeros(1,pn+1)];
%! Pvref=expm(M)*[zeros(pn,1);1];
%! Pvref = Pvref(1:pn);
%! assert(PV2(:),Pvref,10*tau^2)
