function [U,err_ss] = exprk2_phisplit(U0,m,tstar,Kcell,g,Uss)
% EXPRK2_PHISPLIT ETD2RK integrator using PHISPLIT.
%
% See also EXAMPLE_LQC, EXAMPLE_ADR.

  tau = tstar/m;
  U = U0;
  t = 0;
  ndir = length(size(U0));
  phi_cache.phi = cell(2,ndir);
  for mu = 1:ndir % we could just one mu for linear quadratic control,
                  % but does not matter
    [phi_cache.phi{1:2,mu}] = phiquad(tau*Kcell{mu},2);
  end

  if (nargin == 6)
    normUss = norm(Uss,'fro');
    err_ss(1) = norm(Uss-U,'fro')/normUss;
    for i = 1:m
      gn = g(t,U);
      Fn = kronsumv(U,Kcell)+ gn;
      aux = phisplit(tau,Kcell,Fn,1,phi_cache);
      U2 = U+tau*aux;
      D2 = g(t+tau,U2) - gn;
      U = U + tau*(aux+phisplit(tau,Kcell,D2,2,phi_cache));
      err_ss(i+1) = norm(Uss-U,'fro')/normUss;
      t = t + tau;
    end
  else
    for i = 1:m
      gn = g(t,U);
      Fn = kronsumv(U,Kcell)+ gn;
      aux = phisplit(tau,Kcell,Fn,1,phi_cache);
      U2 = U+tau*aux;
      D2 = g(t+tau,U2) - gn;
      U = U + tau*(aux+phisplit(tau,Kcell,D2,2,phi_cache));
      t  = t + tau;
    end
  end
end
