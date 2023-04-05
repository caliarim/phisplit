function [U,err_ss] = expeuler_phisplit(U0,m,tstar,Kcell,F,Uss)
% EXPEULER_PHISPLIT Exponential Euler integrator using PHISPLIT.
%
% See also EXAMPLE_ADR.

  tau = tstar/m;
  U = U0;
  t = 0;
  phi_cache = [];

  if (nargin == 6)
    normUss = norm(Uss,'fro');
    err_ss(1) = norm(Uss-U,'fro')/normUss;
    for i = 1:m
      [aux,phi_cache] = phisplit(tau,Kcell,F(t,U),1,phi_cache);
      U = U + tau*aux;
      err_ss(i+1) = norm(Uss-U,'fro')/normUss;
      t = t + tau;
    end
  else
    for i = 1:m
      [aux,phi_cache] = phisplit(tau,Kcell,F(t,U),1,phi_cache);
      U = U + tau*aux;
      t = t + tau;
    end
  end
end
