function [U,err_ss]=lawsoneuler(U0,m,tstar,Kcell,G,Uss)
% LAWSONEULER Lawson--Euler integrator.
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
      [U,phi_cache] = phisplit(tau,Kcell,U + tau*G(t,U),0,phi_cache);
      err_ss(i+1) = norm(Uss-U,'fro')/normUss;
      t = t + tau;
    end
  else
    for i = 1:m
      [U,phi_cache] = phisplit(tau,Kcell,U + tau*G(t,U),0,phi_cache);
      t = t + tau;
    end
  end
end
