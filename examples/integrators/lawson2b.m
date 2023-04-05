function [U,err_ss]=lawson2b(U0,m,tstar,Kcell,G,Uss)
% LAWSON2B Lawson2b integrator.
%
% See also EXAMPLE_ADR.

  tau = tstar/m; 

  U = U0;
  t = 0;
  phi_cache = [];

  if (nargin == 6)
    normUss = norm(Uss(:),inf);
    err_ss(1) = norm(Uss(:)-U(:),inf)/normUss;
    for i = 1:m
      [aux, phi_cache] = phisplit(tau,Kcell,U,0,phi_cache);
      aux2 = phisplit(tau,Kcell,G(t,U),0,phi_cache);
      U2 = aux + tau*aux2;
      U = aux + tau/2*(aux2+G(t+tau,U2));
      err_ss(i+1) = norm(Uss(:)-U(:),inf)/normUss;
      t = t + tau;
    end
  else
    for i = 1:m
      [aux, phi_cache] = phisplit(tau,Kcell,U,0,phi_cache);
      aux2 = phisplit(tau,Kcell,G(t,U),0,phi_cache);
      U2 = aux + tau*aux2;
      U = aux + tau/2*(aux2+G(t+tau,U2));
      t = t + tau;
    end
  end
end
