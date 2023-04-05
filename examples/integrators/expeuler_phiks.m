function [U,err_ss]=expeuler_phiks(U0,m,tstar,Kcell,F,kappa,Uss)
% EXPEULER_PHIKS Exponential Euler integrator using PHIKS.
%
% See also EXAMPLE_ADR.

  tau = tstar/m;

  U = U0;
  t = 0;

  if (nargin == 7)
    normUss = norm(Uss,'fro');
    err_ss(1) = norm(Uss-U,'fro')/normUss;
    for i = 1:m
      tol = 1e-14 + kappa*tau^2*norm(U(:));
      aux = phiks(tau,Kcell,F(t,U),1,tol);
      U = U + tau*aux{2};
      err_ss(i+1) = norm(Uss-U,'fro')/normUss;
      t = t + tau;
    end
  else
    for i = 1:m
      tol = 1e-14 + kappa*tau^2*norm(U(:));
      aux = phiks(tau,Kcell,F(t,U),1,tol);
      U = U + tau*aux{2};
      t = t + tau;
    end
  end
end
