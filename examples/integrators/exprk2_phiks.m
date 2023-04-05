function [U,err_ss] = exprk2_phiks(U0,m,tstar,Kcell,g,kappa,Uss)
% EXPRK2_PHIKS ETD2RK integrator using PHIKS.
%
% See also EXAMPLE_LQC, EXAMPLE_ADR.

  tau =tstar/m;
  U = U0;
  t = 0;
  if (nargin == 7)
    normUss = norm(Uss,'fro');
    err_ss(1) = norm(Uss-U,'fro')/normUss;
    for i = 1:m
      tol = 1e-14 + kappa*tau^3*norm(U(:));
      gn = g(t,U);
      Fn = kronsumv(U,Kcell)+ gn;
      aux = phiks(tau,Kcell,Fn,1,tol);
      U2 = U + tau*aux{2};
      D2 = g(t+tau,U2) - gn;
      aux2 = phiks(tau,Kcell,D2,2,tol);
      U = U + tau*(aux{2}+aux2{3});
      err_ss(i+1) = norm(Uss-U,'fro')/normUss;
      t = t + tau;
    end
  else
    for i = 1:m
      tol = 1e-14 + kappa*tau^3*norm(U(:));
      gn = g(t,U);
      Fn = kronsumv(U,Kcell)+ gn;
      aux = phiks(tau,Kcell,Fn,1,tol);
      U2 = U + tau*aux{2};
      D2 = g(t+tau,U2) - gn;
      aux2 = phiks(tau,Kcell,D2,2,tol);
      U = U + tau*(aux{2}+aux2{3});
      t  = t + tau;
    end
  end
end
