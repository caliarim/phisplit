function [U,err_ss]=exprbeuler_phiks(U0,m,tstar,F,JF,kappa,Uss)
% EXPRBEULER_PHIKS Exponential Rosenbrock--Euler integrator using PHIKS.
%
% See also EXAMPLE_LQC.

  tau =tstar/m;
  U = U0;
  t = 0;
  if (nargin == 7)
    normUss = norm(Uss,'fro');
    err_ss(1) = norm(Uss-U,'fro')/normUss;
    for i = 1:m
      tol = 1e-14 + kappa*tau^3*norm(U(:)); % norm(U0(:)) = 0
      aux = phiks(tau,{JF{1}(t,U),JF{2}(t,U)},F(t,U),1,tol);
      U = U + tau*aux{2};
      err_ss(i+1) = norm(Uss-U,'fro')/normUss;
      t = t + tau;
    end
  else
    ndir = length(size(U0));
    for i = 1:m
      for mu = 1:ndir
        Jmat{mu} = JF{mu}(t,U);
      end
      tol = 1e-14 + kappa*tau^3*norm(U(:));
      aux = phiks(tau,Jmat,F(t,U),1,tol);
      U = U + tau*aux{2};
      t = t + tau;
    end
  end
end
