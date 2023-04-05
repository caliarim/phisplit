function [U,err_ss]=exprbeuler_lqc_phisplit(U0,m,tstar,F,JF,Uss)
% EXPRBEULER_LQC_PHISPLIT Exponential Rosenbrock--Euler integrator
% for linear quadratic control example using PHISPLIT.
%
% See also EXAMPLE_LQC.

tau = tstar/m;
  U = U0;
  t = 0;
  phi_cache.phi = cell(1,2);
  if (nargin == 6)
    normUss = norm(Uss,'fro');
    err_ss(1) = norm(Uss-U,'fro')/normUss;
    for i = 1:m
      Jmat{1} = JF{1}(t,U);
      Jmat{2} = Jmat{1};
      phi_cache.phi{1,1} = phiquad(tau*Jmat{1},1,3);
      phi_cache.phi{1,2} = phi_cache.phi{1,1};
      U = U + tau*phisplit(tau,Jmat,F(t,U),1,phi_cache);
      err_ss(i+1) = norm(Uss-U,'fro')/normUss;
      t = t + tau;
    end
  else
    for i = 1:m
      Jmat{1} = JF{1}(t,U);
      Jmat{2} = Jmat{1};
      phi_cache.phi{1,1} = phiquad(tau*Jmat{1},1,3);
      phi_cache.phi{1,2} = phi_cache.phi{1,1};
      U = U + tau*phisplit(tau,Jmat,F(t,U),1,phi_cache);
      t = t + tau;
    end
  end
end
