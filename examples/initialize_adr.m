function [U0,A,g,F,odefun,odeJac]=initialize_adr(n,epsilon,alpha)
% INITIALIZE_ADR Auxiliary function for example of
% advection--diffusion--reaction equation.
%
% See also EXAMPLE_ADR.

  d = length(n);
  for mu = 1:d
    x{mu} = linspace(0,1,n(mu)+2).';
    x{mu} = x{mu}(2:n(mu)+1);
    h(mu) = 1/(n(mu)+1);
    D2{mu} = spdiags(ones(n(mu),1)*[1,-2,1]/h(mu)^2,-1:1,n(mu),n(mu));
    D1{mu} = spdiags(ones(n(mu),1)*[-1,1]/2/h(mu),[-1,1],n(mu),n(mu));
    A{mu} = epsilon*D2{mu}+alpha(mu)*D1{mu};
  end

  K = kronsum(A);
  for mu = 1:d
    A{mu} = full(A{mu});
  end
  [X{1:d}] = ndgrid(x{1:d});
  U0 = ones(n);
  for mu = 1:d
    U0 = U0.*X{mu}.*(1-X{mu})*4;
  end
  Phi1 = zeros(size(X{1}));
  Phi2 = zeros(size(X{1}));
  for mu = 1:d
    aux1 = alpha(mu)*(1-2*X{mu})*4;
    aux2 = epsilon*(-2)*4;
    for eta = [1:mu-1,mu+1:d]
      aux1 = aux1.*X{eta}.*(1-X{eta})*4;
      aux2 = aux2.*X{eta}.*(1-X{eta})*4;
    end
    Phi1 = Phi1+aux1;
    Phi2 = Phi2+aux2;
  end
  g = @(t,U) 1./(1+U.^2)+(U0-Phi2-Phi1)*exp(t)-1./(1+U0.^2*exp(2*t));
  dgdu = @(t,U) -2*U(:)./(1+U(:).^2).^2;
  F = @(t,U) kronsumv(U,A) + g(t,U);
  gv = @(t,U) 1./(1+U(:).^2)+(U0(:)-Phi2(:)-Phi1(:))*exp(t)-...
       1./(1+U0(:).^2*exp(2*t));
  odefun = @(t,u) K*u + gv(t,u);
  prodn = prod(n);
  odeJac = @(t,u) K + spdiags(dgdu(t,u),0,prodn,prodn);
end
