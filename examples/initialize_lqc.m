function [U0,Atcell,G,F,JFcell,odefun,odeJac,Uss]=initialize_lqc(nhat)
% INITIALIZE_LQC Auxiliary function for example of linear quadratic
% control.
%
% See also EXAMPLE_LQC.

  xl = 0;
  xr = 1;
  alpha = 100;
  d = nhat*nhat;

  U0 = zeros(d,d);

  x = linspace(xl,xr,nhat+2).';
  x = x(2:nhat+1);
  y = x;

  h = (xr-xl)/(nhat+1);

  D1 = spdiags(ones(nhat,1)*([-1,0,1]/(2*h)),-1:1,nhat,nhat);
  D2 = spdiags(ones(nhat,1)*([1,-2,1]/(h^2)),-1:1,nhat,nhat);

  Ax  = D2 + spdiags(-10*x,0,nhat,nhat)*D1;
  Ay  = D2 + spdiags(-100*y,0,nhat,nhat)*D1;

  b = repmat((x>0.1).*(x<=0.3),nhat,1);
  c = repmat(((x>0.7).*(x<=0.9)).',1,nhat);

  A = kronsum(Ax,Ay);
  Afull = full(A);
  At = A.';
  Atfull = full(At);
  Atcell = {Atfull,Atfull};

  C = alpha*(c.'*c);
  Cvec = C(:);
  B = -b*b.';

  G = @(t,U) C + U*B*U;
  F = @(t,U) Atfull*U + U*Afull + G(t,U);

  JF1 = @(t,U) Atfull + U*B;
  JF2 = JF1;
  JFcell = {JF1,JF2};

  KK = kronsum(At,At);
  odefun = @(t,u) KK*u + Cvec + ...
           reshape(reshape(u,d,d)*(B*reshape(u,d,d)),d*d,1);
  odeJac = @(t,u) KK + ...
           kronsum(reshape(u,d,d)*B,reshape(u,d,d).'*B.'); % it is sparse

  if nhat == 20
    load('lqc_steadystate_20.mat','Uss')
  end
end
