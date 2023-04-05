function [res_ode,err_ss_ode]=solver_matlab(U0,tstar,odefun,options,method)
% SOLVER_MATLAB Matlab solver from ODE suite.
%
% See also EXAMPLE_LQC, EXAMPLE_ADR.

  global res_ode;
  options.OutputFcn = @(t,u,flag) myoutfcn(t,u,flag,tstar);
  feval(method,odefun,[0,tstar],U0(:),options);
end

function status = myoutfcn(t,app,flag,tstar)
% MYOUTFCN Output function for matlab solver from ODE suite.
%
% See also SOLVER_MATLAB.

  switch flag
    case 'init'
      status = 0;
    case ''
      if (t(end) == tstar)
        global res_ode;
        res_ode = app(:,end);
      end
      status = 0;
    case 'done'
      status = 1;
  end
end
