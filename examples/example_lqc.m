% Example of linear quadratic control (see [CC23, section 4.1])
%
% Equation:
% U'(t) = A^T U(t) + U(t) A + C + U(t) B U(t),
% U(0) = 0.
%
% Time integration method: exponential Rosenbrock--Euler and ETD2RK
%
% [CC23] M. Caliari and F. Cassini.
%        Direction splitting of phi-functions in exponential integrators
%        for d-dimensional problems in Kronecker form, Submitted 2023

clear all
close all

addpath('../')
addpath('../extern/KronPACK/src')
addpath('../extern/phiks')
addpath('integrators')

disp('### Example linear quadratic control ###')

%%% Steady state check %%%
disp('-- Convergence to the steady state --')

nhat = 20;
tstar = 0.25;
m = 200;

[U0,Atcell,G,F,JFcell,~,~,Uss]=initialize_lqc(nhat);

disp('Exponential Rosenbrock Euler phisplit')
[~,err_ss_rb_phisplit]=exprbeuler_lqc_phisplit(U0,m,tstar,F,JFcell,Uss);

disp('Exponential Rosenbrock Euler phiks')
[~,err_ss_rb_phiks]=exprbeuler_phiks(U0,m,tstar,F,JFcell,1e-14,Uss);

disp('Exponential Runge-Kutta 2 phisplit')
[~,err_ss_rk2_phisplit]=exprk2_phisplit(U0,m,tstar,Atcell,G,Uss);

disp('Exponential Runge-Kutta 2 phiks')
[~,err_ss_rk2_phiks]=exprk2_phiks(U0,m,tstar,Atcell,G,1e-14,Uss);

tspan = linspace(0,tstar,m+1);
figure;
semilogy(tspan,err_ss_rb_phisplit,'o--r',...
         tspan,err_ss_rb_phiks,'o-r',...
         tspan,err_ss_rk2_phisplit,'d--b',...
         tspan,err_ss_rk2_phiks,'d-b')
legend('Exponential Rosenbrock Euler phisplit',...
       'Exponential Rosenbrock Euler phiks',...
       'ETD2RK phisplit',...
       'ETD2RK phiks')
title('Convergence to the steady state')
xlabel('time')
ylabel('Relative error Frobenius norm')
drawnow

%%% Simulations %%%
tstar = 0.025;
nhatspan = [30, 40];
countersim = 0;
for nhat = nhatspan
  fprintf('-- Simulation with nhat=%i --\n',nhat)
  countersim = countersim + 1;
  [U0,Atcell,G,F,JFcell,odefun]=initialize_lqc(nhat);

  load(sprintf('lqc_ref_%i.mat',nhat))
  normref = norm(Uref,'fro');

  disp('Matlab Solver ode23')
  if nhat == 30
    tolspan = [5e-3, 1e-3, 1e-4, 5.5e-5, 5e-5];
  elseif nhat == 40
    tolspan = [3e-3, 4e-4, 2.3e-4, 9.5e-5, 8.7e-5];
  end
  for i = 1:length(tolspan)
    fprintf('Simulation number %i\n',i)
    options.RelTol = tolspan(i);
    options.AbsTol = tolspan(i);
    tic
    res_odesolver=solver_matlab(U0,tstar,odefun,options,'ode23');
    cpu_odesolver(countersim,i)=toc;
    err_odesolver(countersim,i) = ...
    norm(Uref-reshape(res_odesolver,size(Uref)),'fro')/normref;
  end

  disp('Exponential Rosenbrock Euler phisplit')
  if nhat == 30
    mspanrb_phisplit = 30:35:170;
  elseif nhat == 40
    mspanrb_phisplit = 30:35:170;
  end
  counter = 0;
  for m = mspanrb_phisplit
    fprintf('Simulation with m=%i\n',m)
    counter = counter + 1;
    tic
    Urb_phisplit = exprbeuler_lqc_phisplit(U0,m,tstar,F,JFcell);
    cpu_rb_phisplit(countersim,counter) = toc;
    err_rb_phisplit(countersim,counter) = norm(Uref-Urb_phisplit,'fro')/normref;
  end

  disp('Exponential Rosenbrock Euler phiks')
  if nhat == 30
    mspanrb_phiks = 10:10:50;
    kappa_rb = 2^31;
  elseif nhat == 40
    mspanrb_phiks = 15:15:75;
    kappa_rb = 2^30;
  end
  counter = 0;
  for m = mspanrb_phiks
    fprintf('Simulation with m=%i\n',m)
    counter = counter + 1;
    tic
    Urb_phiks = exprbeuler_phiks(U0,m,tstar,F,JFcell,kappa_rb);
    cpu_rb_phiks(countersim,counter) = toc;
    err_rb_phiks(countersim,counter) = norm(Uref-Urb_phiks,'fro')/normref;
  end

  disp('Exponential Runge--Kutta 2 phisplit')
  if nhat == 30
    mspanrk2_phisplit = 30:35:170;
  elseif nhat == 40
    mspanrk2_phisplit = 30:35:170;
  end
  counter = 0;
  for m = mspanrk2_phisplit
    fprintf('Simulation with m=%i\n',m)
    counter = counter + 1;
    tic
    Urk2_phisplit = exprk2_phisplit(U0,m,tstar,Atcell,G);
    cpu_rk2_phisplit(countersim,counter) = toc;
    err_rk2_phisplit(countersim,counter) = ...
    norm(Uref-Urk2_phisplit,'fro')/normref;
  end

  disp('Exponential Runge-Kutta 2 phiks')
  if nhat == 30
    mspanrk2_phiks = 7:7:35;
    kappa_rk2 = 2^31;
  elseif nhat == 40
    mspanrk2_phiks = 10:10:50;
    kappa_rk2 = 2^27;
  end
  counter = 0;
  for m = mspanrk2_phiks
    fprintf('Simulation with m=%i\n',m)
    counter = counter + 1;
    tic
    Urk2_phiks = exprk2_phiks(U0,m,tstar,Atcell,G,kappa_rk2);
    cpu_rk2_phiks(countersim,counter) = toc;
    err_rk2_phiks(countersim,counter) = norm(Uref-Urk2_phiks,'fro')/normref;
  end

  disp('Order Exponential Rosenbrock Euler phisplit')
  ln=length(mspanrb_phisplit);
  den=log(mspanrb_phisplit(1:ln-1))-log(mspanrb_phisplit(2:ln));
  disp(-(log(err_rb_phisplit(countersim,1:ln-1))-...
         log(err_rb_phisplit(countersim,2:ln)))./den)
  disp('Order Exponential Rosenbrock Euler phiks')
  ln=length(mspanrb_phiks);
  den=log(mspanrb_phiks(1:ln-1))-log(mspanrb_phiks(2:ln));
  disp(-(log(err_rb_phiks(countersim,1:ln-1))-...
         log(err_rb_phiks(countersim,2:ln)))./den)
  disp('Order Exponential Runge-Kutta 2 phisplit')
  ln=length(mspanrk2_phisplit);
  den=log(mspanrk2_phisplit(1:ln-1))-log(mspanrk2_phisplit(2:ln));
  disp(-(log(err_rk2_phisplit(countersim,1:ln-1))-...
            log(err_rk2_phisplit(countersim,2:ln)))./den)
  disp('Order Exponential Runge-Kutta 2 phiks')
  ln=length(mspanrk2_phiks);
  den=log(mspanrk2_phiks(1:ln-1))-log(mspanrk2_phiks(2:ln));
  disp(-(log(err_rk2_phiks(countersim,1:ln-1))-...
         log(err_rk2_phiks(countersim,2:ln)))./den)

  figure
  semilogy(cpu_rk2_phisplit(countersim,:),err_rk2_phisplit(countersim,:),'d--b')
  hold on
  semilogy(cpu_rk2_phiks(countersim,:),err_rk2_phiks(countersim,:),'d-b')
  semilogy(cpu_rb_phisplit(countersim,:),err_rb_phisplit(countersim,:),'o--r')
  semilogy(cpu_rb_phiks(countersim,:),err_rb_phiks(countersim,:),'o-r')
  semilogy(cpu_odesolver(countersim,:),err_odesolver(countersim,:),'*--m')

  title(sprintf('Precision diagram order 2 nhat=%i',nhat))
  legend('Exponential Runge-Kutta 2 phisplit',...
         'Exponential Runge-Kutta 2 phiks',...
         'Exponential Rosenbrock Euler phisplit',...
         'Exponential Rosenbrock Euler phiks',...
         'matlab solver ode23')
  xlabel('CPU')
  ylabel('Error')
  drawnow

end

rmpath('integrators')
rmpath('../extern/phiks')
rmpath('../extern/KronPACK/src')
rmpath('../')
