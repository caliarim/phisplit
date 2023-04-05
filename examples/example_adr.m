% Example of advection--diffusion--reaction (see [CC23, section 4.2])
%
% Equation:
% \partial_t u(t,x) = \epsilon\Delta u(t,x)
%                     +\alpha(\sum_\mu \partial_{x_\mu})
%                     + 1/(1+u(t,x)^2) + \Psi(t,x),
% u_0(x) = 64\sum_\mu x_\mu(1-x_\mu).
%
% with homogeneous Dirichlet boundary conditions.
%
% Time integration method: exponential Lawson, exponential Euler,
%                          Lawson2b and ETD2RK
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


disp('### Example ADR ###')

nspan = [40,80];
d = 3;
epsilon = 0.75;
alpha = 0.1*ones(1,d);
tstar = 1;

countersim = 0;
for n = nspan
  fprintf('-- Simulation with n=%i --\n',n)
  countersim = countersim + 1;
  nvec = n*ones(1,d) + (0:2);
  [U0,Acell,g,F,odefun,odeJac]=initialize_adr(nvec,epsilon,alpha);
  Uref = U0*exp(tstar);
  normref = norm(Uref(:),inf);

  disp('Exponential Euler phisplit')
  counter = 0;
  mspan_ee_phisplit = 50:400:1650;
  for m = mspan_ee_phisplit
    fprintf('Simulation with m=%i\n',m)
    counter = counter + 1;
    tic
    Uee_phisplit=expeuler_phisplit(U0,m,tstar,Acell,F);
    cpu_ee_phisplit(countersim,counter) = toc;
    err_ee_phisplit(countersim,counter) = ...
    norm(Uref(:)-Uee_phisplit(:),inf)/normref;
  end

  disp('Exponential Euler phiks')
  counter = 0;
  mspan_ee_phiks = 50:400:1650;
  if n == 40
    kappa_ee = 2^7;
  elseif n == 80
    kappa_ee = 2^6;
  end
  for m = mspan_ee_phiks
    fprintf('Simulation with m=%i\n',m)
    counter = counter + 1;
    tic
    Uee_phiks=expeuler_phiks(U0,m,tstar,Acell,F,kappa_ee);
    cpu_ee_phiks(countersim,counter) = toc;
    err_ee_phiks(countersim,counter) = norm(Uref(:)-Uee_phiks(:),inf)/normref;
  end

  disp('Lawson Euler')
  counter = 0;
  mspan_le = 800:8000:32800;
  for m = mspan_le
    fprintf('Simulation with m=%i\n',m)
    counter = counter + 1;
    tic
    Ule=lawsoneuler(U0,m,tstar,Acell,g);
    cpu_le(countersim,counter) = toc;
    err_le(countersim,counter) = norm(Uref(:)-Ule(:),inf)/normref;
  end

  disp('ETD2RK phisplit')
  counter = 0;
  mspan_rk2_phisplit = 40:100:440;
  for m = mspan_rk2_phisplit
    fprintf('Simulation with m=%i\n',m)
    counter = counter + 1;
    tic
    Urk2_phisplit=exprk2_phisplit(U0,m,tstar,Acell,g);
    cpu_rk2_phisplit(countersim,counter) = toc;
    err_rk2_phisplit(countersim,counter) = ...
    norm(Uref(:)-Urk2_phisplit(:),inf)/normref;
  end

  disp('ETD2RK phiks')
  counter = 0;
  mspan_rk2_phiks = 20:60:260;
  if n == 40
    kappa_rk2 = 2^7;
  elseif n == 80
    kappa_rk2 = 2^5;
  end
  for m = mspan_rk2_phiks
    fprintf('Simulation with m=%i\n',m)
    counter = counter + 1;
    tic
    Urk2_phiks=exprk2_phiks(U0,m,tstar,Acell,g,kappa_rk2);
    cpu_rk2_phiks(countersim,counter) = toc;
    err_rk2_phiks(countersim,counter) = norm(Uref(:)-Urk2_phiks(:),inf)/normref;
  end

  disp('Lawson2b')
  counter = 0;
  if n == 40
    mspan_l2b = 1500:4000:17500;
  elseif n == 80
    mspan_l2b = 3000:1500:9000;
  end
  for m = mspan_l2b
    fprintf('Simulation with m=%i\n',m)
    counter = counter + 1;
    tic
    Ul2b=lawson2b(U0,m,tstar,Acell,g);
    cpu_l2b(countersim,counter) = toc;
    err_l2b(countersim,counter) = norm(Uref(:)-Ul2b(:),inf)/normref;
  end

  if n == 40
    disp('Matlab Solver ode23t')
    counter = 0;
    tolrange = [8e-3,4e-5,1e-5,5e-6];
    for jj = 1:length(tolrange)
    fprintf('Simulation number %i\n',jj)
      options.RelTol = tolrange(jj);
      options.AbsTol = tolrange(jj);
      options.Jacobian = odeJac;
      tic
      res_odesolver=solver_matlab(U0,tstar,odefun,options,'ode23t');
      cpu_odesolver(countersim,jj)=toc;
      err_odesolver(countersim,jj) = norm(Uref(:)-res_odesolver(:),inf)/normref;
    end
  end

  disp('Order Exponential Euler phisplit')
  ln=length(mspan_ee_phisplit);
  den=log(mspan_ee_phisplit(1:ln-1))-log(mspan_ee_phisplit(2:ln));
  disp(-(log(err_ee_phisplit(countersim,1:ln-1))-...
         log(err_ee_phisplit(countersim,2:ln)))./den)

  disp('Order Exponential Euler phiks')
  ln=length(mspan_ee_phiks);
  den=log(mspan_ee_phiks(1:ln-1))-log(mspan_ee_phiks(2:ln));
  disp(-(log(err_ee_phiks(countersim,1:ln-1))-...
         log(err_ee_phiks(countersim,2:ln)))./den)

  disp('Order Lawson Euler')
  ln=length(mspan_le);
  den=log(mspan_le(1:ln-1))-log(mspan_le(2:ln));
  disp(-(log(err_le(countersim,1:ln-1))-log(err_le(countersim,2:ln)))./den)

  disp('Order ETD2RK phisplit')
  ln=length(mspan_rk2_phisplit);
  den=log(mspan_rk2_phisplit(1:ln-1))-log(mspan_rk2_phisplit(2:ln));
  disp(-(log(err_rk2_phisplit(countersim,1:ln-1))-...
         log(err_rk2_phisplit(countersim,2:ln)))./den)

  disp('Order ETD2RK phiks')
  ln=length(mspan_rk2_phiks);
  den=log(mspan_rk2_phiks(1:ln-1))-log(mspan_rk2_phiks(2:ln));
  disp(-(log(err_rk2_phiks(countersim,1:ln-1))-...
         log(err_rk2_phiks(countersim,2:ln)))./den)

  disp('Order Lawson2b')
  ln=length(mspan_l2b);
  den=log(mspan_l2b(1:ln-1))-log(mspan_l2b(2:ln));
  disp(-(log(err_l2b(countersim,1:ln-1))-log(err_l2b(countersim,2:ln)))./den)

  figure
  semilogy(cpu_ee_phisplit(countersim,:),err_ee_phisplit(countersim,:),'^--c')
  hold on
  semilogy(cpu_ee_phiks(countersim,:),err_ee_phiks(countersim,:),'^-c')
  semilogy(cpu_le(countersim,:),err_le(countersim,:),'x--g')
  title(sprintf('Precision diagram order 1 n=%i',n))
  legend('Exponential Euler phisplit',...
         'Exponential Euler phiks',...
         'Lawson Euler')
  xlabel('CPU')
  ylabel('Error')
  drawnow

  figure
  semilogy(cpu_rk2_phisplit(countersim,:),err_rk2_phisplit(countersim,:),'d--b')
  hold on
  semilogy(cpu_rk2_phiks(countersim,:),err_rk2_phiks(countersim,:),'d-b')
  semilogy(cpu_l2b(countersim,:),err_l2b(countersim,:),'s--g')
  if n == 40
    semilogy(cpu_odesolver(countersim,:),err_odesolver(countersim,:),'*--m')
    legend('ETD2RK phisplit',...
           'ETD2RK phiks',...
           'Lawson2b','matlab solver ode23t')
  elseif n == 80
    legend('ETD2RK phisplit',...
           'ETD2RK phiks',...
           'Lawson2b')
  end
  title(sprintf('Precision diagram order 2 n=%i',n))
  xlabel('CPU')
  ylabel('Error')
  drawnow

end

rmpath('integrators')
rmpath('../extern/phiks')
rmpath('../extern/KronPACK/src')
rmpath('../')
