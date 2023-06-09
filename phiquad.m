function varargout = phiquad (A, p, q, expm_fun)
%PHIQUAD phi-function of a matrix using a Gauss-Legendre-Lobatto
%   quadrature formula.
%   PHIQUAD(A,P) computes PHI_P(A), where PHI_P is the
%   exponential-like function defined by the integral
%
%   PHI_P(A) = \int_0^1 \exp((1-\theta)A)\theta^(P-1)/factorial(P-1)d\theta
%
%   When P == 0, the matrix exponential is returned.
%
%   PHIQUAD(A,P,Q) computes PHI_P(A) using Q quadrature nodes (the
%   default is nine). Q must be larger than 2 and smaller than 13.
%
%   [PHI_1, PHI_2,..., PHI_P] = PHIQUAD (A,P,___) compute all the
%   phi-functions up to order P.
%
%   [PHI_1, PHI_2,..., PHI_P, PHI_0] = PHIQUAD (A,P,___) compute all the
%   phi-functions up to order P and return the matrix exponential as
%   last output argument.
%
%   [___] = PHIQUAD (___, EXPM_FUN) employs the function EXPM_FUN to computer
%   the required matrix exponentials.
%
%   [CC23] M. Caliari and F. Cassini.
%          Direction splitting of phi-functions in exponential integrators
%          for d-dimensional problems in Kronecker form, Submitted 2023
  if ((nargout > 1 && nargout < p) || (nargout > p + 1))
    error ('Inconsistent number of output arguments, see documentation.')
  end
  if (nargin < 4)
    expm_fun = @(A) expm(A);
  end
  if (p == 0)
    varargout{1} = expm_fun (A);
    return
  end
  s = max(0, ceil (log2 (norm(A,1)))); % scaling
  A = A / 2 ^ s;
  if ((nargin == 2) || isempty (q))
      q = 9;
  end
  [x, w] = XW (q);
  fact = w(q);
  Ei = eye (size (A));
  varargout = cell(1,max(nargout,1));
  for l = 1:p % quadrature from theta=1
    varargout{l} = fact * Ei;
    fact = fact / l;
  end
  for i = q - 1:-1:ceil (q / 2)
    Eitemp{q - i} = expm_fun ((x(i + 1) - x(i)) * A);
    Ei = Eitemp{q - i} * Ei;
    fact = w(i);
    for l = 1:p
      varargout{l} = varargout{l} + fact * Ei;
      fact = fact * x(i) / l;
    end
  end
  for i = floor ((q - 1) / 2):-1:2
    Ei = Eitemp{i} * Ei;
    fact = w(i);
    for l = 1:p
      varargout{l} = varargout{l} + fact * Ei;
      fact = fact * x(i) / l;
    end
  end
  Ei = Eitemp{1} * Ei;
  varargout{1} = varargout{1} + w(1) * Ei; % theta=0 for PHI_1
  for j = 1:s % squaring
    for l = p:-1:1
      varargout{l} = varargout{l} + Ei * varargout{l};
      fact = factorial (l - 1);
      for k = 1:l - 1
        varargout{l} = varargout{l} + varargout{k} / fact;
        fact = fact / (l - k);
      end
      varargout{l} = varargout{l} / 2 ^ l;
    end
    Ei = Ei * Ei;
  end
  if (nargout <= 1)
    varargout{1} = varargout{p};
  end
  if (nargout == p + 1)
    varargout{p + 1} = Ei;
  end
end

function [x,w] = XW(q) % Gauss-Legendre-Lobatto quadrature
  X = [0,...
       0,...
       0,...
       0,...
       0,...
       0,...
       0,...
       0,...
       0,...
       0;...
       5e-01,...
       2.763932022500211e-01,...
       1.726731646460114e-01,...
       1.174723380352676e-01,...
       8.488805186071652e-02,...
       6.412992574519671e-02,...
       5.012100229426991e-02,...
       4.023304591677057e-02,...
       3.299928479597042e-02,...
       2.755036388855892e-02;...
       1,...
       7.236067977499789e-01,...
       5e-01,...
       3.573842417596774e-01,...
       2.655756032646429e-01,...
       2.041499092834289e-01,...
       1.614068602446311e-01,...
       1.306130674472475e-01,...
       1.077582631684278e-01,...
       9.036033917799668e-02;...
       NaN,...
       1,...
       8.273268353539887e-01,...
       6.426157582403226e-01,...
       5e-01,...
       3.953503910487606e-01,...
       3.184412680869109e-01,...
       2.610375250947777e-01,...
       2.173823365018975e-01,...
       1.835619234840696e-01;...
       NaN,...
       NaN,...
       1,...
       8.825276619647324e-01,...
       7.344243967353571e-01,...
       6.046496089512394e-01,...
       5e-01,...
       4.173605211668065e-01,...
       3.521209322065303e-01,...
       3.002345295173255e-01;...
       NaN,...
       NaN,...
       NaN,...
       1,...
       9.151119481392835e-01,...
       7.958500907165711e-01,...
       6.815587319130891e-01,...
       5.826394788331936e-01,...
       5e-01,...
       4.317235335725362e-01;...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       1,...
       9.358700742548033e-01,...
       8.385931397553689e-01,...
       7.389624749052223e-01,...
       6.478790677934697e-01,...
       5.682764664274638e-01;...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       1,...
       9.498789977057300e-01,...
       8.693869325527526e-01,...
       7.826176634981026e-01,...
       6.997654704826745e-01;...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       1,...
       9.597669540832294e-01,...
       8.922417368315723e-01,...
       8.164380765159304e-01;...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       1,...
       9.670007152040296e-01,...
       9.096396608220033e-01;...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       1,...
       9.724496361114411e-01;...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       1];
  W = [1.666666666666667e-01,...
       6.666666666666666e-01,...
       1.666666666666667e-01,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN;...
       8.333333333333333e-02,...
       4.166666666666667e-01,...
       4.166666666666667e-01,...
       8.333333333333333e-02,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN;...
       5e-02,...
       2.722222222222222e-01,...
       3.555555555555556e-01,...
       2.722222222222222e-01,...
       5e-02,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN;...
       3.333333333333333e-02,...
       1.892374781489235e-01,...
       2.774291885177432e-01,...
       2.774291885177432e-01,...
       1.892374781489235e-01,...
       3.333333333333333e-02,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN;...
       2.380952380952381e-02,...
       1.384130236807830e-01,...
       2.158726906049313e-01,...
       2.438095238095238e-01,...
       2.158726906049313e-01,...
       1.384130236807830e-01,...
       2.380952380952381e-02,...
       NaN,...
       NaN,...
       NaN,...
       NaN,...
       NaN;...
       1.785714285714286e-02,...
       1.053521135717530e-01,...
       1.705613462417522e-01,...
       2.062293973293519e-01,...
       2.062293973293519e-01,...
       1.705613462417522e-01,...
       1.053521135717530e-01,...
       1.785714285714286e-02,...
       NaN,...
       NaN,...
       NaN,...
       NaN;...
       1.388888888888889e-02,...
       8.274768078040276e-02,...
       1.372693562500809e-01,...
       1.732142554865232e-01,...
       1.857596371882086e-01,...
       1.732142554865232e-01,...
       1.372693562500809e-01,...
       8.274768078040276e-02,...
       1.388888888888889e-02,...
       NaN,...
       NaN,...
       NaN;...
       1.111111111111111e-02,...
       6.665299542553506e-02,...
       1.124446710315632e-01,...
       1.460213418398419e-01,...
       1.637698805919487e-01,...
       1.637698805919487e-01,...
       1.460213418398419e-01,...
       1.124446710315632e-01,...
       6.665299542553506e-02,...
       1.111111111111111e-02,...
       NaN,...
       NaN;...
       9.090909090909090e-03,...
       5.480613663349743e-02,...
       9.358494089015260e-02,...
       1.240240521320142e-01,...
       1.434395623895040e-01,...
       1.501087977278454e-01,...
       1.434395623895040e-01,...
       1.240240521320142e-01,...
       9.358494089015260e-02,...
       5.480613663349743e-02,...
       9.090909090909090e-03,...
       NaN;...
       7.575757575757576e-03,...
       4.584225870659807e-02,...
       7.898735278218506e-02,...
       1.062542088805106e-01,...
       1.256378015996006e-01,...
       1.357026204553481e-01,...
       1.357026204553481e-01,...
       1.256378015996006e-01,...
       1.062542088805106e-01,...
       7.898735278218506e-02,...
       4.584225870659807e-02,...
       7.575757575757576e-03];
  x = X(1:q,q-2);
  w = W(q-2,1:q);
end
%!test
%! A = toeplitz([-2,1,0,0]);
%! assert(phiquad(A,0),expm(A),1e-10)
%!test
%! A = toeplitz([-2,1,0,0]);
%! res = phiquad(A,0);
%! assert(res,expm(A),1e-10)
%!test
%! A = toeplitz([-2,1,0,0]);
%! E = expm(A);
%! I = eye(size(A));
%! assert(phiquad(A,1),(E-I)/A,1e-10)
%!test
%! A = toeplitz([-2,1,0,0]);
%! E = expm(A);
%! I = eye(size(A));
%! assert(phiquad(A,2),(E-I-A)/A/A,1e-10)
%!test
%! A = toeplitz([-2,1,0,0]);
%! E = expm(A);
%! I = eye(size(A));
%! res1 = phiquad(A,1);
%! assert(res1,(E-I)/A,1e-10)
%!test
%! A = toeplitz([-2,1,0,0]);
%! E = expm(A);
%! I = eye(size(A));
%! [res1,res2] = phiquad(A,2);
%! assert(res1,(E-I)/A,1e-10)
%! assert(res2,(E-I-A)/A/A,1e-10)
%!test
%! A = toeplitz([-2,1,0,0]);
%! E = expm(A);
%! I = eye(size(A));
%! [res1,res2] = phiquad(A,1);
%! assert(res1,(E-I)/A,1e-10)
%! assert(res2,expm(A),1e-10)
%!test
%! A = toeplitz([-2,1,0,0]);
%! E = expm(A);
%! I = eye(size(A));
%! [res1,res2,res3] = phiquad(A,3);
%! assert(res1,(E-I)/A,1e-10)
%! assert(res2,(E-I-A)/A/A,1e-10)
%! assert(res3,(E-I-A-A*A/2)/A/A/A,1e-10)
%!error
%! A = toeplitz([-2,1,0,0]);
%! [res1,res2] = phiquad(A,3);
%!test
%! A = toeplitz([-2,1,0,0]);
%! E = expm(A);
%! I = eye(size(A));
%! [res1,res2,res3,res4] = phiquad(A,4);
%! assert(res1,(E-I)/A,1e-10)
%! assert(res2,(E-I-A)/A/A,1e-10)
%! assert(res3,(E-I-A-A*A/2)/A/A/A,1e-10)
%! assert(res4,(E-I-A-A*A/2-A*A*A/6)/A/A/A/A,1e-10)
%!test
%! A = toeplitz([-2,1,0,0]);
%! E = expm(A);
%! I = eye(size(A));
%! [res1,res2,res3,res4] = phiquad(A,4,8);
%! assert(res1,(E-I)/A,1e-10)
%! assert(res2,(E-I-A)/A/A,1e-10)
%! assert(res3,(E-I-A-A*A/2)/A/A/A,1e-10)
%! assert(res4,(E-I-A-A*A/2-A*A*A/6)/A/A/A/A,1e-10)
%!test
%! A = toeplitz([-2,1,0,0]);
%! E = expm(A);
%! I = eye(size(A));
%! [res1,res2,res3,res4] = phiquad(A,4,7);
%! assert(res1,(E-I)/A,1e-10)
%! assert(res2,(E-I-A)/A/A,1e-10)
%! assert(res3,(E-I-A-A*A/2)/A/A/A,1e-10)
%! assert(res4,(E-I-A-A*A/2-A*A*A/6)/A/A/A/A,1e-10)
%!test
%! A = toeplitz([-2,1,0,0]);
%! E = expm(A);
%! I = eye(size(A));
%! [res1,res2,res3,res4] = phiquad(A,4,3);
%! assert(res1,(E-I)/A,1e-2)
%! assert(res2,(E-I-A)/A/A,1e-2)
%! assert(res3,(E-I-A-A*A/2)/A/A/A,1e-2)
%! assert(res4,(E-I-A-A*A/2-A*A*A/6)/A/A/A/A,1e-2)
