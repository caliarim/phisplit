# PHISPLIT -- phi-functions of Kronecker sums with direction splitting

The **PHISPLIT** algorithm approximates actions of phi-functions for
d-dimensional Kronecker sums of matrices using a mu-mode approach.
The technique is based on a direction splitting of the involved matrix
functions, which generates an approximation error compatible with
an exponential integrator up to second order.
For more details, see the reference manuscript

[M. Caliari and F. Cassini. Direction splitting of phi-functions in exponential
integrators for d-dimensional problems in Kronecker form,
arXiv preprint arXiv:2304.02327, 2023](https://arxiv.org/abs/2304.02327)

This GitHub repository contains:
- ```phisplit.m``` (fully compatible with GNU Octave) which approximates
action of phi-functions of Kronecker sums on a tensor
- ```phiquad.m``` (fully compatible with GNU Octave) whith computes
matrix phi-functions using a quadrature rule
-  all the functions and the scripts needed to reproduce the numerical
examples contained in the reference manuscript (in the folder ```examples```)

## Software

You can directly download the zip archive of the latest complete release
[here](https://github.com/caliarim/phisplit/releases/download/v0.1/phisplit-0.1.zip).

Alternatively, since the repository makes use of git submodules, you can
obtain all the relevant sources by executing

```sh
git clone --recursive https://github.com/caliarim/phisplit
```

## Installing

The software is all consituted by MATLAB scripts.
In order use **PHISPLIT** it is needed to have in path the function
```tucker.m``` from the package [KronPACK](https://github.com/caliarim/KronPACK), already
available in the repository (see the section **Software**).
In MATLAB/GNU Octave, this can be achieved by running the command

```sh
addpath('extern/KronPACK/src')
```

inside the folder ```phisplit```.

## Testing in GNU Octave

**PHISPLIT** contains several GNU Octave built-in self-tests that can be
executed by running the command

```sh
test phisplit
```

For each test, the first argument in the GNU Octave assert command is
compared to the second argument, up to the tolerance possibly given
in the third argument. Therefore, it is possible to perform the same
tests in MATLAB by slightly modifying the syntax of the test.
